
%Code modified from LAP_5 code created on 09/28 - this just has a different definition and values -
%looking at the sequences not transitions - reaction time is stimulus onset - speech offset

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Define the consonant-vowel sequences and their phonotactic probabilities
phonotacticProbabilities = struct(...
    'ghoo', 0.0004, 'voo', 0, 'soo', 0.0023, 'too', 0.0021, ...
    'ghah', 0.002, 'vah', 0.001, 'sah', 0.0022, 'tah', 0.0022, ...
    'ghee', 0.0001, 'vee', 0.001, 'see', 0.0028, 'tee', 0.0012);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);
subtable.triplet_tables = cell(nsubs_to_analyze, 1);

% Arrays to store significant results
significant_seq_pp_rt_dbsIDs = {};

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));

    % Load phoneme and triplet tables if they exist
    if exist(phonemeFilePath, 'file')
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

        % Add empty columns for RT and phonotactic probabilities for each sequence
        subtable.phoneme_tables{i_sub}.seq1_RT = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.seq2_RT = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.seq3_RT = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.pp1 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.pp2 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.pp3 = NaN(height(subtable.phoneme_tables{i_sub}), 1);

        % Get unique trial_ids
        unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);

        % Check if 'onset' exists, else use 'starts'
        if ismember('onset', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
            timing_field = 'onset';
        elseif ismember('starts', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
            timing_field = 'starts';
        else
            error('Neither onset nor starts column found.');
        end

        % Initialize arrays for correlation analysis for the current subject
        RT_values = [];
        PP_values = [];

        % Loop through each trial
        for trial = 1:length(unique_trials)
            trial_id = unique_trials(trial);

            % Get the phonemes for this trial
            trial_phonemes = subtable.phoneme_tables{i_sub}(subtable.phoneme_tables{i_sub}.trial_id == trial_id, :);

            % Ensure there are enough phonemes to form the 3 sequences (seq1, seq2, seq3)
            if height(trial_phonemes) >= 6
                % Calculate reaction times for each sequence
                seq1_RT = trial_phonemes.(timing_field)(2) - trial_phonemes.(timing_field)(1); % first consonant-vowel sequence
                seq2_RT = trial_phonemes.(timing_field)(4) - trial_phonemes.(timing_field)(3); % second consonant-vowel sequence
                seq3_RT = trial_phonemes.(timing_field)(6) - trial_phonemes.(timing_field)(5); % third consonant-vowel sequence

                % Skip trials with NaN values
                if isnan(seq1_RT) || isnan(seq2_RT) || isnan(seq3_RT)
                    continue;
                end

                % Get the sequences (consonant + vowel)
                seq1 = strcat(trial_phonemes.stim{1}, trial_phonemes.stim{2});
                seq2 = strcat(trial_phonemes.stim{3}, trial_phonemes.stim{4});
                seq3 = strcat(trial_phonemes.stim{5}, trial_phonemes.stim{6});

                % Get phonotactic probabilities for each sequence
                if isfield(phonotacticProbabilities, seq1)
                    pp1 = phonotacticProbabilities.(seq1);
                else
                    pp1 = NaN;
                end

                if isfield(phonotacticProbabilities, seq2)
                    pp2 = phonotacticProbabilities.(seq2);
                else
                    pp2 = NaN;
                end

                if isfield(phonotacticProbabilities, seq3)
                    pp3 = phonotacticProbabilities.(seq3);
                else
                    pp3 = NaN;
                end

                % Store reaction times and phonotactic probabilities
                row_indices = find(subtable.phoneme_tables{i_sub}.trial_id == trial_id);
                subtable.phoneme_tables{i_sub}.seq1_RT(row_indices(1)) = seq1_RT;
                subtable.phoneme_tables{i_sub}.seq2_RT(row_indices(1)) = seq2_RT;
                subtable.phoneme_tables{i_sub}.seq3_RT(row_indices(1)) = seq3_RT;
                subtable.phoneme_tables{i_sub}.pp1(row_indices(1)) = pp1;
                subtable.phoneme_tables{i_sub}.pp2(row_indices(1)) = pp2;
                subtable.phoneme_tables{i_sub}.pp3(row_indices(1)) = pp3;

                % Collect data for correlation analysis
                RT_values = [RT_values; seq1_RT; seq2_RT; seq3_RT];
                PP_values = [PP_values; pp1; pp2; pp3];
            end
        end

        % Perform correlation analysis for the current subject
        if ~isempty(RT_values) && ~isempty(PP_values)
            [corr_value, p_value] = corr(PP_values, RT_values, 'Rows', 'complete');
            disp(['Correlation for ', dbsID, ' between phonotactic probability and reaction time: ', num2str(corr_value), ' (p-value: ', num2str(p_value), ')']);

            % If significant correlation, store DBS ID
            if p_value < 0.05
                significant_seq_pp_rt_dbsIDs = [significant_seq_pp_rt_dbsIDs; {dbsID}];
            end
        end

        % Write the updated phoneme table to a file
        writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
    end
end

disp('Data processing and correlation analysis complete.');

% Output the significant DBS IDs
disp('Significant correlations found in:');
disp(significant_seq_pp_rt_dbsIDs);
%% Calculating means and STD 

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Define the consonant-vowel sequences and their phonotactic probabilities
phonotacticProbabilities = struct(...
    'ghoo', 0.0004, 'voo', 0, 'soo', 0.0023, 'too', 0.0021, ...
    'ghah', 0.002, 'vah', 0.001, 'sah', 0.0022, 'tah', 0.0022, ...
    'ghee', 0.0001, 'vee', 0.001, 'see', 0.0028, 'tee', 0.0012);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);
subtable.triplet_tables = cell(nsubs_to_analyze, 1);

% Initialize containers to collect all reaction times per sequence
sequences = fieldnames(phonotacticProbabilities);  % Get all sequence names
all_RT = struct();  % To store all reaction times for each sequence

% Initialize the containers for mean and std
meanRT = zeros(length(sequences), 1);
stdRT = zeros(length(sequences), 1);

for i = 1:length(sequences)
    all_RT.(sequences{i}) = [];  % Initialize empty arrays for each sequence
end

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));

    % Load phoneme and triplet tables if they exist
    if exist(phonemeFilePath, 'file')
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

        % Get unique trial_ids
        unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);

        % Check if 'onset' exists, else use 'starts'
        if ismember('onset', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
            timing_field = 'onset';
        elseif ismember('starts', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
            timing_field = 'starts';
        else
            error('Neither onset nor starts column found.');
        end

        % Loop through each trial
        for trial = 1:length(unique_trials)
            trial_id = unique_trials(trial);

            % Get the phonemes for this trial
            trial_phonemes = subtable.phoneme_tables{i_sub}(subtable.phoneme_tables{i_sub}.trial_id == trial_id, :);

            % Ensure there are enough phonemes to form the 3 sequences (seq1, seq2, seq3)
            if height(trial_phonemes) >= 6
                % Calculate reaction times for each sequence
                seq1_RT = trial_phonemes.(timing_field)(2) - trial_phonemes.(timing_field)(1); % first consonant-vowel sequence
                seq2_RT = trial_phonemes.(timing_field)(4) - trial_phonemes.(timing_field)(3); % second consonant-vowel sequence
                seq3_RT = trial_phonemes.(timing_field)(6) - trial_phonemes.(timing_field)(5); % third consonant-vowel sequence

                % Skip trials with NaN values
                if isnan(seq1_RT) || isnan(seq2_RT) || isnan(seq3_RT)
                    continue;
                end

                % Get the sequences (consonant + vowel)
                seq1 = strcat(trial_phonemes.stim{1}, trial_phonemes.stim{2});
                seq2 = strcat(trial_phonemes.stim{3}, trial_phonemes.stim{4});
                seq3 = strcat(trial_phonemes.stim{5}, trial_phonemes.stim{6});

                % Append the reaction times to the corresponding sequence array
                if isfield(all_RT, seq1)
                    all_RT.(seq1) = [all_RT.(seq1); seq1_RT];
                end
                if isfield(all_RT, seq2)
                    all_RT.(seq2) = [all_RT.(seq2); seq2_RT];
                end
                if isfield(all_RT, seq3)
                    all_RT.(seq3) = [all_RT.(seq3); seq3_RT];
                end
            end
        end
    end
end

% Now calculate the mean and std of reaction times for each sequence
for i = 1:length(sequences)
    if ~isempty(all_RT.(sequences{i}))  % Ensure there are data points
        meanRT(i) = mean(all_RT.(sequences{i}));
        stdRT(i) = std(all_RT.(sequences{i}));
    end
end

% Sort the phonotactic probabilities and corresponding data
[sortedPhonotacticProb, sortIdx] = sort(phonotacticProb);  % Sort phonotactic probabilities in ascending order
sortedMeanRT = meanRT(sortIdx);  % Sort mean reaction times based on the sorting order
sortedStdRT = stdRT(sortIdx);    % Sort standard deviations similarly

% Create the scatter plot with error bars
figure;
errorbar(sortedPhonotacticProb, sortedMeanRT, sortedStdRT, 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time');
grid on;

% Add a trend line (linear fit)
hold on;
p = polyfit(sortedPhonotacticProb, sortedMeanRT, 1); % Linear fit
yfit = polyval(p, sortedPhonotacticProb);
plot(sortedPhonotacticProb, yfit, '-r', 'LineWidth', 2); % Plot the fitted line
hold off;

% Customize the plot for better readability
xticks(sortedPhonotacticProb); % Set x-ticks at the sorted phonotactic probability values
xtickformat('%.4f'); % Format x-tick labels for better readability
xtickangle(45); % Rotate the x-tick labels for better visibility

%% 
%Trying Better Plotting
% Define colors for each sequence
colors = lines(length(sequences));  % Generate distinct colors for each sequence

% Create the scatter plot with error bars
figure;
hold on;
for i = 1:length(sequences)
    errorbar(sortedPhonotacticProb(i), sortedMeanRT(i), sortedStdRT(i), 'o', 'MarkerSize', 8, 'LineWidth', 1.5, ...
             'CapSize', 10, 'MarkerEdgeColor', colors(i, :), 'MarkerFaceColor', colors(i, :));  % Color each point
end

% Add labels and grid
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time');
grid on;

% Add a trend line (linear fit)
p = polyfit(sortedPhonotacticProb, sortedMeanRT, 1); % Linear fit
yfit = polyval(p, sortedPhonotacticProb);
plot(sortedPhonotacticProb, yfit, '-r', 'LineWidth', 2); % Plot the fitted line

% Add a legend to map colors to sequences
legend(sequences, 'Location', 'bestoutside');  % Add a legend showing sequence names

hold off;

%%
% Create the scatter plot with error bars
figure;
hold on;
for i = 1:length(sequences)
    errorbar(sortedPhonotacticProb(i), sortedMeanRT(i), sortedStdRT(i), 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);
    
    % Add a label for each sequence next to its point
    text(sortedPhonotacticProb(i), sortedMeanRT(i), sequences{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

% Add labels and grid
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time)');
grid on;

% Add a trend line (linear fit)
p = polyfit(sortedPhonotacticProb, sortedMeanRT, 1); % Linear fit
yfit = polyval(p, sortedPhonotacticProb);
plot(sortedPhonotacticProb, yfit, '-r', 'LineWidth', 2); % Plot the fitted line

hold off;
%% Overlap between phonotactic probabilities and VC sequenc - plotting based on initial consonant 
% Define the sequences starting with 's'
s_sequences = {'soo', 'see', 'sah'};

% Create the figure for the 's' sequences
figure;
hold on;

% Loop over each 's' sequence to plot phonotactic probability vs. mean reaction time
for seqIdx = 1:length(s_sequences)
    sequenceName = s_sequences{seqIdx};
    
    % Find the index of the sequence in the sorted list
    seqIndex = find(strcmp(sequences, sequenceName));
    
    % Plot the data for the current sequence with error bars
    errorbar(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sortedStdRT(seqIndex), 'o', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);  
    
    % Label the data point with the sequence name
    text(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sequenceName, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

% Add labels and title
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time for "s" Sequences');
grid on;

hold off;

%V
% Define the sequences starting with 'v'
v_sequences = {'voo', 'vee', 'vah'};

% Create the figure for the 'v' sequences
figure;
hold on;

% Loop over each 'v' sequence to plot phonotactic probability vs. mean reaction time
for seqIdx = 1:length(v_sequences)
    sequenceName = v_sequences{seqIdx};
    
    % Find the index of the sequence in the sorted list
    seqIndex = find(strcmp(sequences, sequenceName));
    
    % Plot the data for the current sequence with error bars
    errorbar(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sortedStdRT(seqIndex), 'o', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);  
    
    % Label the data point with the sequence name
    text(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sequenceName, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

% Add labels and title
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time for "v" Sequences');
grid on;

hold off;

%G
% Define the sequences starting with 'g'
g_sequences = {'ghoo', 'ghee', 'ghah'};

% Create the figure for the 'g' sequences
figure;
hold on;

% Loop over each 'g' sequence to plot phonotactic probability vs. mean reaction time
for seqIdx = 1:length(g_sequences)
    sequenceName = g_sequences{seqIdx};
    
    % Find the index of the sequence in the sorted list
    seqIndex = find(strcmp(sequences, sequenceName));
    
    % Plot the data for the current sequence with error bars
    errorbar(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sortedStdRT(seqIndex), 'o', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);  
    
    % Label the data point with the sequence name
    text(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sequenceName, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

% Add labels and title
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time for "g" Sequences');
grid on;

hold off;

%%
%T
% Define the sequences starting with 't'
t_sequences = {'too', 'tee', 'tah'};

% Create the figure for the 't' sequences
figure;
hold on;

% Loop over each 't' sequence to plot phonotactic probability vs. mean reaction time
for seqIdx = 1:length(t_sequences)
    sequenceName = t_sequences{seqIdx};
    
    % Find the index of the sequence in the sorted list
    seqIndex = find(strcmp(sequences, sequenceName));
    
    % Plot the data for the current sequence with error bars
    errorbar(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sortedStdRT(seqIndex), 'o', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 10);  
    
    % Label the data point with the sequence name
    text(sortedPhonotacticProb(seqIndex), sortedMeanRT(seqIndex), sequenceName, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

% Add labels and title
xlabel('Phonotactic Probability');
ylabel('Mean Reaction Time (ms)');
title('Phonotactic Probability vs. Reaction Time for "t" Sequences');
grid on;

hold off;
%%


