%This code correctly assigns trans1,trans2,rt1,rt2,pp1,pp2 - then performs corr analysis - finally


clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Define the phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'oot', 0, ...
    'eeg', 0.0005, ...
    'oog', 0, ...
    'ahs', 0.0007, ...
    'aht', 0.0001, ...
    'ahv', 0.0001, ...
    'ees', 0.0003, ...
    'eev', 0.0007);

% List of DBS IDs
% % % % % % % % % % % % % % % % % % % % % % % % % dbsIDs = {'DBS3002', 'DBS3004', 'DBS3006', 'DBS3010', 'DBS3011', 'DBS3014', ...
% % % % % % % % % % % % % % % % % % % % % % % % %           'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
% % % % % % % % % % % % % % % % % % % % % % % % %           'DBS3032', 'DBS4079', 'DBS4080'};

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.triplet_tables = cell(nsubs_to_analyze,1);

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));
   
    % Load phoneme and triplet tables if they exist
    if exist(phonemeFilePath, 'file')
        % Load phoneme table
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

        % Add empty columns for trans1, trans2, rt1, rt2, pp1, and pp2
        subtable.phoneme_tables{i_sub}.trans1 = cell(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.trans2 = cell(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.rt1 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.rt2 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.pp1 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.pp2 = NaN(height(subtable.phoneme_tables{i_sub}), 1);

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
        rt1_values = [];
        pp1_values = [];
        rt2_values = [];
        pp2_values = [];

        % Loop through each trial
        for trial = 1:length(unique_trials)
            trial_id = unique_trials(trial);
            
            % Get the phonemes for this trial
            trial_phonemes = subtable.phoneme_tables{i_sub}(subtable.phoneme_tables{i_sub}.trial_id == trial_id, :);
            
            % Ensure there are enough phonemes to form trans1 and trans2
            if height(trial_phonemes) >= 4
                % trans1: first vowel + first consonant
                trans1 = strcat(trial_phonemes.stim{2}, trial_phonemes.stim{3});
                % trans2: second vowel + second consonant
                trans2 = strcat(trial_phonemes.stim{4}, trial_phonemes.stim{5});

                % Get the row indices for the current trial
                row_indices = find(subtable.phoneme_tables{i_sub}.trial_id == trial_id);

                % Calculate rt1: time difference between consonant 1 and vowel 1
                rt1 = trial_phonemes.(timing_field)(3) - trial_phonemes.(timing_field)(2);  % consonant - vowel

                % Calculate rt2: time difference between consonant 2 and vowel 2
                rt2 = trial_phonemes.(timing_field)(5) - trial_phonemes.(timing_field)(4);  % consonant - vowel

                % Phonotactic probabilities for trans1 and trans2
                if isfield(phonotacticProbabilities, trans1)
                    pp1 = phonotacticProbabilities.(trans1);
                else
                    pp1 = NaN; % Default to NaN if not found
                end

                if isfield(phonotacticProbabilities, trans2)
                    pp2 = phonotacticProbabilities.(trans2);
                else
                    pp2 = NaN; % Default to NaN if not found
                end

                % Assign trans1, trans2, rt1, rt2, pp1, and pp2 only to the first row of the trial
                subtable.phoneme_tables{i_sub}.trans1(row_indices(1)) = {trans1};
                subtable.phoneme_tables{i_sub}.trans2(row_indices(1)) = {trans2};
                subtable.phoneme_tables{i_sub}.rt1(row_indices(1)) = rt1;
                subtable.phoneme_tables{i_sub}.rt2(row_indices(1)) = rt2;
                subtable.phoneme_tables{i_sub}.pp1(row_indices(1)) = pp1;
                subtable.phoneme_tables{i_sub}.pp2(row_indices(1)) = pp2;

                % Collect data for correlation analysis (where data is available)
                if ~isnan(rt1) && ~isnan(pp1)
                    rt1_values = [rt1_values; rt1];
                    pp1_values = [pp1_values; pp1];
                end
                if ~isnan(rt2) && ~isnan(pp2)
                    rt2_values = [rt2_values; rt2];
                    pp2_values = [pp2_values; pp2];
                end

                % Set trans1, trans2, rt1, rt2, pp1, and pp2 to NaN for all other rows in the same trial
                subtable.phoneme_tables{i_sub}.trans1(row_indices(2:end)) = {NaN};
                subtable.phoneme_tables{i_sub}.trans2(row_indices(2:end)) = {NaN};
                subtable.phoneme_tables{i_sub}.rt1(row_indices(2:end)) = NaN;
                subtable.phoneme_tables{i_sub}.rt2(row_indices(2:end)) = NaN;
                subtable.phoneme_tables{i_sub}.pp1(row_indices(2:end)) = NaN;
                subtable.phoneme_tables{i_sub}.pp2(row_indices(2:end)) = NaN;
            end
        end

        % Perform correlation analysis between rt1 and pp1, rt2 and pp2 for the current subject
        if ~isempty(rt1_values) && ~isempty(pp1_values)
            [corr_rt1_pp1, pval_rt1_pp1] = corr(rt1_values, pp1_values, 'Rows', 'complete');
            disp(['Correlation for ', dbsID, ' between rt1 and pp1: ', num2str(corr_rt1_pp1), ' (p-value: ', num2str(pval_rt1_pp1), ')']);
        end
        if ~isempty(rt2_values) && ~isempty(pp2_values)
            [corr_rt2_pp2, pval_rt2_pp2] = corr(rt2_values, pp2_values, 'Rows', 'complete');
            disp(['Correlation for ', dbsID, ' between rt2 and pp2: ', num2str(corr_rt2_pp2), ' (p-value: ', num2str(pval_rt2_pp2), ')']);
        end

        % Write the updated phoneme table with trans1, trans2, rt1, rt2, pp1, and pp2 to a file
        writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
        writetable(subtable.triplet_tables{i_sub}, fullfile(outputDir, sprintf('%s_triplet_table_processed.txt', dbsID)));
    end
end

disp('Data processing and correlation analysis complete.');

%%
% Initialize variables to store mean rt1 and rt2 values per DBS ID
mean_rt1_per_subject = [];
mean_rt2_per_subject = [];
subject_ids = {};

% Initialize counters
rt1_greater_count = 0;
rt2_greater_count = 0;

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    % Collect the current subject's ID
    dbsID = subtable.subject{i_sub};
    
    % Ensure that subtable.phoneme_tables{i_sub} is a table
    phoneme_table = subtable.phoneme_tables{i_sub};
    
    % Check if phoneme_table is a table and contains rt1 and rt2 columns
    if istable(phoneme_table) && ismember('rt1', phoneme_table.Properties.VariableNames) && ismember('rt2', phoneme_table.Properties.VariableNames)
        % Collect rt1 and rt2 from the current subject's phoneme table
        current_rt1 = phoneme_table.rt1;
        current_rt2 = phoneme_table.rt2;

        % Append valid (non-NaN) values to the overall list
        valid_idx = ~isnan(current_rt1) & ~isnan(current_rt2); % Only keep pairs where both rt1 and rt2 exist
        valid_rt1 = current_rt1(valid_idx);
        valid_rt2 = current_rt2(valid_idx);

        % Calculate mean rt1 and rt2 for the current subject
        mean_rt1 = mean(valid_rt1);
        mean_rt2 = mean(valid_rt2);
        
        % Store the results
        mean_rt1_per_subject = [mean_rt1_per_subject; mean_rt1];
        mean_rt2_per_subject = [mean_rt2_per_subject; mean_rt2];
        subject_ids = [subject_ids; {dbsID}];
    else
        warning(['Subject ', num2str(i_sub), ' does not have the expected columns rt1 and rt2.']);
    end
end

% Display the comparison results for each subject and count occurrences
for i = 1:length(subject_ids)
    dbsID = subject_ids{i};
    mean_rt1 = mean_rt1_per_subject(i);
    mean_rt2 = mean_rt2_per_subject(i);
    
    % Display the mean values
    disp(['DBS ID: ', dbsID]);
    disp(['Mean rt1: ', num2str(mean_rt1)]);
    disp(['Mean rt2: ', num2str(mean_rt2)]);
    
    % Compare rt1 and rt2 and count occurrences
    if mean_rt2 > mean_rt1
        disp('rt2 is greater than rt1 on average.');
        rt2_greater_count = rt2_greater_count + 1;  % Increment rt2 > rt1 count
    else
        disp('rt1 is greater than or equal to rt2 on average.');
        rt1_greater_count = rt1_greater_count + 1;  % Increment rt1 >= rt2 count
    end
    
    disp('--------------------------------------------');
end

% Display the final counts
disp(['Number of times rt2 is greater than rt1: ', num2str(rt2_greater_count)]);
disp(['Number of times rt1 is greater than or equal to rt2: ', num2str(rt1_greater_count)]);
%%
