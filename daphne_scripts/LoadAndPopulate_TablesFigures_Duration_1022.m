clear;

setpaths_dbs_triplet();

% Define the phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'oot', 0, ...
    'oov', 0, ...
    'oog', 0, ...
    'oos', 0, ...
    'aht', 0.0001, ...
    'ahv', 0.0001, ...
    'ahg', 0, ...
    'ahs', 0.0007, ...
    'eet', 0.0002, ...
    'eev', 0.0007, ...
    'eeg', 0.0005, ...
    'ees', 0.0003);

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.syllable_tables = cell(nsubs_to_analyze,1);
subtable.transition_tables = cell(nsubs_to_analyze,1); % New column for storing transitions

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    syllableFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_syllable.txt', dbsID));
   
    % Load syllable table if it exists
    if exist(syllableFilePath, 'file')
        % Load syllable table
        subtable.syllable_tables{i_sub} = readtable(syllableFilePath);

        % Create a table to store transitions, durations, and phonotactic probabilities for this subject
        transTable = table('Size', [0, 8], ... % Start with an empty table and dynamically add rows
            'VariableTypes', {'double', 'double', 'cell', 'cell', 'double', 'double', 'double', 'double'}, ...
            'VariableNames', {'session_id', 'trial_id', 'trans1', 'trans2', 'trans1_duration', 'trans2_duration', 'trans1_prob', 'trans2_prob'});

        % Get unique combinations of session_id and trial_id
        unique_trials = unique(subtable.syllable_tables{i_sub}(:, {'session_id', 'trial_id'}), 'rows');

        % Loop through each trial
        for trial = 1:height(unique_trials)
            session_id = unique_trials.session_id(trial);
            trial_id = unique_trials.trial_id(trial);
            
            % Get the syllables for this session_id and trial_id
            trial_syllables = subtable.syllable_tables{i_sub}(subtable.syllable_tables{i_sub}.session_id == session_id & ...
                                                              subtable.syllable_tables{i_sub}.trial_id == trial_id, :);
            
            % Ensure there are at least 3 syllables to form trans1 and trans2
            if height(trial_syllables) >= 3
                % Extract the 'stim' and 'duration' values
                stim_syl1 = trial_syllables.stim{1}; % stim of syl_id 1
                stim_syl2 = trial_syllables.stim{2}; % stim of syl_id 2
                stim_syl3 = trial_syllables.stim{3}; % stim of syl_id 3

                % Extract real components of durations for each transition component
                dur_v_syl1 = real(trial_syllables.duration_V(1)); % vowel duration of syl_id 1
                dur_c_syl2 = real(trial_syllables.duration_C(2)); % consonant duration of syl_id 2
                dur_v_syl2 = real(trial_syllables.duration_V(2)); % vowel duration of syl_id 2
                dur_c_syl3 = real(trial_syllables.duration_C(3)); % consonant duration of syl_id 3

                % Skip if any duration is NaN
                if any(isnan([dur_v_syl1, dur_c_syl2, dur_v_syl2, dur_c_syl3]))
                    continue;
                end

                % Calculate transitions
                trans1 = strcat(stim_syl1(end-1:end), stim_syl2(1));
                trans2 = strcat(stim_syl2(end-1:end), stim_syl3(1));

                % Calculate the durations for trans1 and trans2
                trans1_duration = dur_v_syl1 + dur_c_syl2;
                trans2_duration = dur_v_syl2 + dur_c_syl3;

                % Retrieve phonotactic probabilities for trans1 and trans2
                if isfield(phonotacticProbabilities, trans1)
                    trans1_prob = phonotacticProbabilities.(trans1);
                else
                    trans1_prob = NaN; % Assign NaN if the transition is not found in the mapping
                end

                if isfield(phonotacticProbabilities, trans2)
                    trans2_prob = phonotacticProbabilities.(trans2);
                else
                    trans2_prob = NaN; % Assign NaN if the transition is not found in the mapping
                end

                % Append the new row to the transTable
                transTable = [transTable; {session_id, trial_id, trans1, trans2, trans1_duration, trans2_duration, trans1_prob, trans2_prob}];
            end
        end

        % Store the transitions, durations, and phonotactic probabilities in the subtable
        subtable.transition_tables{i_sub} = transTable;

        % Write the updated transition table to a file
        writetable(transTable, fullfile(outputDir, sprintf('%s_transition_table_processed.txt', dbsID)));
    end
end

disp('Data processing and file writing complete.');

%%
% Perform correlation analysis between trans1_prob and trans1_duration
significant_dbsIDs_trans1 = {};  % Store DBS IDs with significant correlation for trans1

for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    transition_table = subtable.transition_tables{i_sub};

    % Ensure there are enough valid data points for correlation (no NaN values)
    valid_idx = ~isnan(transition_table.trans1_duration) & ~isnan(transition_table.trans1_prob);
    trans1_duration = transition_table.trans1_duration(valid_idx);
    trans1_prob = transition_table.trans1_prob(valid_idx);

    if length(trans1_duration) > 1  % Ensure there are enough valid data points
        % Perform correlation analysis for trans1_duration and trans1_prob
        [r_trans1, p_trans1] = corr(trans1_duration, trans1_prob, 'Type', 'Pearson');
        disp(['Correlation between trans1_duration and trans1_prob for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans1), ', p = ', num2str(p_trans1)]);
        
        % Check for significance
        if p_trans1 < 0.05
            significant_dbsIDs_trans1{end+1} = dbsID;  % Store significant DBS IDs for trans1
        end
    else
        disp(['Not enough valid data for trans1_duration and trans1_prob correlation for DBSID: ', dbsID]);
    end
end

significant_dbsIDs_trans2 = {};  % Store DBS IDs with significant correlation for trans2

for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    transition_table = subtable.transition_tables{i_sub};

 % Ensure there are enough valid data points for correlation (no NaN values)
    valid_idx = ~isnan(transition_table.trans2_duration) & ~isnan(transition_table.trans2_prob);
    trans2_duration = transition_table.trans2_duration(valid_idx);
    trans2_prob = transition_table.trans2_prob(valid_idx);

    % Perform correlation analysis for pp2 and reaction_time
    if length(trans2_duration) > 1  % Ensure there are enough data points
        [r_trans2, p_trans2] = corr(trans2_duration, trans2_prob, 'Type', 'Pearson');
        disp(['Correlation between pp2 and reaction_time for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans2), ', p = ', num2str(p_trans2)]);
        
        % Check for significance
        if p_trans2 < 0.05
            significant_dbsIDs_trans2{end+1} = dbsID;  % Store significant DBS IDs for pp
        end     
    else
        disp(['Not enough valid data for pp2 and reaction_time correlation for DBSID: ', dbsID]);
    end
end

% % Display significant DBS IDs for trans1
 disp('DBS IDs with significant correlation between trans1_duration and trans1_prob:');
disp(significant_dbsIDs_trans1);

disp('DBS IDs with significant correlation between trans2_duration and reaction_time:');
 disp(significant_dbsIDs_trans2);
%%
% Perform correlation analysis between trans1_prob and trans1_duration using Spearman's correlation
significant_dbsIDs_trans1 = {};  % Store DBS IDs with significant correlation for trans1

for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    transition_table = subtable.transition_tables{i_sub};

    % Ensure there are enough valid data points for correlation (no NaN values)
    valid_idx = ~isnan(transition_table.trans1_duration) & ~isnan(transition_table.trans1_prob);
    trans1_duration = transition_table.trans1_duration(valid_idx);
    trans1_prob = transition_table.trans1_prob(valid_idx);

    if length(trans1_duration) > 1  % Ensure there are enough valid data points
        % Perform Spearman correlation analysis for trans1_duration and trans1_prob
        [r_trans1, p_trans1] = corr(trans1_duration, trans1_prob, 'Type', 'Spearman');
        disp(['Spearman correlation between trans1_duration and trans1_prob for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans1), ', p = ', num2str(p_trans1)]);

        % Check for significance
        if p_trans1 < 0.05
            significant_dbsIDs_trans1{end+1} = dbsID;  % Store significant DBS IDs for trans1
        end
    else
        disp(['Not enough valid data for trans1_duration and trans1_prob correlation for DBSID: ', dbsID]);
    end
end

% Perform correlation analysis between trans2_prob and trans2_duration using Spearman's correlation
significant_dbsIDs_trans2 = {};  % Store DBS IDs with significant correlation for trans2

for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    transition_table = subtable.transition_tables{i_sub};

    % Ensure there are enough valid data points for correlation (no NaN values)
    valid_idx = ~isnan(transition_table.trans2_duration) & ~isnan(transition_table.trans2_prob);
    trans2_duration = transition_table.trans2_duration(valid_idx);
    trans2_prob = transition_table.trans2_prob(valid_idx);

    if length(trans2_duration) > 1  % Ensure there are enough valid data points
        % Perform Spearman correlation analysis for trans2_duration and trans2_prob
        [r_trans2, p_trans2] = corr(trans2_duration, trans2_prob, 'Type', 'Spearman');
        disp(['Spearman correlation between trans2_duration and trans2_prob for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans2), ', p = ', num2str(p_trans2)]);

        % Check for significance
        if p_trans2 < 0.05
            significant_dbsIDs_trans2{end+1} = dbsID;  % Store significant DBS IDs for trans2
        end
    else
        disp(['Not enough valid data for trans2_duration and trans2_prob correlation for DBSID: ', dbsID]);
    end
end

% Display significant DBS IDs for trans1 and trans2
disp('DBS IDs with significant correlation between trans1_duration and trans1_prob:');
disp(significant_dbsIDs_trans1);

disp('DBS IDs with significant correlation between trans2_duration and trans2_prob:');
disp(significant_dbsIDs_trans2);
