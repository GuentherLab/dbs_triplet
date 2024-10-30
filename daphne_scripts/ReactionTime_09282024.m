clear;

% Set paths for the project
setpaths_dbs_triplet();

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir);

% List of specific DBS IDs to process
dbsIDs = {
    'DBS4087', 'DBS4086', 'DBS4085', 'DBS4084', 'DBS4083', ...
    'DBS4082', 'DBS4081', 'DBS4080', 'DBS4079', 'DBS4078', ...
    'DBS4077', 'DBS4076', 'DBS4075', 'DBS4074', 'DBS4073', ...
    'DBS4072', 'DBS4071', 'DBS4070', 'DBS4069', 'DBS4068', ...
    'DBS4067', 'DBS4066', 'DBS4062', 'DBS4061', 'DBS4060', ...
    'DBS4058', 'DBS4057', 'DBS3031', 'DBS3030', 'DBS3029', ...
    'DBS3028', 'DBS3027', 'DBS3026', 'DBS3025', 'DBS3024', ...
    'DBS3023', 'DBS3022', 'DBS3021', 'DBS3020', 'DBS3019', ...
    'DBS3018', 'DBS3017', 'DBS3016', 'DBS3015', 'DBS3014', ...
    'DBS3012', 'DBS3011', 'DBS3010', 'DBS3008', 'DBS3006', ...
    'DBS3005', 'DBS3004', 'DBS3003', 'DBS3002', 'DBS3001'
};

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

% Loop through the specified DBS IDs
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    
    % Construct file paths for syllable data (produced_syllable.txt)
    syllableFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_syllable.txt', dbsID));
    
    % Check if the file exists
    if exist(syllableFilePath, 'file')
        % Load the syllable data
        syllableData = readtable(syllableFilePath);
        
        disp(['Loaded syllable data for DBS ID: ', dbsID]);
        disp(head(syllableData));  % Display the first few rows

        % Initialize new columns in the syllable data
        syllableData.trans1 = strings(height(syllableData), 1);  % Transition 1
        syllableData.trans2 = strings(height(syllableData), 1);  % Transition 2
        syllableData.RT_1 = NaN(height(syllableData), 1);        % Reaction Time 1
        syllableData.RT_2 = NaN(height(syllableData), 1);        % Reaction Time 2
        syllableData.phonotactic_prob_trans_1 = NaN(height(syllableData), 1);
        syllableData.phonotactic_prob_trans_2 = NaN(height(syllableData), 1);
        
        % Loop through each trial
        unique_trials = unique(syllableData.trial_id);
        for itrial = 1:length(unique_trials)
            trial_rows = syllableData(syllableData.trial_id == unique_trials(itrial), :);
            
            % Print trial rows for debugging
            disp(['Trial ID: ', num2str(unique_trials(itrial))]);
            disp(trial_rows);  % Show rows for this trial

            % Check if we have exactly 3 rows for processing
            if height(trial_rows) ~= 3
                warning('Trial does not have exactly three rows for DBS ID: %s', dbsID);
                continue;  % Skip to the next trial if not
            end
            
            % Extract the stimuli
            stim1 = trial_rows.stim{1};
            stim2 = trial_rows.stim{2};
            stim3 = trial_rows.stim{3};

            % Calculate transitions
            if ~isempty(stim1) && ~isempty(stim2) && ~isempty(stim3)
                trans1 = strcat(stim1(end-1:end), stim2(1));  % e.g., "oot"
                trans2 = strcat(stim2(end-1:end), stim3(1));  % e.g., "eeg"
            else
                warning('One of the stimuli is empty for DBS ID: %s', dbsID);
                continue;  % Skip to the next trial
            end
            
            disp(['Stimuli: ', stim1, ', ', stim2, ', ', stim3]);
            disp(['Transitions: ', trans1, ', ', trans2]);

            % Store the transitions
            syllableData.trans1(syllableData.trial_id == unique_trials(itrial)) = trans1;
            syllableData.trans2(syllableData.trial_id == unique_trials(itrial)) = trans2;
            
            % Assign phonotactic probabilities
            if isfield(phonotacticProbabilities, trans1)
                syllableData.phonotactic_prob_trans_1(syllableData.trial_id == unique_trials(itrial)) = phonotacticProbabilities.(trans1);
                probTrans1 = syllableData.phonotactic_prob_trans_1(syllableData.trial_id == unique_trials(itrial));
                disp(['Phonotactic Prob Trans1: ', num2str(probTrans1(1))]);  % Display first value
            end
            if isfield(phonotacticProbabilities, trans2)
                syllableData.phonotactic_prob_trans_2(syllableData.trial_id == unique_trials(itrial)) = phonotacticProbabilities.(trans2);
                probTrans2 = syllableData.phonotactic_prob_trans_2(syllableData.trial_id == unique_trials(itrial));
                disp(['Phonotactic Prob Trans2: ', num2str(probTrans2(1))]);  % Display first value
            end
            
            % Calculate RT_1 and RT_2 using the 'starts' and 'ends' columns
            RT_1 = trial_rows.starts(2) - trial_rows.ends(1);  % Reaction time from end of syl_id 1 to start of syl_id 2
            RT_2 = trial_rows.starts(3) - trial_rows.ends(2);  % Reaction time from end of syl_id 2 to start of syl_id 3
            
            disp(['RT_1: ', num2str(RT_1), ', RT_2: ', num2str(RT_2)]);

            % Store reaction times
            syllableData.RT_1(syllableData.trial_id == unique_trials(itrial)) = RT_1;
            syllableData.RT_2(syllableData.trial_id == unique_trials(itrial)) = RT_2;
        end
        
        % Save the processed syllable table with transitions, reaction times, and phonotactic probabilities
        outputFile = fullfile(outputDir, sprintf('%s_syllable_table_processed.txt', dbsID));
        writetable(syllableData, outputFile);
        
        % Correlation analysis between RT_1, RT_2 and phonotactic probabilities
        valid_idx1 = ~isnan(syllableData.RT_1) & ~isnan(syllableData.phonotactic_prob_trans_1);
        valid_idx2 = ~isnan(syllableData.RT_2) & ~isnan(syllableData.phonotactic_prob_trans_2);
        
        if sum(valid_idx1) > 1  % Ensure there are enough data points
            [corr_RT1, p_RT1] = corr(syllableData.RT_1(valid_idx1), syllableData.phonotactic_prob_trans_1(valid_idx1), 'Type', 'Pearson');
            fprintf('DBS ID: %s, Correlation RT1 vs Phonotactic Prob: %.2f, p-value: %.4f\n', dbsID, corr_RT1, p_RT1);
        end
        if sum(valid_idx2) > 1
            [corr_RT2, p_RT2] = corr(syllableData.RT_2(valid_idx2), syllableData.phonotactic_prob_trans_2(valid_idx2), 'Type', 'Pearson');
            fprintf('DBS ID: %s, Correlation RT2 vs Phonotactic Prob: %.2f, p-value: %.4f\n', dbsID, corr_RT2, p_RT2);
        end
    else
        warning('File does not exist for DBS ID: %s', dbsID);
    end
end

disp('Reaction time calculation and correlation analysis complete.');
