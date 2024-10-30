% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir);

% Define the consonant-vowel sequences and their phonotactic probabilities
phonotacticProbabilities = struct(...
    'ghoo', 0.0004, 'voo', 0, 'soo', 0.0023, 'too', 0.0021, ...
    'ghah', 0.002, 'vah', 0.001, 'sah', 0.0022, 'tah', 0.0022, ...
    'ghee', 0.0001, 'vee', 0.001, 'see', 0.0028, 'tee', 0.0012);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.syllable_tables = cell(nsubs_to_analyze, 1);

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    syllableFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_syllable.txt', dbsID));

    % Load syllable table if it exists
    if exist(syllableFilePath, 'file')
        subtable.syllable_tables{i_sub} = readtable(syllableFilePath);

        % Ensure necessary fields exist: starts, duration_C, duration_V
        if all(ismember({'starts', 'duration_C', 'duration_V'}, subtable.syllable_tables{i_sub}.Properties.VariableNames))

            % Add a new column for Reaction Time (initialize with NaN)
            subtable.syllable_tables{i_sub}.Reaction_Time = NaN(height(subtable.syllable_tables{i_sub}), 1);

            % Get unique trial_ids
            unique_trials = unique(subtable.syllable_tables{i_sub}.trial_id);

            % Loop through each trial
            for trial = 1:length(unique_trials)
                trial_id = unique_trials(trial);
                trial_syllable = subtable.syllable_tables{i_sub}(subtable.syllable_tables{i_sub}.trial_id == trial_id, :);

                % Loop through each syllable in the trial
                for idx = 1:height(trial_syllable)
                    % Ensure starts, duration_C, and duration_V values are valid (not NaN or missing)
                    if ~isnan(trial_syllable.starts(idx)) && ~isnan(trial_syllable.duration_C(idx)) && ~isnan(trial_syllable.duration_V(idx))
                        % Calculate the reaction time for each consonant-vowel pair (stimulus onset + duration_C + duration_V)
                        RT = trial_syllable.starts(idx) + trial_syllable.duration_C(idx) + trial_syllable.duration_V(idx);

                        % Store the reaction time in the new 'Reaction_Time' column
                        subtable.syllable_tables{i_sub}.Reaction_Time(subtable.syllable_tables{i_sub}.trial_id == trial_id & ...
                                                                      subtable.syllable_tables{i_sub}.syl_id == trial_syllable.syl_id(idx)) = RT;
                    else
                        % If any value is NaN or missing, issue a warning but don't store NaN
                        warning(['Missing data for trial_id: ', num2str(trial_id), ', syl_id: ', num2str(trial_syllable.syl_id(idx)), ' in DBSID: ', dbsID]);
                    end
                end
            end
        else
            error(['Missing required fields in syllable table for DBSID: ', dbsID]);
        end

        % Save the updated syllable table with the reaction times
        writetable(subtable.syllable_tables{i_sub}, fullfile(outputDir, sprintf('%s_syllable_table_with_RT.txt', dbsID)));
    else
        error(['File does not exist: ', syllableFilePath]);
    end
end

disp('Reaction times calculated and added to syllable tables.');
