% Duration Analysis code here of the transitions 
%This code checks for NaN in the 'duration' column and excludes trials with NaN from the duration_tables.
clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Define the updated phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'oot', 0, ...
    'oov', 0, ...
    'oogh', 0, ...
    'oos', 0, ...
    'aht', 0.0001, ...
    'ahv', 0.0001, ...
    'ahgh', 0, ...
    'ahs', 0.0007, ...
    'eet', 0.0002, ...
    'eev', 0.0007, ...
    'eegh', 0.0005, ...
    'ees', 0.0003);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.duration_tables = cell(nsubs_to_analyze,1); % New table to store transitions, durations, and probabilities

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Load phoneme table if it exists
    if exist(phonemeFilePath, 'file')
        % Load phoneme table
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);

        % Create an empty duration table for this subject
        duration_table = table();

        % Initialize empty arrays for storing trial-specific data
        trial_ids = [];
        trans1_vals = {};
        trans2_vals = {};
        pp1_vals = [];
        pp2_vals = [];
        dur_trans1_vals = [];
        dur_trans2_vals = [];

        % Get unique trial_ids
        unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);

        % Loop through each trial
        for trial = 1:length(unique_trials)
            trial_id = unique_trials(trial);

            % Get the phonemes for this trial
            trial_phonemes = subtable.phoneme_tables{i_sub}(subtable.phoneme_tables{i_sub}.trial_id == trial_id, :);

            % Check if any NaN values are present in the 'duration' column for this trial
            if any(isnan(trial_phonemes.duration))
                % If any NaN is found, skip this trial
                disp(['Skipping trial ', num2str(trial_id), ' due to NaN in duration.']);
                continue;
            end

            % Ensure there are enough phonemes (at least 3) to form trans1 and trans2
            if height(trial_phonemes) >= 3
                % Identify the required syllables for trans1 and trans2
                vowel1_idx = find(strcmp(trial_phonemes.type, 'vowel') & trial_phonemes.syl_id == 1, 1, 'first');
                consonant1_idx = find(strcmp(trial_phonemes.type, 'consonant') & trial_phonemes.syl_id == 2, 1, 'first');

                vowel2_idx = find(strcmp(trial_phonemes.type, 'vowel') & trial_phonemes.syl_id == 2, 1, 'first');
                consonant2_idx = find(strcmp(trial_phonemes.type, 'consonant') & trial_phonemes.syl_id == 3, 1, 'first');

                % Check if the indices are valid and exist
                if ~isempty(vowel1_idx) && ~isempty(consonant1_idx)
                    % Construct trans1 (vowel from syl_id 1 + consonant from syl_id 2)
                    trans1 = strcat(trial_phonemes.stim{vowel1_idx}, trial_phonemes.stim{consonant1_idx});
                    % Calculate dur_trans1: sum of durations for vowel1 and consonant1
                    dur_trans1 = trial_phonemes.duration(vowel1_idx) + trial_phonemes.duration(consonant1_idx);
                    % Assign phonotactic probability for trans1
                    pp1 = NaN;
                    if isfield(phonotacticProbabilities, trans1)
                        pp1 = phonotacticProbabilities.(trans1);
                    end
                end

                if ~isempty(vowel2_idx) && ~isempty(consonant2_idx)
                    % Construct trans2 (vowel from syl_id 2 + consonant from syl_id 3)
                    trans2 = strcat(trial_phonemes.stim{vowel2_idx}, trial_phonemes.stim{consonant2_idx});
                    % Calculate dur_trans2: sum of durations for vowel2 and consonant2
                    dur_trans2 = trial_phonemes.duration(vowel2_idx) + trial_phonemes.duration(consonant2_idx);
                    % Assign phonotactic probability for trans2
                    pp2 = NaN;
                    if isfield(phonotacticProbabilities, trans2)
                        pp2 = phonotacticProbabilities.(trans2);
                    end
                end

                % Store the valid data for this trial
                trial_ids = [trial_ids; trial_id];
                trans1_vals = [trans1_vals; trans1];
                trans2_vals = [trans2_vals; trans2];
                pp1_vals = [pp1_vals; pp1];
                pp2_vals = [pp2_vals; pp2];
                dur_trans1_vals = [dur_trans1_vals; dur_trans1];
                dur_trans2_vals = [dur_trans2_vals; dur_trans2];
            end
        end

        % Only add trials without NaN in duration to the final duration_table
        duration_table.trial_id = trial_ids;
        duration_table.trans1 = trans1_vals;
        duration_table.trans2 = trans2_vals;
        duration_table.pp1 = pp1_vals;
        duration_table.pp2 = pp2_vals;
        duration_table.dur_trans1 = dur_trans1_vals;
        duration_table.dur_trans2 = dur_trans2_vals;

        % Store the duration table for this subject
        subtable.duration_tables{i_sub} = duration_table;

        % Write the updated duration table to a file for this subject
        writetable(duration_table, fullfile(outputDir, sprintf('%s_duration_table.txt', dbsID)));
    end
end

disp('Transition extraction, phonotactic probability assignment, and duration calculation complete, stored in duration_tables for each subject.');
