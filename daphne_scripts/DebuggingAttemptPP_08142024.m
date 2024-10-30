clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% List of specified DBS IDs
dbsIDs = {'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
          'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', ...
          'DBS3029', 'DBS3030', 'DBS3031'};

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);
subtable.triplet_tables = cell(nsubs_to_analyze, 1);
subtable.phonetic_ontarget_scored = false(nsubs_to_analyze, 1); 

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

% Loop through each DBS ID
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    
    % Check if the current dbsID is in the list of specified DBS IDs
    if ismember(dbsID, dbsIDs)
        phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
        tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));

        % Load phoneme and triplet tables if they exist
        if exist(phonemeFilePath, 'file')
            subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
            subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

            if ismember('phonetic_ontarget', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
                subtable.phonetic_ontarget_scored(i_sub) = true;

                % Initialize new columns in the phoneme table
                subtable.phoneme_tables{i_sub}.phonotactic_prob_orig1 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
                subtable.phoneme_tables{i_sub}.phonotactic_prob_orig2 = NaN(height(subtable.phoneme_tables{i_sub}), 1);

                % Loop through each trial
                unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);
                for itrial = 1:length(unique_trials)
                    trial_id = unique_trials(itrial);
                    match_rows = subtable.phoneme_tables{i_sub}.trial_id == trial_id;

                    % Extract the relevant rows for each syllable
                    syl1_vowel_row = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == 1) & strcmp(subtable.phoneme_tables{i_sub}.type, 'vowel');
                    syl2_consonant_row = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == 2) & strcmp(subtable.phoneme_tables{i_sub}.type, 'consonant');
                    syl2_vowel_row = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == 2) & strcmp(subtable.phoneme_tables{i_sub}.type, 'vowel');
                    syl3_consonant_row = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == 3) & strcmp(subtable.phoneme_tables{i_sub}.type, 'consonant');

                    % Combine the stimulus strings for transition 1 and transition 2
                    if any(syl1_vowel_row) && any(syl2_consonant_row)
                        combined_str_orig1 = strcat(subtable.phoneme_tables{i_sub}.stim{syl1_vowel_row}, subtable.phoneme_tables{i_sub}.stim{syl2_consonant_row});
                        if isfield(phonotacticProbabilities, combined_str_orig1)
                            subtable.phoneme_tables{i_sub}.phonotactic_prob_orig1(syl2_consonant_row) = phonotacticProbabilities.(combined_str_orig1);
                        else
                            subtable.phoneme_tables{i_sub}.phonotactic_prob_orig1(syl2_consonant_row) = NaN;
                        end
                    end
                    if any(syl2_vowel_row) && any(syl3_consonant_row)
                        combined_str_orig2 = strcat(subtable.phoneme_tables{i_sub}.stim{syl2_vowel_row}, subtable.phoneme_tables{i_sub}.stim{syl3_consonant_row});
                        if isfield(phonotacticProbabilities, combined_str_orig2)
                            subtable.phoneme_tables{i_sub}.phonotactic_prob_orig2(syl3_consonant_row) = phonotacticProbabilities.(combined_str_orig2);
                        else
                            subtable.phoneme_tables{i_sub}.phonotactic_prob_orig2(syl3_consonant_row) = NaN;
                        end
                    end
                end

                % Write the processed phoneme tables to separate files
                writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
            end
        end
    end
end

disp('Data processing and file writing complete.');
