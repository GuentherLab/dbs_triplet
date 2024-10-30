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
subtable =  readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.triplet_tables = cell(nsubs_to_analyze,1);
subtable.phonetic_ontarget_scored = false(nsubs_to_analyze,1); 

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

                % Calculate errors and add to the triplet table
                ntrials = height(subtable.triplet_tables{i_sub});
                for itrial = 1:ntrials
                    match_rows = subtable.phoneme_tables{i_sub}.trial_id == subtable.triplet_tables{i_sub}.trial_id(itrial);

                    % Initialize error counts for the three syllables
                    errors_syl1 = 0;
                    errors_syl2 = 0;
                    errors_syl3 = 0;

                    % Loop through each syllable (1 to 3)
                    for isyl = 1:3
                        % Get the rows corresponding to this syl_id
                        syl_rows = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == isyl);

                        % Count errors for consonant and vowel
                        consonant_error = nnz(subtable.phoneme_tables{i_sub}.phonetic_ontarget(syl_rows & strcmp(subtable.phoneme_tables{i_sub}.type, 'consonant')) == 0);
                        vowel_error = nnz(subtable.phoneme_tables{i_sub}.phonetic_ontarget(syl_rows & strcmp(subtable.phoneme_tables{i_sub}.type, 'vowel')) == 0);

                        % Sum errors for the syllable
                        total_errors = consonant_error + vowel_error;

                        % Store the error count in the appropriate variable
                        if isyl == 1
                            errors_syl1 = total_errors;
                        elseif isyl == 2
                            errors_syl2 = total_errors;
                        elseif isyl == 3
                            errors_syl3 = total_errors;
                        end
                    end

                    % Sum the errors across the three syllables for this trial
                    total_errors_per_trial = errors_syl1 + errors_syl2 + errors_syl3;
                    
                    % Store these values in the triplet table
                    subtable.triplet_tables{i_sub}.errors_syl1(itrial) = errors_syl1;
                    subtable.triplet_tables{i_sub}.errors_syl2(itrial) = errors_syl2;
                    subtable.triplet_tables{i_sub}.errors_syl3(itrial) = errors_syl3;
                    subtable.triplet_tables{i_sub}.num_error_phonemes(itrial) = total_errors_per_trial; % total errors for the trial
                end

                % Write the processed phoneme and triplet tables to separate files
                writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
                writetable(subtable.triplet_tables{i_sub}, fullfile(outputDir, sprintf('%s_triplet_table_processed.txt', dbsID)));
            end
        end
    end
end

disp('Data processing and file writing complete.');
