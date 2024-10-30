clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% List of DBS IDs
% % % % % % % % % % % % % % % % % % % % % % % % % dbsIDs = {'DBS3002', 'DBS3004', 'DBS3006', 'DBS3010', 'DBS3011', 'DBS3014', ...
% % % % % % % % % % % % % % % % % % % % % % % % %           'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
% % % % % % % % % % % % % % % % % % % % % % % % %           'DBS3032', 'DBS4079', 'DBS4080'};

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.triplet_tables = cell(nsubs_to_analyze,1);
subtable.transition_tables = cell(nsubs_to_analyze,1); % New transition table

% Create empty tables with the specified columns
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));
   
    % Load phoneme and triplet tables if they exist
    if exist(phonemeFilePath, 'file')
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);
        
        % Initialize the transition table for this subject
        transition_table = table(cell(0,1), cell(0,1), 'VariableNames', {'trans1', 'trans2'});

        % Get unique trial_ids
        unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);

        % Loop through each trial
        for trial = 1:length(unique_trials)
            trial_id = unique_trials(trial);
            
            % Get the phonemes for this trial
            trial_phonemes = subtable.phoneme_tables{i_sub}(subtable.phoneme_tables{i_sub}.trial_id == trial_id, :);
            
            % Identify consonants and vowels
            consonants = trial_phonemes(strcmp(trial_phonemes.stim, 'consonant'), :);
            vowels = trial_phonemes(strcmp(trial_phonemes.stim, 'vowel'), :);

            % Ensure we have enough consonants and vowels
            if height(consonants) >= 3 && height(vowels) >= 2
                % trans1: first vowel + second consonant
                trans1 = strcat(vowels.stim{1}, consonants.stim{2});
                disp(['trans1: ', trans1]);

                % trans2: second vowel + third consonant
                trans2 = strcat(vowels.stim{2}, consonants.stim{3});
                disp(['trans2: ', trans2]);

                % Add to the transition table
                transition_table = [transition_table; {trans1, trans2}];
            end
        end

        % Store the transition table for the current subject
        subtable.transition_tables{i_sub} = transition_table;

        % Write the processed phoneme, triplet, and transition tables to separate files
        writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
        writetable(subtable.triplet_tables{i_sub}, fullfile(outputDir, sprintf('%s_triplet_table_processed.txt', dbsID)));
        writetable(subtable.transition_tables{i_sub}, fullfile(outputDir, sprintf('%s_transition_table_processed.txt', dbsID))); % Write transition table
    end
end

disp('Data processing and file writing complete.');
