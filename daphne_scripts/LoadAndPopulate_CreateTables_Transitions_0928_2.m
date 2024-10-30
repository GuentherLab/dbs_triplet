% % %This code correctly assigns trans1 and trans2 and rt1 and rt2

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

        % Add empty columns for trans1, trans2, rt1, and rt2
        subtable.phoneme_tables{i_sub}.trans1 = cell(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.trans2 = cell(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.rt1 = NaN(height(subtable.phoneme_tables{i_sub}), 1);
        subtable.phoneme_tables{i_sub}.rt2 = NaN(height(subtable.phoneme_tables{i_sub}), 1);

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

                % Assign trans1, trans2, rt1, and rt2 only to the first row of the trial
                subtable.phoneme_tables{i_sub}.trans1(row_indices(1)) = {trans1};
                subtable.phoneme_tables{i_sub}.trans2(row_indices(1)) = {trans2};
                subtable.phoneme_tables{i_sub}.rt1(row_indices(1)) = rt1;
                subtable.phoneme_tables{i_sub}.rt2(row_indices(1)) = rt2;

                % Set trans1, trans2, rt1, and rt2 to NaN for all other rows in the same trial
                subtable.phoneme_tables{i_sub}.trans1(row_indices(2:end)) = {NaN};
                subtable.phoneme_tables{i_sub}.trans2(row_indices(2:end)) = {NaN};
                subtable.phoneme_tables{i_sub}.rt1(row_indices(2:end)) = NaN;
                subtable.phoneme_tables{i_sub}.rt2(row_indices(2:end)) = NaN;
            end
        end

        % Write the updated phoneme table with trans1, trans2, rt1, and rt2 to a file
        writetable(subtable.phoneme_tables{i_sub}, fullfile(outputDir, sprintf('%s_phoneme_table_processed.txt', dbsID)));
        writetable(subtable.triplet_tables{i_sub}, fullfile(outputDir, sprintf('%s_triplet_table_processed.txt', dbsID)));
    end
end

disp('Data processing and file writing complete.');
