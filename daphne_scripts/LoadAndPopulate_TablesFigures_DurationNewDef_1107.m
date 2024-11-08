clear;
setpaths_dbs_triplet();

% Define the phonotactic probabilities mapping
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

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.transition_tables = cell(nsubs_to_analyze,1); % New column for storing transitions

clear;
setpaths_dbs_triplet();

% Define the base directory
baseDir = 'Z:\DBS';

%Code for determining which 'subjects' have 'starts' and 'ends' and do not 
% Initialize lists for DBS IDs with and without 'starts' and 'ends' columns
dbsIDs_with_starts_ends = {};
dbsIDs_without_starts_ends = {};

% Load the main table with DBS IDs
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);

% Loop through each DBS ID to check for 'starts' and 'ends' columns
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Check if the file exists
    if exist(phonemeFilePath, 'file')
        phonemeTable = readtable(phonemeFilePath);
        
        % Check if 'starts' and 'ends' columns exist
        if all(ismember({'starts', 'ends'}, phonemeTable.Properties.VariableNames))
            dbsIDs_with_starts_ends{end+1} = dbsID;
        else
            dbsIDs_without_starts_ends{end+1} = dbsID;
        end
    else
        disp(['File not found for DBS ID ', dbsID, '.']);
    end
end

% Display results
disp('DBS IDs with "starts" and "ends" columns:');
disp(dbsIDs_with_starts_ends);

disp('DBS IDs without "starts" and "ends" columns:');
disp(dbsIDs_without_starts_ends);
%% Using Duration instead 

% Loop through each DBS ID
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Load phoneme data if it exists
    if exist(phonemeFilePath, 'file')
        phonemeTable = readtable(phonemeFilePath);
        subtable.phoneme_tables{i_sub} = phonemeTable;

        % Initialize table to store transitions for this subject
        transTable = table('Size', [0, 8], ...
                           'VariableTypes', {'double', 'double', 'cell', 'cell', 'double', 'double', 'double', 'double'}, ...
                           'VariableNames', {'session_id', 'trial_id', 'trans1', 'trans2', 'trans1_duration', 'trans2_duration', 'trans1_prob', 'trans2_prob'});

        % Get unique session and trial IDs
        unique_trials = unique(phonemeTable.trial_id);
        unique_sessions = unique(phonemeTable.session_id);

        % Loop through each session and trial
        for session = unique_sessions'
            for trial = unique_trials'
                % Get relevant rows for this session and trial
                match_rows = (phonemeTable.session_id == session) & (phonemeTable.trial_id == trial);
                phonemes = phonemeTable(match_rows, :);

                % Ensure we have at least five rows to form trans1 and trans2
                if height(phonemes) >= 5
                    % Initialize transitions, durations, and probabilities as empty or NaN
                    trans1 = ''; trans1_duration = NaN; trans1_prob = NaN;
                    trans2 = ''; trans2_duration = NaN; trans2_prob = NaN;

                    % Define trans1 as the concatenation of 'stim' from rows 2 and 3
                    % Define trans2 as the concatenation of 'stim' from rows 4 and 5
                    trans1 = strcat(phonemes.stim{2}, phonemes.stim{3});
                    trans1_duration = phonemes.duration(2) + phonemes.duration(3); % Sum durations from rows 2 and 3

                    trans2 = strcat(phonemes.stim{4}, phonemes.stim{5});
                    trans2_duration = phonemes.duration(4) + phonemes.duration(5); % Sum durations from rows 4 and 5

                    % Retrieve phonotactic probabilities
                    if isfield(phonotacticProbabilities, trans1)
                        trans1_prob = phonotacticProbabilities.(trans1);
                    else
                        trans1_prob = NaN; % Assign NaN if probability is missing
                    end

                    if isfield(phonotacticProbabilities, trans2)
                        trans2_prob = phonotacticProbabilities.(trans2);
                    else
                        trans2_prob = NaN; % Assign NaN if probability is missing
                    end

                    % Append this trial's transitions, durations, and probabilities to the transTable
                    transTable = [transTable; {session, trial, trans1, trans2, trans1_duration, trans2_duration, trans1_prob, trans2_prob}];
                    
                    % Debug output
                    fprintf('DBS ID: %s, Session ID: %d, Trial ID: %d, Trans1: %s, Trans1 Duration: %f, Trans1 Probability: %f, Trans2: %s, Trans2 Duration: %f, Trans2 Probability: %f\n', ...
                            dbsID, session, trial, trans1, trans1_duration, trans1_prob, trans2, trans2_duration, trans2_prob);
                end
            end
        end

        % Store and save the transition table if it contains data
        if ~isempty(transTable)
            subtable.transition_tables{i_sub} = transTable;
            writetable(transTable, fullfile(outputDir, sprintf('%s_transition_table_processed.txt', dbsID)));
        else
            fprintf('No valid transitions found for DBS ID %s.\n', dbsID);
        end
    else
        disp(['File not found for DBS ID ', dbsID, '.']);
    end
end

disp('Data processing and file writing complete.');

%%
%% Correlation Analysis

% Initialize lists to store significant DBS IDs for trans1 and trans2 correlations
significant_dbsIDs_trans1 = {};  % Store DBS IDs with significant correlation for trans1
significant_dbsIDs_trans2 = {};  % Store DBS IDs with significant correlation for trans2

% Perform correlation analysis between trans1_prob and trans1_duration, trans2_prob and trans2_duration
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    transTable = subtable.transition_tables{i_sub};

    % Check if transTable is empty or not a table, and skip if so
    if isempty(transTable) || ~istable(transTable)
        continue;
    end

    % Filter valid rows for trans1 and trans2
    valid_rows_trans1 = ~isnan(transTable.trans1_duration) & ~isnan(transTable.trans1_prob);
    valid_rows_trans2 = ~isnan(transTable.trans2_duration) & ~isnan(transTable.trans2_prob);

    % Get valid data for trans1
    trans1_durations = transTable.trans1_duration(valid_rows_trans1);
    trans1_probs = transTable.trans1_prob(valid_rows_trans1);

    % Get valid data for trans2
    trans2_durations = transTable.trans2_duration(valid_rows_trans2);
    trans2_probs = transTable.trans2_prob(valid_rows_trans2);

    % Perform correlation analysis for trans1_prob and trans1_duration
    if length(trans1_probs) > 1  % Ensure there are enough data points
        [r_trans1, p_trans1] = corr(trans1_probs, trans1_durations, 'Type', 'Pearson');
        disp(['Correlation between trans1_prob and trans1_duration for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans1), ', p = ', num2str(p_trans1)]);
        
        % Check for significance
        if p_trans1 < 0.05
            significant_dbsIDs_trans1{end+1} = dbsID;  % Store significant DBS IDs for trans1
        end
    else
        disp(['Not enough valid data for trans1_prob and trans1_duration correlation for DBSID: ', dbsID]);
    end

    % Perform correlation analysis for trans2_prob and trans2_duration
    if length(trans2_probs) > 1  % Ensure there are enough data points
        [r_trans2, p_trans2] = corr(trans2_probs, trans2_durations, 'Type', 'Pearson');
        disp(['Correlation between trans2_prob and trans2_duration for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_trans2), ', p = ', num2str(p_trans2)]);
        
        % Check for significance
        if p_trans2 < 0.05
            significant_dbsIDs_trans2{end+1} = dbsID;  % Store significant DBS IDs for trans2
        end
    else
        disp(['Not enough valid data for trans2_prob and trans2_duration correlation for DBSID: ', dbsID]);
    end
end

% Display significant DBS IDs for trans1 and trans2
disp('DBS IDs with significant correlation between trans1_prob and trans1_duration:');
disp(significant_dbsIDs_trans1);

disp('DBS IDs with significant correlation between trans2_prob and trans2_duration:');
disp(significant_dbsIDs_trans2);
