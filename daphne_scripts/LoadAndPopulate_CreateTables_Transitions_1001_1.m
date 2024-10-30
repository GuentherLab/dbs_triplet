clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);

% Create fields to store data for each subject
subtable.stimulus_tables = cell(nsubs_to_analyze,1);  % For stimulus syllable or triplet data
subtable.producedsyllables_tables = cell(nsubs_to_analyze,1); % For produced syllable data
subtable.reactiontime_tables = cell(nsubs_to_analyze,1); % Placeholder for reaction time analysis later

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    
    % Paths to stimulus syllable and triplet text files
    stimulusSyllableFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_stimulus_syllable.txt', dbsID)); 
    stimulusTripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_stimulus_triplet.txt', dbsID)); 
    producedsyllablesFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_syllable.txt', dbsID)); 
   
    % Check if the stimulus syllable or triplet file exists and load the appropriate file
    if exist(stimulusSyllableFilePath, 'file')
        % Load stimulus syllable table if it exists
        subtable.stimulus_tables{i_sub} = readtable(stimulusSyllableFilePath);
        disp(['Loaded stimulus syllable file for DBSID: ', dbsID]);
    elseif exist(stimulusTripletFilePath, 'file')
        % Load stimulus triplet table if the syllable file doesn't exist
        subtable.stimulus_tables{i_sub} = readtable(stimulusTripletFilePath);
        disp(['Loaded stimulus triplet file for DBSID: ', dbsID]);
    else
        % If neither file exists, display a warning and skip to the next subject
        warning(['No stimulus file (syllable or triplet) found for DBSID: ', dbsID]);
        continue;
    end
    
    % Check if the produced syllable file exists and load it
    if exist(producedsyllablesFilePath, 'file')
        subtable.producedsyllables_tables{i_sub} = readtable(producedsyllablesFilePath);
        disp(['Loaded produced syllable file for DBSID: ', dbsID]);
    else
        % If the file doesn't exist, display a warning
        warning(['Produced syllable file not found for DBSID: ', dbsID]);
        continue;
    end
    
    % Initialize a reaction time table for further analysis (this can be used later)
    subtable.reactiontime_tables{i_sub} = table(cell(0,1), cell(0,1), 'VariableNames', {'reactiontime', 'trans2'});
    
    % Get unique trial_ids from the stimulus syllable or triplet data
    unique_trials = unique(subtable.stimulus_tables{i_sub}.trial_id);

    % Loop through each trial (to be used for future processing if needed)
    for trial = 1:length(unique_trials)
        trial_id = unique_trials(trial);
        
        % Filter out the data for the current trial from both tables
        trial_stimulus = subtable.stimulus_tables{i_sub}(subtable.stimulus_tables{i_sub}.trial_id == trial_id, :);
        trial_syllables = subtable.producedsyllables_tables{i_sub}(subtable.producedsyllables_tables{i_sub}.trial_id == trial_id, :);

        % Further analysis can be added here (e.g., reaction time calculations, phoneme processing)
    end
end

disp('Data loading complete.');
