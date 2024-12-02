clear;
setpaths_dbs_triplet();

% Define the base directory and output directory
baseDir = 'Z:\DBS';
outputDir = fullfile(baseDir, 'Outputs');
mkdir(outputDir);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);

% Initialize columns for storing phoneme tables
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);

% Loop through each DBS ID
for i_sub = 1:nsubs_to_analyze
    % Get DBS ID and corresponding file path
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Check if the phoneme file exists
    if exist(phonemeFilePath, 'file')
        % Load the phoneme table
        phonemeTable = readtable(phonemeFilePath);

        % Check for missing 'starts' and 'ends' and calculate if needed
        if ~all(ismember({'starts', 'ends'}, phonemeTable.Properties.VariableNames))
            phonemeTable.starts = phonemeTable.onset; % Starts are the same as onsets
            phonemeTable.ends = phonemeTable.onset + phonemeTable.duration; % Ends are onset + duration
            disp(['Calculated "starts" and "ends" for DBS ID: ', dbsID]);
        else
            disp(['"Starts" and "Ends" already exist for DBS ID: ', dbsID]);
        end

        % Store the updated phoneme table in the subtable
        subtable.phoneme_tables{i_sub} = phonemeTable;

        % Save the updated phoneme table to a .txt file
        outputFilePath = fullfile(outputDir, sprintf('%s_updated_phoneme_table.txt', dbsID));
        writetable(phonemeTable, outputFilePath, 'Delimiter', '\t');
        disp(['Updated phoneme table saved for DBS ID: ', dbsID, ' at ', outputFilePath]);
    else
        disp(['Phoneme file not found for DBS ID: ', dbsID]);
    end
end

% Save the updated subtable with phoneme data
subtableFilePath = fullfile(outputDir, 'Updated_Subtable.txt');
writetable(subtable, subtableFilePath, 'Delimiter', '\t');
disp(['Subtable saved with updated phoneme data at ', subtableFilePath]);
