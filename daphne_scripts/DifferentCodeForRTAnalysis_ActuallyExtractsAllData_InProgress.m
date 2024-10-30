clear;

% Define the base directory containing DBS folders
baseDir = 'Z:\DBS';
outputDir = 'C:\Users\amsmeier\dbs_triplet\daphne_scripts';  % Set your desired output directory
mkdir(outputDir);

% Get a list of all DBS folders in the base directory
dbsFolders = dir(fullfile(baseDir, 'DBS*'));  % Adjust this if needed for folder naming convention

% Initialize a cell array to hold the valid DBS IDs and their corresponding data
validDBSIDs = {};
subtable_phonemes = {};  % To hold the valid DBS IDs and their corresponding data

% Define headers for the output table
headers = {'onset', 'duration', 'session_id', 'trial_id', 'syl_id', 'type', 'stim'};

% Loop through each DBS folder
for i_sub = 1:length(dbsFolders)
    if dbsFolders(i_sub).isdir  % Check if it is a directory
        dbsID = dbsFolders(i_sub).name;  % Get the name of the folder
        
        % Construct the file path for the phoneme data
        phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
        
        % Check if the file exists
        if exist(phonemeFilePath, 'file')
            % Load the phoneme data
            data = readtable(phonemeFilePath);
            
            % Check for required variables
            if ismember('onset', data.Properties.VariableNames) && ismember('duration', data.Properties.VariableNames)
                validDBSIDs{end+1} = dbsID;  % Add to valid DBS IDs
                
                % Convert the first 7 columns of data to a cell array
                dataSubset = table2cell(data(:, 1:7));  % Convert to cell array
                
                % Create a new cell array for the combined data
                combinedData = [headers; dataSubset];  % Stack headers on top of data

                % Store the DBS ID and its corresponding combined data
                subtable_phonemes(end+1, :) = [{dbsID}, {combinedData}];  % Wrap combinedData in a cell
            end
        end
    end
end

% Create a new table with headers
subtable_phonemes_table = cell2table(subtable_phonemes, 'VariableNames', {'DBS_ID', 'Data'});

% Optionally save the subtable to a file
if ~isempty(subtable_phonemes_table)
    writetable(subtable_phonemes_table, fullfile(outputDir, 'valid_phoneme_data_with_headers.txt'));
end

disp('Data processing complete.');
disp('Valid DBS IDs:');
disp(validDBSIDs);  % Display valid DBS IDs
