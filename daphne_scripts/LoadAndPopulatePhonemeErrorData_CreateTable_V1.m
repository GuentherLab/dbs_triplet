% Define the base directory
baseDir = 'Z:\DBS';

% List of DBS IDs
dbsIDs = {'DBS3002', 'DBS3004', 'DBS3006', 'DBS3010', 'DBS3011', 'DBS3014', ...
          'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3032', 'DBS4079', 'DBS4080'};

% Output file path
outputFilePath = fullfile(baseDir, 'trialtable_produced_phoneme_table.txt');
fileID = fopen(outputFilePath, 'w');  % Open the file for writing
if fileID == -1
    error('File could not be opened for writing.');
end

% Loop through each DBS ID
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Check if the phoneme file exists
    if exist(phonemeFilePath, 'file')
        % Read the content of the phoneme file
        fileContent = fileread(phonemeFilePath);
        % Write the content to the output file
        fprintf(fileID, 'Data from %s:\n%s\n\n', dbsID, fileContent);
    else
        % Print a message to the output file if the phoneme file does not exist
        fprintf(fileID, '%s not available\n\n', sprintf('%s_produced_phoneme.txt', dbsID));
        disp(fprintf('%s_produced_phoneme.txt not available', dbsID));  % Also display message in the command window
    end
end

% Close the output file
fclose(fileID);

% Confirm completion
disp('Data compilation complete. Check the output file.');