%This MATLAB script that performs all the tasks you described: reading an input file, checking for rows where the "Inaccurate" value is greater than zero, and writing the corresponding "DBS#" values to a new text file.

% Define the file paths
inputFilePath = 'Z:\DBS\accuracy_results.txt';  % Path to the input data file
outputFilePath = 'Z:\DBS\filtered_DBS_numbers.txt';  % Path for the output file

% Read the data from the input file
opts = detectImportOptions(inputFilePath, 'FileType', 'text');
accuracyResults = readtable(inputFilePath, opts);

% Display the data to ensure it's read correctly (optional)
disp(head(accuracyResults));

% Filter rows where the 'Inaccurate' column is greater than 0
filteredDBS = accuracyResults(accuracyResults.Inaccurate > 0, :);

% Display filtered DBS# for verification (optional)
disp('Filtered DBS Numbers:');
disp(filteredDBS.DBS_);

% Open the file for writing
fileID = fopen(outputFilePath, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('File could not be opened for writing.');
end

% Write a header to the new file (optional)
fprintf(fileID, 'Filtered DBS#\n');

% Write each filtered DBS# to the file
for i = 1:height(filteredDBS)
    fprintf(fileID, '%s\n', filteredDBS.DBS_{i});
end

% Close the file
fclose(fileID);

% Confirm completion
disp('Filtered DBS# have been written to the new file.');
