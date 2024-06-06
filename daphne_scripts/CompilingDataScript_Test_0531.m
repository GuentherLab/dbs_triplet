% Loop through each subject in the subtable
for i = 1:height(subtable)
    subjectID = subtable.subject{i};
    phonemeData = subtable.phoneme_tables{i};
    tripletData = subtable.triplet_tables{i};

    % Define the filename for output
    filename = fullfile(baseDir, sprintf('%s_CompiledDataForTTA_DT.txt', subjectID));

    % Create a temporary file to write details
    fileID = fopen(filename, 'w');
    if fileID == -1
        error('File could not be opened for writing. Check permissions and path.');
    end
   
    % Write subject ID
    fprintf(fileID, 'Subject: %s\n\n', subjectID);
    fclose(fileID); % Close to ensure data is written

    % Write phoneme table data
    writetable(phonemeData, filename, 'FileType', 'text', 'WriteMode', 'append', 'WriteVariableNames', true, 'Delimiter', '\t');
   
    % Append triplet table data
    % Since MATLAB doesn't support direct append using writetable, we need to manually append
    % Open the file again to append
    fileID = fopen(filename, 'a');
    fprintf(fileID, '\nTriplet Table Data:\n');
    fclose(fileID);
    writetable(tripletData, filename, 'FileType', 'text', 'WriteMode', 'append', 'WriteVariableNames', true, 'Delimiter', '\t');
end

disp('Data compilation complete. Check the output files.');
