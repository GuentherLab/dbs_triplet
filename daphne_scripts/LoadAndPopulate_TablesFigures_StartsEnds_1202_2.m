for i_sub = 1:nsubs_to_analyze
    % Get DBS ID and corresponding file path
    dbsID = subtable.subject{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));

    % Debugging: Print file paths
    disp(['Processing DBS ID: ', dbsID]);
    disp(['Phoneme File Path: ', phonemeFilePath]);

    % Check if the phoneme file exists
    if exist(phonemeFilePath, 'file')
        % Load the phoneme table
        phonemeTable = readtable(phonemeFilePath);

        % Debugging: Check the contents of phonemeTable
        disp(['Loaded phoneme table for DBS ID: ', dbsID]);
        disp(head(phonemeTable)); % Show the first few rows

        % Check for missing 'starts' and 'ends' and calculate if needed
        if ~all(ismember({'starts', 'ends'}, phonemeTable.Properties.VariableNames))
            phonemeTable.starts = phonemeTable.onset; % Starts = Onset
            phonemeTable.ends = phonemeTable.onset + phonemeTable.duration; % Ends = Onset + Duration
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
        % If file does not exist, display a warning
        warning(['Phoneme file not found for DBS ID: ', dbsID]);
    end
end

% Save the updated subtable
subtableFilePath = fullfile(outputDir, 'Updated_Subtable_with_PhonemeTables.txt');
writetable(subtable, subtableFilePath, 'Delimiter', '\t');
disp(['Subtable saved with updated phoneme data at ', subtableFilePath]);
