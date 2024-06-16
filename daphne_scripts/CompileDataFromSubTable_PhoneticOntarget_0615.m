% triplet_tables only - comparing the phonotactic frequency  

% Initialize the Excel filename
excelFileName = 'CompiledDataFromSubTable.xlsx';

% Get the number of subjects
numSubjects = size(subtable, 1);

% Loop through each subject
for i = 1:numSubjects
    % Get the subject ID (assuming it is in the first column)
    subject = subtable.subject{i};
    
    % Extract the phonemes_table and triplet_tables for the current subject
    phonemeTables = subtable.phoneme_tables{i};
    tripletsTable = subtable.triplet_tables{i};
    
    % Convert tables to cell arrays for writing to Excel
    phonemesCell = [phonemeTables.Properties.VariableNames; table2cell(phonemeTables)];
    tripletsCell = [tripletsTable.Properties.VariableNames; table2cell(tripletsTable)];
    
    % Write the phonemes table to the Excel file
    sheetNamePhonemes = sprintf('%s_phonemes', subject);
    writecell(phonemesCell, excelFileName, 'Sheet', sheetNamePhonemes);
    
    % Write the triplets table to the Excel file
    sheetNameTriplets = sprintf('%s_triplets', subject);
    writecell(tripletsCell, excelFileName, 'Sheet', sheetNameTriplets);
end

fprintf('Data has been compiled and saved to %s\n', excelFileName);
