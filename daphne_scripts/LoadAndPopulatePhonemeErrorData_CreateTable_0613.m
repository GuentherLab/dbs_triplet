% Clear previous variables and data
clear

% Define the base directory
baseDir = 'Z:\DBS';

% List of DBS IDs
dbsIDs = {'DBS3001', 'DBS3002', 'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', ...
          'DBS3011', 'DBS3012', 'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', ...
          'DBS3019', 'DBS3020', 'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', ...
          'DBS3027', 'DBS3028', 'DBS3029', 'DBS3030', 'DBS3031', 'DBS3032'};

% Initialize the subtable
subtable = table(dbsIDs', 'VariableNames', {'subject'});
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);
subtable.triplet_tables = cell(nsubs_to_analyze, 1);

% Output file path
outputFilePath = fullfile(baseDir, 'trialtable_produced_phoneme_table.txt');
fileID = fopen(outputFilePath, 'w');  % Open the file for writing
if fileID == -1
    error('File could not be opened for writing. Check permissions and path.');
end

fprintf(fileID, 'Compiled Data from Various DBS Phoneme Files\n\n'); % Header

tic;
% Loop through each DBS ID
for i_sub = 1:length(dbsIDs)
    dbsID = dbsIDs{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));

    % Check if the phoneme file exists
    if exist(phonemeFilePath, 'file')
        fprintf(fileID, '\n--------------------------------------------------\n');
        fprintf(fileID, 'DBS ID: %s\n', dbsID);
        fprintf(fileID, '--------------------------------------------------\n');

        % Read and write content
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

        ntrials = height(subtable.triplet_tables{i_sub});

        for itrial = 1:ntrials 
            match_rows = subtable.phoneme_tables{i_sub}.trial_id == subtable.triplet_tables{i_sub}.trial_id(itrial);
            match_rows = match_rows & subtable.phoneme_tables{i_sub}.session_id == subtable.triplet_tables{i_sub}.session_id(itrial);
            nerrors = nnz(subtable.phoneme_tables{i_sub}.phonetic_ontarget(match_rows) ~= 1);
            subtable.triplet_tables{i_sub}.num_error_phonemes(itrial) = nerrors;
        end
    else
        fprintf(fileID, '\n--------------------------------------------------\n');
        fprintf(fileID, 'Notice: %s_produced_phoneme.txt not available\n', dbsID);
        fprintf(fileID, '--------------------------------------------------\n');
    end
end
time_elapsed = toc;

fclose(fileID); % Close the output file

disp('Data compilation complete. Check the output file.');
