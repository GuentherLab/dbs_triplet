clear;
% Define the base directory
baseDir = 'Z:\DBS';
% List of DBS IDs
dbsIDs = {'DBS3002', 'DBS3004', 'DBS3006', 'DBS3010', 'DBS3011', 'DBS3014', ...
          'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3032', 'DBS4079', 'DBS4080'};
% Initialize the main table
subtable = table(dbsIDs','VariableNames',{'subject'});
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.triplet_tables = cell(nsubs_to_analyze,1);

% Loop through each DBS ID
for i_sub = 1:length(dbsIDs)
    dbsID = dbsIDs{i_sub};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
    tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));
   
    % Load phoneme and triplet tables if they exist
    if exist(phonemeFilePath, 'file')
        subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
        subtable.triplet_tables{i_sub} = readtable(tripletFilePath);
       
        % Calculate errors and add to the triplet table
        ntrials = height(subtable.triplet_tables{i_sub});
        for itrial = 1:ntrials
            match_rows = subtable.phoneme_tables{i_sub}.trial_id == subtable.triplet_tables{i_sub}.trial_id(itrial);
            nerrors = nnz(strcmp(subtable.phoneme_tables{i_sub}.accuracy(match_rows),'inaccurate'));
            subtable.triplet_tables{i_sub}.num_error_phonemes(itrial) = nerrors;
        end

        % Write the processed phoneme and triplet tables to separate files
        writetable(subtable.phoneme_tables{i_sub}, fullfile(baseDir, 'Outputs', sprintf('%s_phoneme_table_processed.txt', dbsID)));
        writetable(subtable.triplet_tables{i_sub}, fullfile(baseDir, 'Outputs', sprintf('%s_triplet_table_processed.txt', dbsID)));
    end
end
disp('Data processing and file writing complete.');
