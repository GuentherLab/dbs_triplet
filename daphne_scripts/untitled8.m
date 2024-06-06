% Define the base directory
baseDir = 'Z:\DBS';

% List of DBS IDs to process
dbsIDs = {'DBS3002', 'DBS3004', 'DBS3006', 'DBS3010', 'DBS3011', 'DBS3014', ...
          'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3032', 'DBS4079', 'DBS4080'};

% Initialize an empty table with predefined variable names and types
trialtable = table('Size', [0 3], 'VariableTypes', {'string', 'double', 'double'}, 'VariableNames', {'DBSID', 'TrialID', 'NumErrorPhonemes'});

% Process each DBS ID
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
   
    % Check if the phoneme file exists
    if exist(phonemeFilePath, 'file')
        % Load phoneme data for the current DBS ID
        produced_phoneme = readtable(phonemeFilePath);

        % Calculate errors for each trial
        ntrials = max(produced_phoneme.trial_id);
        for itrial = 1:ntrials
            matching_rows = produced_phoneme.trial_id == itrial;
            nerrors = sum(strcmp(produced_phoneme.accuracy(matching_rows), 'inaccurate'));

            % Append results to the trial table
            trialtable = [trialtable; {dbsID, itrial, nerrors}];
        end
    else
        disp(['Phoneme file not found for ' dbsID]);
    end
end

% Save the trial table to a CSV file
outputFilePath = fullfile(baseDir, 'trialtable_produced_phoneme_table.csv');
writetable(trialtable, outputFilePath);

% Confirm completion
disp('Error analysis complete. Results saved.');