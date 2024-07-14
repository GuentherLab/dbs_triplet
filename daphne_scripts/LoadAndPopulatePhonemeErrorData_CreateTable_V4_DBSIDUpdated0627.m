%This V4 and modified script contains checks and corrections since V3 fprintf was not working..
%Recommendation = tablevariables, can nest them can do something like - need to rename variables to the ones i used before 
% dont use forloops
% now we have a metric for difficulty of trials - can regress based on phonotactic 
% next steps = answer question = does phonotactic probability correlate with errors on a trial by trial basis - errors now will have to compare where we have PP the resp_table

clear

% Define the base directory
baseDir = 'Z:\DBS';

% List of DBS IDs
dbsIDs = {'DBS3002', 'DBS3004', 'DBS3010', 'DBS3011', 'DBS3014', 'DBS3017', ...
     'DBS3018', 'DBS3019', 'DBS3020', 'DBS3032'};

%SubTable
subtable = table(dbsIDs','VariableNames',{'subject'});
nsubs_to_analyze= height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze,1);
subtable.triplet_tables = cell(nsubs_to_analyze,1);

%Idea is now that load phoneme table, put in cell in second column

% Output file path
outputFilePath = fullfile(baseDir, 'trialtable_produced_phoneme_table.txt');
fileID = fopen(outputFilePath, 'w');  % Open the file for writing
if fileID == -1
    error('File could not be opened for writing. Check permissions and path.');
end

fprintf(fileID, 'Compiled Data from Various DBS Phoneme Files\n\n'); % Header

tic;
% Loop through each DBS ID %get the accuracies 
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
%         fprintf(fileID, '%s\n', fileContent);
    else
%        % fprintf(fileID, '\n--------------------------------------------------\n');
%        % fprintf(fileID, 'Notice: %s_produced_phoneme.txt not available\n', dbsID);
%         fprintf(fileID, '--------------------------------------------------\n');
    end
           
            ntrials = height(subtable.triplet_tables{i_sub});

    for itrial = 1:ntrials % 
                        match_rows = subtable.phoneme_tables{i_sub}.trial_id == subtable.triplet_tables{i_sub}.trial_id(itrial); % logic = subtable.everything = the name of the trial we're on & matching the indices - phoneme tables that have the same ID
                       nerrors = nnz(strcmp(subtable.phoneme_tables{i_sub}.accuracy(match_rows),'inaccurate'));

%                         for iphoneme = 1:length(match_rows) % try to write this in a single or three lines - get the matching rows in the phoneme table, tell me how many of these are errors, save that onto 
%                         thisrow = match_rows(iphoneme);
%                         
%                         nerrors = 0;
%                   
%                         
%                                  if strcmp(subtable.phoneme_tables{i_sub}.accuracy(thisrow), 'inaccurate')
%                                     nerrors = errors+1;
%                                  end
%                    
%                     
%                          end
                    
                   subtable.triplet_tables{i_sub}.num_error_phonemes(itrial) = nerrors; % {} indexes into cells, () index into an array - get used to this 
                    
   end
end
time_elapsed = toc;

% 
% fclose(fileID); % Close the output file
% 
% disp('Data compilation complete. Check the output file.');
