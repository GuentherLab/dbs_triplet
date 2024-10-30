
%Leaving some DBS ID's blank 
clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% List of specified DBS IDs
dbsIDs = {'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
          'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', ...
          'DBS3029', 'DBS3030', 'DBS3031'};

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);
subtable.phoneme_tables = cell(nsubs_to_analyze, 1);
subtable.vowel_consonant_transitions = cell(nsubs_to_analyze, 1);  % New cell for transitions

% Loop through each DBS ID
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    
    % Check if the current dbsID is in the list of specified DBS IDs
    if ismember(dbsID, dbsIDs)
        phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
       
        % Load phoneme data if it exists
        if exist(phonemeFilePath, 'file')
            subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
            
            % Initialize new cell array for vowel-consonant transitions
            vowel_consonant_data = [];

            % Get unique session and trial IDs
            unique_trials = unique(subtable.phoneme_tables{i_sub}.trial_id);
            unique_sessions = unique(subtable.phoneme_tables{i_sub}.session_id);

            % Loop through each session and trial
            for session = unique_sessions'
                for trial = unique_trials'
                    % Get relevant rows for this session and trial
                    match_rows = (subtable.phoneme_tables{i_sub}.session_id == session) & ...
                                 (subtable.phoneme_tables{i_sub}.trial_id == trial);
                    phonemes = subtable.phoneme_tables{i_sub}(match_rows, :);
                    
                    % Ensure we have at least two rows to work with
                    if height(phonemes) >= 2
                        % Find the first vowel and first consonant
                        vowel1_row = phonemes(strcmp(phonemes.type, 'vowel'), :);
                        consonant1_row = phonemes(strcmp(phonemes.type, 'consonant'), :);

                        if ~isempty(vowel1_row) && ~isempty(consonant1_row)
                            trans1 = strcat(vowel1_row.stim{1}, consonant1_row.stim{1}); % First vowel and consonant
                        else
                            trans1 = ''; % Default value if not applicable
                        end

                        % Find the second vowel and next consonant
                        vowel2_row = phonemes(strcmp(phonemes.type, 'vowel'), :);
                        consonant2_row = phonemes(strcmp(phonemes.type, 'consonant'), :);
                        
                        if height(vowel2_row) > 1 && height(consonant2_row) > 1
                            trans2 = strcat(vowel2_row.stim{2}, consonant2_row.stim{1}); % Second vowel and next consonant
                        else
                            trans2 = ''; % Default value if not applicable
                        end

                        % Store the transitions
                        vowel_consonant_data = [vowel_consonant_data; {session}, {trial}, {trans1}, {trans2}];
                    end
                end
            end
            
            % Save vowel and consonant transition data into the subtable
            subtable.vowel_consonant_transitions{i_sub} = cell2table(vowel_consonant_data, ...
                'VariableNames', {'Session_ID', 'Trial_ID', 'Trans1', 'Trans2'});
        end
    end
end

disp('Data processing complete.'); 
