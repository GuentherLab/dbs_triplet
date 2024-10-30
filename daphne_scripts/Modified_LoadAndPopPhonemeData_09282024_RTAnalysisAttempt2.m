clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir)

% List of specified DBS IDs
dbsIDs = { 
    'DBS3001', 'DBS3002', 'DBS3003', 'DBS3004', 'DBS3005', ...
    'DBS3006', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
    'DBS3014', 'DBS3015', 'DBS3016', 'DBS3017', 'DBS3018', ...
    'DBS3019', 'DBS3020', 'DBS3021', 'DBS3022', 'DBS3023', ...
    'DBS3024', 'DBS3025', 'DBS3026', 'DBS3027', 'DBS3028', ...
    'DBS3029', 'DBS3030', 'DBS3031', 'DBS4057', 'DBS4058', ...
    'DBS4060', 'DBS4061', 'DBS4062', 'DBS4066', 'DBS4067', ...
    'DBS4068', 'DBS4069', 'DBS4070', 'DBS4071', 'DBS4072', ...
    'DBS4073', 'DBS4074', 'DBS4075', 'DBS4076', 'DBS4077', ...
    'DBS4078', 'DBS4079', 'DBS4080', 'DBS4081', 'DBS4082', ...
    'DBS4083', 'DBS4084', 'DBS4085', 'DBS4086', 'DBS4087'
};


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
                            rt_trans1 = (consonant1_row.onset(1)) - (vowel1_row.onset(1)); % Reaction time for trans1
                            onset_trans1 = vowel1_row.onset(1); % Onset for trans1
                            duration_trans1 = vowel1_row.duration(1) + consonant1_row.duration(1); % Duration for trans1
                        else
                            trans1 = ''; % Default value if not applicable
                            rt_trans1 = NaN;
                            onset_trans1 = NaN;
                            duration_trans1 = NaN;
                        end

                        % Find the second vowel and next consonant
                        vowel2_row = phonemes(strcmp(phonemes.type, 'vowel'), :);
                        consonant2_row = phonemes(strcmp(phonemes.type, 'consonant'), :);
                        
                        if height(vowel2_row) > 1 && height(consonant2_row) > 1
                            trans2 = strcat(vowel2_row.stim{2}, consonant2_row.stim{1}); % Second vowel and next consonant
                            rt_trans2 = (consonant2_row.onset(1)) - (vowel2_row.onset(2)); % Reaction time for trans2
                            onset_trans2 = vowel2_row.onset(2); % Onset for trans2
                            duration_trans2 = vowel2_row.duration(2) + consonant2_row.duration(1); % Duration for trans2
                        else
                            trans2 = ''; % Default value if not applicable
                            rt_trans2 = NaN;
                            onset_trans2 = NaN;
                            duration_trans2 = NaN;
                        end

                        % Store the transitions and corresponding duration/onset values
                        vowel_consonant_data = [vowel_consonant_data; {session}, {trial}, {trans1}, {onset_trans1}, {duration_trans1}, {rt_trans1}, {trans2}, {onset_trans2}, {duration_trans2}, {rt_trans2}];
                    end
                end
            end
            
            % Save vowel and consonant transition data into the subtable
            subtable.vowel_consonant_transitions{i_sub} = cell2table(vowel_consonant_data, ...
                'VariableNames', {'Session_ID', 'Trial_ID', 'Trans1', 'Onset_Trans1', 'Duration_Trans1', 'RT_Trans1', ...
                                  'Trans2', 'Onset_Trans2', 'Duration_Trans2', 'RT_Trans2'});
        end
    end
end

disp('Data processing complete.'); 
