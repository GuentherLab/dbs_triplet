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
subtable.triplet_tables = cell(nsubs_to_analyze, 1);
subtable.phonetic_ontarget_scored = false(nsubs_to_analyze, 1); 

% Define the phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'ahsah', 0.0007, ...
    'ahsee', 0.0008, ...
    'asoo', 0.0008, ...
    'ahtah', 0.0018, ...
    'atee', 0.0007, ...
    'ahtoo', 0.001, ...
    'ahvah', 0.0002, ...
    'ahvee', 0.0004, ...
    'ahvoo', 0.0001, ...
    'ahghah', 0.0001, ...
    'aghee', 0, ...
    'aghoo', 0, ...
    'eesah', 0.0003, ...
    'eesee', 0.0004, ...
    'eesoo', 0.0004, ...
    'eetah', 0.0019, ...
    'eetee', 0.0008, ...
    'eetoo', 0.0011, ...
    'eevah', 0.0008, ...
    'eevee', 0.001, ...
    'eevoo', 0.0007, ...
    'eeghah', 0.0006, ...
    'eeghee', 0.0005, ...
    'eeghoo', 0.0005, ...
    'oosah', 0, ...
    'oosee', 0.0001, ...
    'oosoo', 0.0001, ...
    'ootah', 0.0017, ...
    'ootee', 0.0006, ...
    'ootoo', 0.0009, ...
    'oovah', 0.0001, ...
    'oovee', 0.0003, ...
    'oovoo', 0, ...
    'ooghah', 0.0001, ...
    'ooghee', 0, ...
    'ooghoo', 0);

% Loop through each DBS ID
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    
    % Check if the current dbsID is in the list of specified DBS IDs
    if ismember(dbsID, dbsIDs)
        phonemeFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_phoneme.txt', dbsID));
        tripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID));

        % Load phoneme and triplet tables if they exist
        if exist(phonemeFilePath, 'file')
            subtable.phoneme_tables{i_sub} = readtable(phonemeFilePath);
            subtable.triplet_tables{i_sub} = readtable(tripletFilePath);

            if ismember('phonetic_ontarget', subtable.phoneme_tables{i_sub}.Properties.VariableNames)
                subtable.phonetic_ontarget_scored(i_sub) = true;

                % Initialize new columns in the triplet table
                subtable.triplet_tables{i_sub}.trans_1 = strings(height(subtable.triplet_tables{i_sub}), 1);
                subtable.triplet_tables{i_sub}.trans_2 = strings(height(subtable.triplet_tables{i_sub}), 1);
                subtable.triplet_tables{i_sub}.phonotactic_prob_trans_1 = NaN(height(subtable.triplet_tables{i_sub}), 1);
                subtable.triplet_tables{i_sub}.phonotactic_prob_trans_2 = NaN(height(subtable.triplet_tables{i_sub}), 1);
                subtable.triplet_tables{i_sub}.error_trans1 = zeros(height(subtable.triplet_tables{i_sub}), 1);
                subtable.triplet_tables{i_sub}.error_trans2 = zeros(height(subtable.triplet_tables{i_sub}), 1);

                % Loop through each trial in the triplet table
                for itrial = 1:height(subtable.triplet_tables{i_sub})
                    stim_string = subtable.triplet_tables{i_sub}.stim{itrial};

                    % Extract characters for trans_1 (positions 2, 3, 4, 5, 6)
                    if length(stim_string) >= 6
                        trans_1 = stim_string(2:6);
                        subtable.triplet_tables{i_sub}.trans_1(itrial) = trans_1;

                        % Assign the corresponding phonotactic probability
                        if isfield(phonotacticProbabilities, trans_1)
                            subtable.triplet_tables{i_sub}.phonotactic_prob_trans_1(itrial) = phonotacticProbabilities.(trans_1);
                        else
                            subtable.triplet_tables{i_sub}.phonotactic_prob_trans_1(itrial) = NaN;
                        end
                    end

                    % Extract characters for trans_2 based on the first consonant ("t", "v", "gh", "s") after which take the rest of the string
                    trans_2 = '';
                    consonants = {'t', 'v', 'gh', 's'};
                    for cons = consonants
                        idx = strfind(stim_string, cons{1});
                        if ~isempty(idx)
                            if strcmp(cons{1}, 'gh')
                                idx = idx(1);  % Take the first occurrence of 'gh'
                            end
                            if idx + 1 <= length(stim_string)
                                trans_2 = stim_string(idx:end);
                                break;
                            end
                        end
                    end
                    
                    subtable.triplet_tables{i_sub}.trans_2(itrial) = trans_2;

                    % Assign the corresponding phonotactic probability
                    if isfield(phonotacticProbabilities, trans_2)
                        subtable.triplet_tables{i_sub}.phonotactic_prob_trans_2(itrial) = phonotacticProbabilities.(trans_2);
                    else
                        subtable.triplet_tables{i_sub}.phonotactic_prob_trans_2(itrial) = NaN;
                    end

                    % Calculate errors for transitions (simplified)
                    match_rows = subtable.phoneme_tables{i_sub}.trial_id == subtable.triplet_tables{i_sub}.trial_id(itrial);
                    errors_syl1 = 0;
                    errors_syl2 = 0;
                    errors_syl3 = 0;

                    % Loop through each syllable (1 to 3)
                    for isyl = 1:3
                        syl_rows = match_rows & (subtable.phoneme_tables{i_sub}.syl_id == isyl);

                        consonant_error = nnz(subtable.phoneme_tables{i_sub}.phonetic_ontarget(syl_rows & strcmp(subtable.phoneme_tables{i_sub}.type, 'consonant')) == 0);
                        vowel_error = nnz(subtable.phoneme_tables{i_sub}.phonetic_ontarget(syl_rows & strcmp(subtable.phoneme_tables{i_sub}.type, 'vowel')) == 0);

                        total_errors = consonant_error + vowel_error;

                        if isyl == 1
                            errors_syl1 = total_errors;
                        elseif isyl == 2
                            errors_syl2 = total_errors;
                        elseif isyl == 3
                            errors_syl3 = total_errors;
                        end
                    end

                    % Store the calculated error values in the triplet table
                    subtable.triplet_tables{i_sub}.error_trans1(itrial) = errors_syl1 + errors_syl2;
                    subtable.triplet_tables{i_sub}.error_trans2(itrial) = errors_syl2 + errors_syl3;
                end

                % Write the processed triplet tables to separate files
                writetable(subtable.triplet_tables{i_sub}, fullfile(outputDir, sprintf('%s_triplet_table_processed.txt', dbsID)));
            end
        end
    end
end

disp('Data processing and file writing complete.');

%%
%% Correlation Analysis B/W Trans_1 + 2 errors w/ Phonotactic Probability

% Initialize arrays to store correlation results
correlationResults_trans_1 = [];
correlationResults_trans_2 = [];
pValues_trans_1 = [];
pValues_trans_2 = [];
validDBSIDs = {};

% Loop through each DBS ID
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    
    % Find the index of the current DBS ID in the subtable
    idx = find(strcmp(subtable.subject, dbsID));
    
    if isempty(idx)
        warning('DBS ID %s not found in subtable', dbsID);
        continue;
    end
    
    % Get the triplet table for this subject
    tripletTable = subtable.triplet_tables{idx};
    
    % Check if the necessary columns are present
    if istable(tripletTable) && ismember('error_trans1', tripletTable.Properties.VariableNames) && ...
       ismember('phonotactic_prob_trans_1', tripletTable.Properties.VariableNames)
        
        % Extract the relevant data
        X1 = tripletTable.phonotactic_prob_trans_1;
        Y1 = tripletTable.error_trans1;
        X2 = tripletTable.phonotactic_prob_trans_2;
        Y2 = tripletTable.error_trans2;
        
        % Remove rows where either X1 or Y1 or X2 or Y2 have NaNs
        valid_idx1 = ~isnan(X1) & ~isnan(Y1);
        valid_idx2 = ~isnan(X2) & ~isnan(Y2);
        X1 = X1(valid_idx1);
        Y1 = Y1(valid_idx1);
        X2 = X2(valid_idx2);
        Y2 = Y2(valid_idx2);
        
        % Skip if there's insufficient data
        if length(X1) < 2 || length(Y1) < 2 || length(X2) < 2 || length(Y2) < 2
            warning('Not enough valid data points for DBS ID %s', dbsID);
            continue;
        end
        
        % Calculate the correlation coefficient and p-value for trans1
        [r_trans1, p_trans1] = corr(X1, Y1, 'Type', 'Pearson');
        
        % Calculate the correlation coefficient and p-value for trans2
        [r_trans2, p_trans2] = corr(X2, Y2, 'Type', 'Pearson');
        
        % Store the results
        correlationResults_trans_1 = [correlationResults_trans_1; r_trans1];
        pValues_trans_1 = [pValues_trans_1; p_trans1];
        correlationResults_trans_2 = [correlationResults_trans_2; r_trans2];
        pValues_trans_2 = [pValues_trans_2; p_trans2];
        validDBSIDs = [validDBSIDs; {dbsID}];
        
        fprintf('DBS ID: %s, Correlation Trans1: %.2f, p-value: %.4f, Correlation Trans2: %.2f, p-value: %.4f\n', dbsID, r_trans1, p_trans1, r_trans2, p_trans2);
    else
        warning('Required columns not found in triplet table for DBS ID %s', dbsID);
    end
end

% Path to save the results
outputPath = 'Z:\DBS\Analysis\triplet_results_am\daphne_analysis\AnalysisResults_Correlation_ErrorTrans_Phonotactic_0815.txt';

% Write the results to a text file
fileID = fopen(outputPath, 'w');
fprintf(fileID, 'DBS_ID\tCorrelation_Trans1\tP_Value_Trans1\tCorrelation_Trans2\tP_Value_Trans2\n');
for i = 1:length(validDBSIDs)
    fprintf(fileID, '%s\t%.2f\t%.4f\t%.2f\t%.4f\n', validDBSIDs{i}, correlationResults_trans_1(i), pValues_trans_1(i), correlationResults_trans_2(i), pValues_trans_2(i));
end
fclose(fileID);

disp('Correlation analysis results saved successfully.');
