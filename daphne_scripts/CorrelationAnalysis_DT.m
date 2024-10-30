% Load the subs table
filename = 'Z:\DBS\Analysis\triplet_results_am\resp_all-subjects_hg_ar-E_ref-CTAR_denoised';
load(filename, 'subs');

% List of specified DBS IDs
dbsIDs = {'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
          'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', ...
          'DBS3029', 'DBS3030', 'DBS3031'};

% Path to the directory containing the processed triplet tables
processedDataDir = 'Z:\DBS\Analysis\triplet_results_am\daphne_analysis\Outputs\Actual_Outputs';

% Loop through each DBS ID to calculate the average Phonotactic Probabilities and align with num_error_phoneme
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    
    % Find the index of the current DBS ID in the subs table
    idx = find(strcmp(subs.subject, dbsID));
    
    if isempty(idx)
        warning('DBS ID %s not found in subs table', dbsID);
        continue;
    end
    
    if isempty(subs.trials{idx}) || isempty(subs.num_error_phoneme{idx})
        warning('subs.trials{%d} or subs.num_error_phoneme{%d} is empty for subject %s', idx, idx, dbsID);
        continue;
    end
    
    % Ensure the necessary columns exist in subs.trials
    if ismember('PhonotacticProbabilities', subs.trials{idx}.Properties.VariableNames)
        % Calculate the average of the two sets of Phonotactic Probabilities
        avgProb = mean(subs.trials{idx}.PhonotacticProbabilities, 2, 'omitnan');
        
        % Add the average to a new column
        subs.trials{idx}.Avg_PhonotacticProbabilities = avgProb;
        
        % Align num_error_phoneme with trials based on row indices
        num_trials = height(subs.trials{idx});
        error_phonemes = nan(num_trials, 1);
        
        for j = 1:num_trials
            if j <= length(subs.num_error_phoneme{idx})
                error_phonemes(j) = subs.num_error_phoneme{idx}(j);
            end
        end
        
        % Add the aligned num_error_phoneme to the trials table
        subs.trials{idx}.Aligned_num_error_phoneme = error_phonemes;
    else
        warning('PhonotacticProbabilities not found for DBS ID %s', dbsID);
    end
end

% Save the updated subs table
save(filename, 'subs');

% Initialize arrays to store results
correlationResults = [];
pValues = [];
validDBSIDs = {};

% Loop through each DBS ID and perform the correlation analysis
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    
    % Find the index of the current DBS ID in the subs table
    idx = find(strcmp(subs.subject, dbsID));
    
    if isempty(idx)
        warning('DBS ID %s not found in subs table', dbsID);
        continue;
    end
    
    if isempty(subs.trials{idx}) || isempty(subs.num_error_phoneme{idx})
        warning('subs.trials{%d} or subs.num_error_phoneme{%d} is empty for subject %s', idx, idx, dbsID);
        continue;
    end
    
    % Ensure the necessary columns exist in subs.trials and subs
    if ismember('Avg_PhonotacticProbabilities', subs.trials{idx}.Properties.VariableNames) && ...
       ismember('Aligned_num_error_phoneme', subs.trials{idx}.Properties.VariableNames)
       
        X = subs.trials{idx}.Avg_PhonotacticProbabilities;
        Y = subs.trials{idx}.Aligned_num_error_phoneme;
        
        % Remove any rows with NaN values
        valid_idx = ~isnan(X) & ~isnan(Y);
        X = X(valid_idx);
        Y = Y(valid_idx);
        
        if length(X) < 2 || length(Y) < 2
            warning('Not enough valid data points for DBS ID %s', dbsID);
            continue;
        end
        
        % Calculate the correlation coefficient and p-value
        [r, p] = corr(X, Y, 'Type', 'Pearson');
        
        % Store the results
        correlationResults = [correlationResults; r];
        pValues = [pValues; p];
        validDBSIDs = [validDBSIDs; {dbsID}];
        
        fprintf('DBS ID: %s, Correlation: %.2f, p-value: %.4f\n', dbsID, r, p);
    else
        warning('Required columns not found for DBS ID %s', dbsID);
    end
end

% Path to save the results
outputPath = 'Z:\DBS\Analysis\triplet_results_am\daphne_analysis\AnalysisResultsFromCorrelation_AveragePP_0717.txt';

% Write the results to a text file
fileID = fopen(outputPath, 'w');
fprintf(fileID, 'DBS_ID\tCorrelation\tP_Value\n');
for i = 1:length(validDBSIDs)
    fprintf(fileID, '%s\t%.2f\t%.4f\n', validDBSIDs{i}, correlationResults(i), pValues(i));
end
fclose(fileID);

disp('Results saved successfully.');
%next steps: could dig more into behavioral - more granular about individual phonotactic probabilities - have three syllables on which you could have errors - what you could do is for each subject - instead of 444 trials could have 80 elements and look at error rate associated with those 
% taking 3 possible errors for syllables 1 - 2 - 3 - getting one average phonotactic probability - more accurate can have phonotactic probability per trial and have each of those associated with each error surrounding it - if i have low PP b/1 syl 1 + 2 (trans1) are you more likely to get errrors in syl 1 + 2 - instead of asking it at broad trial level - breaking it down first half of the trial, second half of the trial - higher rate of error in syllables 1 + 2 
% will require going back to look at average error rate per trial 