% Load the subs table
filename = 'Z:\DBS\Analysis\triplet_results_am\resp_all-subjects_hg_ar-E_ref-CTAR_denoised';
load(filename, 'subs');

% List of specified DBS IDs
dbsIDs = {'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
          'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', ...
          'DBS3029', 'DBS3030', 'DBS3031'};

% Initialize arrays to store data for mixed-effects model
subjectID = [];
phonotacticProbabilities1 = [];
phonotacticProbabilities2 = [];
numErrorPhonemes = [];

for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    idx = find(strcmp(subs.subject, dbsID));
    
    if ~isempty(idx)
        % Check if the necessary columns exist
        if ismember('stim_trial_id', subs.trials{idx}.Properties.VariableNames) && ...
           ismember('PhonotacticProbabilities', subs.trials{idx}.Properties.VariableNames) && ...
           ~isempty(subs.num_error_phoneme{idx})
            
            % Extract data
            trialIDs = subs.trials{idx}.stim_trial_id;
            phonProb = subs.trials{idx}.PhonotacticProbabilities;
            numErrors = subs.num_error_phoneme{idx};
            
            % Ensure same length for all arrays
            minLength = min([length(trialIDs), length(numErrors), size(phonProb, 1)]);
            trialIDs = trialIDs(1:minLength);
            phonProb = phonProb(1:minLength, :);
            numErrors = numErrors(1:minLength);

            % Check if phonProb has two columns and is numeric
            if size(phonProb, 2) == 2
                % Convert to numeric if necessary
                phonProb1 = str2double(string(phonProb(:, 1)));
                phonProb2 = str2double(string(phonProb(:, 2)));

                % Replace NaNs resulting from non-numeric conversion with zero
                phonProb1(isnan(phonProb1)) = 0;
                phonProb2(isnan(phonProb2)) = 0;

                % Concatenate data for mixed-effects model
                subjectID = [subjectID; repmat({dbsID}, length(phonProb1), 1)];
                phonotacticProbabilities1 = [phonotacticProbabilities1; phonProb1];
                phonotacticProbabilities2 = [phonotacticProbabilities2; phonProb2];
                numErrorPhonemes = [numErrorPhonemes; numErrors];
            else
                warning('PhonotacticProbabilities does not have two columns for DBS ID %s', dbsID);
            end
        else
            warning('Required columns not found for DBS ID %s', dbsID);
        end
    else
        warning('DBS ID %s not found in subs table', dbsID);
    end
end

% Ensure phonotacticProbabilities1 and phonotacticProbabilities2 have the same length
if length(phonotacticProbabilities1) ~= length(phonotacticProbabilities2)
    error('PhonotacticProbabilities columns do not have the same length.');
end

% Check if all columns are numeric
if ~isnumeric(phonotacticProbabilities1) || ~isnumeric(phonotacticProbabilities2) || ~isnumeric(numErrorPhonemes)
    error('Predictor variables must be numeric.');
end

% Create table for mixed-effects model
tbl = table(subjectID, phonotacticProbabilities1, phonotacticProbabilities2, numErrorPhonemes, ...
    'VariableNames', {'SubjectID', 'PhonProb1', 'PhonProb2', 'NumErrors'});

% Convert subjectID to categorical
tbl.SubjectID = categorical(tbl.SubjectID);

% Display the first few rows of the table for inspection
disp('Preview of the table:');
disp(head(tbl));

% Fit linear mixed-effects model
try
    lme = fitlme(tbl, 'NumErrors ~ PhonProb1 + PhonProb2 + (1|SubjectID)');
    % Display results
    disp(lme);

    % Save results to a text file
    resultsFilename = 'Z:\DBS\Analysis\triplet_results_am\daphne_analysis\AnalysisResultsFromCorrelation_AveragePP_0717.txt';
    fileID = fopen(resultsFilename, 'w');
    fprintf(fileID, 'Linear Mixed-Effects Model Results\n\n');
    fprintf(fileID, '%s\n', evalc('disp(lme)'));
    fclose(fileID);
    disp('Results saved successfully.');
catch ME
    disp('Error fitting the linear mixed-effects model:');
    disp(ME.message);
    disp('Please check the preview of the table for any issues with the data.');
end
