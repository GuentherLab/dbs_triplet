% Load the subs table
filename = 'Z:\DBS\Analysis\triplet_results_am\resp_all-subjects_hg_ar-E_ref-CTAR_denoised';
load(filename, 'subs');

% List of specified DBS IDs
dbsIDs = {'DBS3003', 'DBS3004', 'DBS3008', 'DBS3010', 'DBS3011', 'DBS3012', ...
          'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', ...
          'DBS3022', 'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', ...
          'DBS3029', 'DBS3030', 'DBS3031'};

% Initialize arrays to store the data
phonotactic_probabilities = []; %dimensions of phonotactic probabilities should be the same as num_error_phonemes 
%with an error # number should have a paired - want them to match up for the linear model - gives an option to make it subject specific 
% could do multi-level model - which basically just says - you want to take into account the fact that there are multiple subjects 
% run a pearson correlation - for loop - run through each of the subjects - within that subjects vector of phonotactic probabilities and vector of num_error_phonemes 
num_error_phonemes = [];

% Loop through each DBS ID and extract the relevant data
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    
    % Find the index of the current DBS ID in the subs table
    idx = find(strcmp(subs.subject, dbsID)); % Adjust 'subject' if the column name is different
    
    if ~isempty(idx)
        % Extract the phonotactic probabilities from subs.trials
        trials = subs.trials{idx}; % Assuming 'trials' is a column containing the trials table
        phonotactic_probabilities = [phonotactic_probabilities; trials.PhonotacticProbabilities]; % Adjust if necessary
        
        % Extract the number of error phonemes from column 22
        num_error_phonemes = [num_error_phonemes; subs{idx, 22}]; % Column 22 for num_error_phonemes
    end
end

% Convert to numeric arrays if they are not already
phonotactic_probabilities = double(phonotactic_probabilities);
num_error_phonemes = cellfun(@double,num_error_phonemes,'UniformOutput',false);

% Perform linear regression
X = phonotactic_probabilities;
y = num_error_phonemes;

% Fit the linear regression model
mdl = fitlm(X, y);

% Display the model summary
disp(mdl);

% Plot the regression line
figure;
scatter(X, y, 'filled');
hold on;
plot(mdl);
xlabel('Phonotactic Probability');
ylabel('Number of Error Phonemes');
title('Linear Regression: Phonotactic Probability vs. Number of Error Phonemes');
hold off;
