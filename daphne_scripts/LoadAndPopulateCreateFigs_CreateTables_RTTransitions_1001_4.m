
%Reaction time is defined here as stimulus offset - stimulus onset 
%corr analysis is performed here
%mean rt per trans is too
clear;

setpaths_dbs_triplet()

% Define the base directory
baseDir = 'Z:\DBS';
outputDir = [PATH_RESULTS, filesep, 'daphne_analysis', filesep, 'Outputs']; 
mkdir(outputDir);

% Define the phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'oot', 0, ...
    'oov', 0, ...
    'oogh', 0, ...
    'oos', 0, ...
    'aht', 0.0001, ...
    'ahv', 0.0001, ...
    'ahgh', 0, ...
    'ahs', 0.0007, ...
    'eet', 0.0002, ...
    'eev', 0.0007, ...
    'eegh', 0.0005, ...
    'ees', 0.0003);

% Initialize the main table
subtable = readtable([PATH_ARTIFACT, filesep, 'P08_Subjects_to_analyze.txt']);
nsubs_to_analyze = height(subtable);

% Create fields to store data for each subject
subtable.stimulus_tables = cell(nsubs_to_analyze, 1);  % For stimulus triplet data
subtable.producedtriplet_tables = cell(nsubs_to_analyze, 1); % For produced triplet data
subtable.reactiontime_tables = cell(nsubs_to_analyze, 1); % Placeholder for reaction time analysis later

% Loop through each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    
    % Paths to stimulus triplet and produced triplet text files
    stimulusTripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_stimulus_triplet.txt', dbsID)); 
    producedTripletFilePath = fullfile(baseDir, dbsID, 'Preprocessed Data', 'Sync', 'annot', sprintf('%s_produced_triplet.txt', dbsID)); 
   
    % Check if the stimulus triplet file exists and load it
    if exist(stimulusTripletFilePath, 'file')
        % Load stimulus triplet table
        subtable.stimulus_tables{i_sub} = readtable(stimulusTripletFilePath);
        disp(['Loaded stimulus triplet file for DBSID: ', dbsID]);
    else
        % If the file doesn't exist, display a warning and skip to the next subject
        warning(['No stimulus triplet file found for DBSID: ', dbsID]);
        continue;
    end
    
    % Check if the produced triplet file exists and load it
    if exist(producedTripletFilePath, 'file')
        subtable.producedtriplet_tables{i_sub} = readtable(producedTripletFilePath);
        disp(['Loaded produced triplet file for DBSID: ', dbsID]);
    else
        % If the file doesn't exist, display a warning
        warning(['Produced triplet file not found for DBSID: ', dbsID]);
        continue;
    end
    
    % Initialize a reaction time table for further analysis (add 'trans1', 'trans2', 'pp1', 'pp2' columns)
    subtable.reactiontime_tables{i_sub} = table([], [], [], [], [], [], 'VariableNames', {'trial_id', 'reaction_time', 'trans1', 'trans2', 'pp1', 'pp2'});

    % Get unique trial_ids from the stimulus triplet data
    unique_trials = unique(subtable.stimulus_tables{i_sub}.trial_id);

    % Loop through each trial to calculate reaction time for stim1 only
    for trial = 1:length(unique_trials)
        trial_id = unique_trials(trial);
        
        % Filter out the data for the current trial from both tables
        trial_stimulus = subtable.stimulus_tables{i_sub}(subtable.stimulus_tables{i_sub}.trial_id == trial_id, :);
        trial_triplet = subtable.producedtriplet_tables{i_sub}(subtable.producedtriplet_tables{i_sub}.trial_id == trial_id, :);

        % Ensure both tables have data for this trial
        if isempty(trial_stimulus) || isempty(trial_triplet)
            warning(['Missing data for trial_id: ', num2str(trial_id), ' in DBSID: ', dbsID]);
            continue;
        end

        % Determine the correct timing field ('starts' or 'onset')
        if ismember('starts', trial_triplet.Properties.VariableNames)
            timing_field = 'starts';
        elseif ismember('onset', trial_triplet.Properties.VariableNames)
            timing_field = 'onset';
        else
            warning(['No "starts" or "onset" field found for trial_id: ', num2str(trial_id), ' in DBSID: ', dbsID]);
            continue;
        end

        % Calculate reaction time as stimulus offset ('ends') minus speech onset for stim1 only
        if ismember('ends', trial_stimulus.Properties.VariableNames)
            stimulus_offset = trial_stimulus.ends(1); % The end of stimulus presentation for stim1
            speech_onset = trial_triplet.(timing_field)(1); % The onset of the first syllable (stim1)

            % Calculate reaction time
            RT = stimulus_offset - speech_onset;

            % Extract last 2 letters of stim1, first letter of stim2 for trans1
            stim1 = trial_stimulus.stim1{1}; 
            stim2 = trial_stimulus.stim2{1}; 
            last_two_stim1 = stim1(end-1:end);
            first_letter_stim2 = stim2(1);
            trans1 = [last_two_stim1, first_letter_stim2];
            % Append 'h' after 'g' in trans1
            trans1 = regexprep(trans1, 'g', 'gh');

            % Extract last 2 letters of stim2, first letter of stim3 for trans2
            stim3 = trial_stimulus.stim3{1}; 
            last_two_stim2 = stim2(end-1:end);
            first_letter_stim3 = stim3(1);
            trans2 = [last_two_stim2, first_letter_stim3];
            % Append 'h' after 'g' in trans2
            trans2 = regexprep(trans2, 'g', 'gh');

            % Get phonotactic probabilities for trans1 and trans2
            if isfield(phonotacticProbabilities, trans1)
                pp1 = phonotacticProbabilities.(trans1);
            else
                pp1 = NaN;  % Set as NaN if not found
            end

            if isfield(phonotacticProbabilities, trans2)
                pp2 = phonotacticProbabilities.(trans2);
            else
                pp2 = NaN;  % Set as NaN if not found
            end

            % Add the calculated reaction time, transitions, and phonotactic probabilities to the reaction time table
            new_entry = {trial_id, RT, trans1, trans2, pp1, pp2};  % Keep RT as numeric
            subtable.reactiontime_tables{i_sub} = [subtable.reactiontime_tables{i_sub}; new_entry];
        else
            warning(['No "ends" field found for trial_id: ', num2str(trial_id), ' in DBSID: ', dbsID]);
        end
    end
end
% Loop to save reaction time tables for each subject
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};

    % Access the reaction time table for the current subject
    reactiontime_table = subtable.reactiontime_tables{i_sub};

    % Check if the reactiontime_table is empty
    if isempty(reactiontime_table)
        warning(['Reaction time table for subject ' dbsID ' is empty. Skipping this subject.']);
        continue;
    end

    % Check if the current subject's reactiontime_table is a table or cell array
    if iscell(reactiontime_table)
        reactiontime_table = reactiontime_table{1};  % Access the first element if it's in a cell
    end

    % Now check if it's a table
    if ~istable(reactiontime_table)
        warning(['Reaction time table for subject ' dbsID ' is not a table. Skipping this subject.']);
        continue;
    end

    % Check if each column is a cell array or numeric and handle accordingly
    % For trial_id
    if iscell(reactiontime_table.trial_id)
        valid_trial_id = cellfun(@isnumeric, reactiontime_table.trial_id) & ~cellfun(@isempty, reactiontime_table.trial_id);
    else
        valid_trial_id = ~isnan(reactiontime_table.trial_id);
    end

    % For reaction_time
    if iscell(reactiontime_table.reaction_time)
        valid_reaction_time = cellfun(@isnumeric, reactiontime_table.reaction_time) & ~cellfun(@isempty, reactiontime_table.reaction_time);
    else
        valid_reaction_time = ~isnan(reactiontime_table.reaction_time);
    end

    % For pp1
    if iscell(reactiontime_table.pp1)
        valid_pp1 = cellfun(@isnumeric, reactiontime_table.pp1) & ~cellfun(@isempty, reactiontime_table.pp1);
    else
        valid_pp1 = ~isnan(reactiontime_table.pp1);
    end

    % For pp2
    if iscell(reactiontime_table.pp2)
        valid_pp2 = cellfun(@isnumeric, reactiontime_table.pp2) & ~cellfun(@isempty, reactiontime_table.pp2);
    else
        valid_pp2 = ~isnan(reactiontime_table.pp2);
    end

    % Combine validity checks for all columns
    valid_rows = valid_trial_id & valid_reaction_time & valid_pp1 & valid_pp2;

    % Filter the table to only keep valid rows
    reactiontime_table = reactiontime_table(valid_rows, :);

    % Convert cell arrays to numeric arrays if necessary
    if iscell(reactiontime_table.trial_id)
        reactiontime_table.trial_id = cell2mat(reactiontime_table.trial_id);
    end

    if iscell(reactiontime_table.reaction_time)
        reactiontime_table.reaction_time = cell2mat(reactiontime_table.reaction_time);
    end

    if iscell(reactiontime_table.pp1)
        reactiontime_table.pp1 = cell2mat(reactiontime_table.pp1);
    end

    if iscell(reactiontime_table.pp2)
        reactiontime_table.pp2 = cell2mat(reactiontime_table.pp2);
    end

    % Now save the filtered reaction time table for each subject
    writetable(reactiontime_table, fullfile(outputDir, sprintf('%s_reaction_time_table.txt', dbsID)), 'Delimiter', '\t');
end

disp('Reaction time calculation, phonotactic probabilities assignment, and table writing complete.');

%%
% Perform correlation analysis between pp1 and reaction_time, pp2 and reaction_time
significant_dbsIDs_pp1 = {};  % Store DBS IDs with significant correlation for pp1
significant_dbsIDs_pp2 = {};  % Store DBS IDs with significant correlation for pp2

for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};

    % Check if reactiontime_table is empty or not a table, and skip if so
    if isempty(reactiontime_table) || ~istable(reactiontime_table)
        continue;
    end

    % Filter valid rows for numeric fields
    valid_rows = ~isnan(reactiontime_table.reaction_time) & ~isnan(reactiontime_table.pp1) & ~isnan(reactiontime_table.pp2);

    % Get valid data
    rt_values = reactiontime_table.reaction_time(valid_rows);
    pp1_values = reactiontime_table.pp1(valid_rows);
    pp2_values = reactiontime_table.pp2(valid_rows);

    % Perform correlation analysis for pp1 and reaction_time
    if length(pp1_values) > 1  % Ensure there are enough data points
        [r_pp1, p_pp1] = corr(pp1_values, rt_values, 'Type', 'Pearson');
        disp(['Correlation between pp1 and reaction_time for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_pp1), ', p = ', num2str(p_pp1)]);
        
        % Check for significance
        if p_pp1 < 0.05
            significant_dbsIDs_pp1{end+1} = dbsID;  % Store significant DBS IDs for pp1
        end
    else
        disp(['Not enough valid data for pp1 and reaction_time correlation for DBSID: ', dbsID]);
    end

    % Perform correlation analysis for pp2 and reaction_time
    if length(pp2_values) > 1  % Ensure there are enough data points
        [r_pp2, p_pp2] = corr(pp2_values, rt_values, 'Type', 'Pearson');
        disp(['Correlation between pp2 and reaction_time for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_pp2), ', p = ', num2str(p_pp2)]);
        
        % Check for significance
        if p_pp2 < 0.05
            significant_dbsIDs_pp2{end+1} = dbsID;  % Store significant DBS IDs for pp2
        end
    else
        disp(['Not enough valid data for pp2 and reaction_time correlation for DBSID: ', dbsID]);
    end
end

% Display significant DBS IDs for pp1 and pp2
disp('DBS IDs with significant correlation between pp1 and reaction_time:');
disp(significant_dbsIDs_pp1);

disp('DBS IDs with significant correlation between pp2 and reaction_time:');
disp(significant_dbsIDs_pp2);

%% Correlation Analysis Part Looking at PP1+RT1 & PP2+RT2

% Initialize lists to store significant DBS IDs for pp1 and pp2
significant_dbsIDs_pp1 = {};  % Store DBS IDs with significant correlation for pp1
significant_dbsIDs_pp2 = {};  % Store DBS IDs with significant correlation for pp2

% Perform correlation analysis between pp1 and reaction_time, pp2 and reaction_time
for i_sub = 1:nsubs_to_analyze
    dbsID = subtable.subject{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};

    % Check if reactiontime_table is empty or not a table, and skip if so
    if isempty(reactiontime_table) || ~istable(reactiontime_table)
        continue;
    end

    % Filter valid rows for numeric fields
    valid_rows = ~isnan(reactiontime_table.reaction_time) & ~isnan(reactiontime_table.pp1) & ~isnan(reactiontime_table.pp2);

    % Get valid data
    rt_values = reactiontime_table.reaction_time(valid_rows);
    pp1_values = reactiontime_table.pp1(valid_rows);
    pp2_values = reactiontime_table.pp2(valid_rows);

    % Perform correlation analysis for pp1 and reaction_time
    if length(pp1_values) > 1  % Ensure there are enough data points
        [r_pp1, p_pp1] = corr(pp1_values, rt_values, 'Type', 'Pearson');
        disp(['Correlation between pp1 and reaction_time for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_pp1), ', p = ', num2str(p_pp1)]);
        
        % Check for significance
        if p_pp1 < 0.05
            significant_dbsIDs_pp1{end+1} = dbsID;  % Store significant DBS IDs for pp1
        end
    else
        disp(['Not enough valid data for pp1 and reaction_time correlation for DBSID: ', dbsID]);
    end

    % Perform correlation analysis for pp2 and reaction_time
    if length(pp2_values) > 1  % Ensure there are enough data points
        [r_pp2, p_pp2] = corr(pp2_values, rt_values, 'Type', 'Pearson');
        disp(['Correlation between pp2 and reaction_time for DBSID: ', dbsID, ...
            ' is r = ', num2str(r_pp2), ', p = ', num2str(p_pp2)]);
        
        % Check for significance
        if p_pp2 < 0.05
            significant_dbsIDs_pp2{end+1} = dbsID;  % Store significant DBS IDs for pp2
        end
    else
        disp(['Not enough valid data for pp2 and reaction_time correlation for DBSID: ', dbsID]);
    end
end

% Display significant DBS IDs for pp1 and pp2
disp('DBS IDs with significant correlation between pp1 and reaction_time:');
disp(significant_dbsIDs_pp1);

disp('DBS IDs with significant correlation between pp2 and reaction_time:');
disp(significant_dbsIDs_pp2);

%% Plotting Everything

% Define the transition strings and their phonotactic probabilities
phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                  'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                  'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
transitions = fieldnames(phonotacticProbabilities);  % List of transitions

% DBS IDs with significant correlation between pp1 and reaction_time
significant_dbsIDs_pp1 = {'DBS3024', 'DBS4073'};

% DBS IDs with significant correlation between pp2 and reaction_time
significant_dbsIDs_pp2 = {'DBS3024', 'DBS4073', 'DBS4084'};

% Loop through each significant DBS ID and plot transitions for pp1 and pp2
for i_sub = 1:numel(significant_dbsIDs_pp1)
    dbsID = significant_dbsIDs_pp1{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};
    
    % Plot phonotactic probabilities vs mean RT for trans1 (scatter plot with Parula color scheme and x-axis buffer)
    plot_phonotactic_vs_rt(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans1');
end

for i_sub = 1:numel(significant_dbsIDs_pp2)
    dbsID = significant_dbsIDs_pp2{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};
    
    % Plot phonotactic probabilities vs mean RT for trans2 (scatter plot with Parula color scheme and x-axis buffer)
    plot_phonotactic_vs_rt(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans2');
end

%% AM Plotting Suggestion - Log as X Axis


% Define the transition strings and their phonotactic probabilities
phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                  'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                  'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
transitions = fieldnames(phonotacticProbabilities);  % List of transitions

% DBS IDs with significant correlation between pp1 and reaction_time
significant_dbsIDs_pp1 = {'DBS3024', 'DBS4073'};

% DBS IDs with significant correlation between pp2 and reaction_time
significant_dbsIDs_pp2 = {'DBS3024', 'DBS4073', 'DBS4084'};

% Loop through each significant DBS ID and plot transitions for pp1 and pp2
for i_sub = 1:numel(significant_dbsIDs_pp1)
    dbsID = significant_dbsIDs_pp1{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};
    
    % Plot phonotactic probabilities vs mean RT for trans1 (scatter plot with Parula color scheme and x-axis buffer)
   plot_phonotactic_vs_rt_by_vowel_with_subplots(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans1'); %trying out jitter
end

for i_sub = 1:numel(significant_dbsIDs_pp2)
    dbsID = significant_dbsIDs_pp2{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};
    
    % Plot phonotactic probabilities vs mean RT for trans2 (scatter plot with Parula color scheme and x-axis buffer)
    plot_phonotactic_vs_rt_by_vowel_with_subplots(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans2');
end

%% Second Kind of Plotting 10182024

% Define the transition strings and their phonotactic probabilities
phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                  'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                  'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
transitions = fieldnames(phonotacticProbabilities);  % List of transitions

% DBS IDs with significant correlation between pp1 and reaction_time
significant_dbsIDs_pp1 = {'DBS3024', 'DBS4073'};

% DBS IDs with significant correlation between pp2 and reaction_time
significant_dbsIDs_pp2 = {'DBS3024', 'DBS4073', 'DBS4084'};

% Call the function for each DBS ID with significant correlations
for i_sub = 1:numel(significant_dbsIDs_pp1)
    dbsID = significant_dbsIDs_pp1{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans1 and trans2
    plot_phonotactic_vs_rt_per_dbs2(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans1');
end

% Similarly, for pp2 transitions:
for i_sub = 1:numel(significant_dbsIDs_pp2)
    dbsID = significant_dbsIDs_pp2{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans2
    plot_phonotactic_vs_rt_per_dbs2(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans2');
end
%%
%% Second Kind of Plotting 10182024

% Define the transition strings and their phonotactic probabilities
phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                  'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                  'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
transitions = fieldnames(phonotacticProbabilities);  % List of transitions

% DBS IDs with significant correlation between pp1 and reaction_time
significant_dbsIDs_pp1 = {'DBS3024', 'DBS4073'};

% DBS IDs with significant correlation between pp2 and reaction_time
significant_dbsIDs_pp2 = {'DBS3024', 'DBS4073', 'DBS4084'};

% Call the function for each DBS ID with significant correlations
for i_sub = 1:numel(significant_dbsIDs_pp1)
    dbsID = significant_dbsIDs_pp1{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans1 and trans2
    plot_phonotactic_vs_rt_per_dbs3(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans1');
end

% Similarly, for pp2 transitions:
for i_sub = 1:numel(significant_dbsIDs_pp2)
    dbsID = significant_dbsIDs_pp2{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans2
    plot_phonotactic_vs_rt_per_dbs3(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans2');
end

%%
%% Second Kind of Plotting 10182024

% Define the transition strings and their phonotactic probabilities
phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                  'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                  'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
transitions = fieldnames(phonotacticProbabilities);  % List of transitions

% DBS IDs with significant correlation between pp1 and reaction_time
significant_dbsIDs_pp1 = {'DBS3024', 'DBS4073'};

% DBS IDs with significant correlation between pp2 and reaction_time
significant_dbsIDs_pp2 = {'DBS3024', 'DBS4073', 'DBS4084'};

% Call the function for each DBS ID with significant correlations
for i_sub = 1:numel(significant_dbsIDs_pp1)
    dbsID = significant_dbsIDs_pp1{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans1 and trans2
    plot_phonotactic_vs_rt_per_dbs5(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans1');
end

% Similarly, for pp2 transitions:
for i_sub = 1:numel(significant_dbsIDs_pp2)
    dbsID = significant_dbsIDs_pp2{i_sub};
    reactiontime_table = subtable.reactiontime_tables{i_sub};  % Ensure this table exists for the subject
    phonotacticProbabilities = struct('ahv', 0.0001, 'ahgh', 0.0008, 'ahs', 0.0006, 'aht', 0.0007, 'oot', 0.0008, 'oov', 0.0007, 'oogh', 0.0001, 'oos', 0.0003, 'eet', 0.0001, 'eev', 0.0005, 'eegh', 0.0002, 'ees', 0.0004);
    transitions = fieldnames(phonotacticProbabilities);  % Transition names
    
    % Call the plot function for trans2
    plot_phonotactic_vs_rt_per_dbs5(reactiontime_table, phonotacticProbabilities, transitions, dbsID, 'trans2');
end

