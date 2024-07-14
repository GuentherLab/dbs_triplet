% Set paths and load necessary data
setpaths_dbs_triplet(); 
set_project_specific_variables();

% Define the list of DBS IDs
dbsIDs = {'DBS3002', 'DBS3004', 'DBS3010', 'DBS3011', 'DBS3014', 'DBS3017', ...
     'DBS3018', 'DBS3019', 'DBS3020', 'DBS3032'};

% Load the phonotactic data
phonotacticFile = 'Z:\DBS\Analysis\triplet_results_am\resp_all-subjects_hg_ar-E_ref-CTAR_denoised';
load(phonotacticFile, 'subs');

% Initialize analysis parameters
op.base_win_sec = [1, 0.3];
op.post_speech_win_sec = 0.5;
op.min_trials_for_good_channel = 4;
op.responsivity_alpha = 0.05;
op.divide_response_by_baseline = 1;
op.rereference_method = 'CTAR';

% Initialize more variables
phonotacticData = cell(length(dbsIDs), 1);

% Loop through each DBS ID and process data
for i = 1:length(dbsIDs)
    dbsID = dbsIDs{i};
    subjectIndex = strcmp({subs.subject}, dbsID);
    
    if any(subjectIndex)
        phonotacticData{i} = subs(subjectIndex).trials;
    else
        warning('No phonotactic data found for DBS ID: %s', dbsID);
    end
end

% Load data
loadfile = [FT_FILE_PREFIX op.resp_signal '_trial_ar-', op.art_crit, '_ref-', op.rereference_method, op.denoise_string, '.mat']; 
load(loadfile, 'D_wavpow_trial');
load_triplet_stim_beh_timing();

% Initialize variables for table construction
nchans = length(D_wavpow_trial.label);
ntrials_stim = height(trials); 
resp = table(D_wavpow_trial.label, repmat({cell(ntrials_stim,1)},nchans,1), 'VariableNames', {'chan', 'base'});

% Initialize more variables
resp.base = repmat({nan(ntrials_stim, 1)}, nchans, 1);
resp.timecourse = repmat({cell(nchans,1)}, nchans, 1);

% Add the following variables to the electrodes response table... use 'electrode'/'chan' as key variable
info_vars_to_copy = {'chan','type','connector','port','strip','comment','target','side','nat_x','nat_y','nat_z',...
    'leadDBS_x','leadDBS_y','leadDBS_z','tkRAS_x','tkRAS_y','tkRAS_z','mni_linear_x','mni_linear_y','mni_linear_z',...
    'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z','fs_anatomy','fs_mni152_cvs_x','fs_mni152_cvs_y','fs_mni152_cvs_z',...
    'DISTAL_label_1','DISTAL_weight_1','DISTAL_label_2','DISTAL_weight_2','DISTAL_label_3','DISTAL_weight_3',...
    'MOREL_label_1','MOREL_weight_1','MOREL_label_2','MOREL_weight_2','MOREL_label_3','MOREL_weight_3',...
    'HCPMMP1_label_1','HCPMMP1_weight_1','HCPMMP1_label_2','HCPMMP1_weight_2'};

% Initialize analysis parameters
field_default('op','base_win_sec', [1, 0.3]); 
field_default('op','post_speech_win_sec',0.5); 
field_default('op','min_trials_for_good_channel', 4); 
field_default('op','responsivity_alpha', 0.05);  
field_default('op','divide_response_by_baseline',1); 
field_default('op','rereference_method','CTAR')

prod_syl1_window_extend_start = 0;  

% Load and organize data
loadfile = [FT_FILE_PREFIX op.resp_signal '_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat']; 
load(loadfile, 'D_wavpow_trial');
load_triplet_stim_beh_timing()

% Get responses in predefined epochs
cfg = [];
cfg.trials = trials_stim_trip; 
cfg.trials.starts = cfg.trials.starts - op.base_win_sec(1);
for itrial = 1:height(trials_stim_trip)
    matchrow = trials_prod_trip.session_id==trials_stim_trip.session_id(itrial) & trials_prod_trip.trial_id==trials_stim_trip.trial_id(itrial);
    if any(matchrow)
        cfg.trials.ends(itrial) = trials_prod_trip.ends(matchrow) + op.post_speech_win_sec; 
    else
        cfg.trials.ends(itrial) = cfg.trials.ends(itrial) + 1 + op.post_speech_win_sec; 
    end
end
cfg.plot_times = 0;
[trials, trials_ft]  = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D_wavpow_trial);

trials.syl = [trials.stim1, trials.stim2, trials.stim3]; 
trials = removevars(trials,{'stim1','stim2','stim3'});

nchans = length(D_wavpow_trial.label);
ntrials_stim = height(trials); 
resp = table(D_wavpow_trial.label, repmat({cell(ntrials_stim,1)},nchans,1), 'VariableNames', {'chan', 'base'});

% Extract epoch-related responses
for itrial = 1:ntrials_stim % Adjust the loop index as necessary
    isess = trials.session_id(itrial); 
    trial_id_in_sess = trials.trial_id(itrial); % session-relative trial number

    % Initialize resp.base{ichan}
    for ichan = 1:nchans
        resp.base{ichan} = nan(ntrials_stim, 1); % Initialize with NaN
    end

    ft_idx = trials.ft_idx(itrial); 
    for isyl = 1:3
        stim_syl_inds = D_wavpow_trial.time{ft_idx} > trials.stim_syl_on(itrial,isyl) & D_wavpow_trial.time{ft_idx} < trials.stim_syl_off(itrial,isyl); 
        if isyl == 1
            prod_syl_inds = D_wavpow_trial.time{ft_idx} > [trials.prod_syl_on(itrial,isyl)-prod_syl1_window_extend_start] & D_wavpow_trial.time{ft_idx} < trials.prod_syl_off(itrial,isyl); 
        elseif ismember(isyl, [2 3])
            prod_syl_inds = D_wavpow_trial.time{ft_idx} > trials.prod_syl_on(itrial,isyl) & D_wavpow_trial.time{ft_idx} < trials.prod_syl_off(itrial,isyl); 
        end
        for ichan = 1:nchans
            resp.stim{ichan}(itrial,isyl) = mean( D_wavpow_trial.trial{ft_idx}(ichan, stim_syl_inds) ) - resp.base{ichan}(itrial); 
            resp.prod{ichan}(itrial,isyl) = mean( D_wavpow_trial.trial{ft_idx}(ichan, prod_syl_inds) ) - resp.base{ichan}(itrial); 
            
            if op.divide_response_by_baseline
                resp.stim{ichan}(itrial,isyl) = resp.stim{ichan}(itrial,isyl) ./  resp.base{ichan}(itrial);
                resp.prod{ichan}(itrial,isyl) = resp.prod{ichan}(itrial,isyl) ./  resp.base{ichan}(itrial);
            end
        end
    end

    % Preparatory responses
    prep_inds = D_wavpow_trial.time{ft_idx} > trials.stim_syl_off(itrial,3) & D_wavpow_trial.time{ft_idx} < [trials.prod_syl_on(itrial,1) - prod_syl1_window_extend_start]; 
    for ichan = 1:nchans
        resp.prep{ichan}(itrial) = mean( D_wavpow_trial.trial{ft_idx}(ichan, prep_inds) , 'includenan' ) - resp.base{ichan}(itrial);
        
        if op.divide_response_by_baseline
            resp.prep{ichan}(itrial) = resp.prep{ichan}(itrial) ./  resp.base{ichan}(itrial);
        end
    end

    % Inter-syllable transition period responses
    for itrans = 1:2
        trials.trans_on(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,1:2) + trials.prod_syl_off(itrial,1:2)]; 
        trials.trans_off(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,2:3) + trials.prod_syl_off(itrial,2:3)]; 
        trans_inds = D_wavpow_trial.time{ft_idx} > trials.trans_on(itrial,itrans) & D_wavpow_trial.time{ft_idx} < trials.trans_off(itrial,itrans); 
        for ichan = 1:nchans
            resp.trans{ichan}(itrial,itrans) = mean( D_wavpow_trial.trial{ft_idx}(ichan, trans_inds) ) - resp.base{ichan}(itrial);
            
            if op.divide_response_by_baseline
                resp.trans{ichan}(itrial,itrans) = resp.trans{ichan}(itrial,itrans) ./  resp.base{ichan}(itrial);
            end
        end
    end
end

% This section determines phonemic transitions between syllables (and phonotactic probabilities) at each trial.... code by Daphne Toglia (DT)
% Add new columns for transitions
trans1 = cell(height(trials), 1);
trans2 = cell(height(trials), 1);

% Loop through each row of the table
for itrial = 1:height(trials)
    % Extract the last two characters of syl1 and the first character of syl2
    syl1 = trials.syl{itrial,1};
    syl2 = trials.syl{itrial,2};
    syl3 = trials.syl{itrial,3};

    if ~isempty(syl1) % if there is phonemeic data for this trial
        % Handle the transition from syl1 to syl2
        trans1{itrial} = strcat(syl1(end-1:end), syl2(1));
    
        % Now handle the transition from syl2 to syl3
        trans2{itrial} =strcat(syl2(end-1:end),syl3(1));
    end
end

% Add new columns to the table
trials.transition_id = [trans1, trans2];

% Define the phonotactic probabilities mapping
phonotacticProbabilities = struct(...
    'oot', 0, ...
    'oov', 0, ...
    'oog', 0, ...
    'oos', 0, ...
    'aht', 0.0001, ...
    'ahv', 0.0001, ...
    'ahg', 0, ...
    'ahs', 0.0007, ...
    'eet', 0.0002, ...
    'eev', 0.0007, ...
    'eeg', 0.0005, ...
    'ees', 0.0003);

% Initialize columns for the probabilities in the trials table
probTrans1 = zeros(height(trials), 1);
probTrans2 = zeros(height(trials), 1);

% add variables for transition tuning to resp table
resp.p_phonotactic_prob = nan(nchans, 2);  % For phonotactic probability tuning
resp.p_trans_id = nan(nchans, 2);  % For transition ID tuning

% Loop through each transition and assign probabilities
for itrial = 1:height(trials)
    this_trans1 = strrep(trials.transition_id{itrial,1}, '''', ''); % Remove the single quote
    this_trans2 = strrep(trials.transition_id{itrial,2}, '''', ''); % Remove the single quote

    % Check if the transition exists in the mapping, and if so, assign its probability
    if isfield(phonotacticProbabilities, this_trans1)
        probTrans1(itrial) = phonotacticProbabilities.(this_trans1);
    else
        warning('No probability found for transition: %s', this_trans1);
    end

    if isfield(phonotacticProbabilities, this_trans2)
        probTrans2(itrial) = phonotacticProbabilities.(this_trans2);
    else
        warning('No probability found for transition: %s', this_trans2);
    end
end

% Set the format to short g for displaying fewer digits after the decimal point
format short g;

% Add the probability columns to the trials table
trials.PhonotacticProbabilities  = [probTrans1, probTrans2];

% Extract additional trial-specific stim information 
trials(:,{'has_speech_timing'}) = [];
ntrials = height(trials); 
get_unq_trialstim = @(x,trials,itrial) cell2mat(unique(trials{itrial,x}));

% Mark with logicals which consonants/vowels/syllables occur on which trials
for itrial = 1:ntrials
    for icons = 1:n_unqcons
        this_cons = unqcons{icons};
        if any(strcmp(this_cons, trials.cons(itrial,:)))
            trials.cons_present(itrial,icons) = true;
        end
    end
    for ivow = 1:n_unqvow
        this_vow = unqvow{ivow};
        if any(strcmp(this_vow, trials.vow(itrial,:)))
            trials.vow_present(itrial,ivow) = true;
        end
    end
    for isyl = 1:n_unqsyl
        this_syl = unqsyl{isyl};
        if any(strcmp(this_syl, trials.syl(itrial,:)))
            trials.syl_present(itrial,isyl) = true;
        end
    end

    % Get the set of unique constituent cons/vow/syl in each trial for later encoding analysis
    trials.cons_constit{itrial} = get_unq_trialstim('cons',trials,itrial);
    trials.vow_constit{itrial} = get_unq_trialstim('vow',trials,itrial);
    trials.syl_constit{itrial} = get_unq_trialstim('syl',trials,itrial);
end

% Test for response types 
resp.elc_info_row = nan(nchans,1); 
for ichan = 1:nchans
    good_trials = ~isnan(resp.base{ichan}); % non-artifactual trials for this channel
    n_good_trials = nnz(good_trials); 
    zeroes_vec = zeros(n_good_trials,1); 
    resp.n_usable_trials(ichan) = nnz(good_trials); 
    if resp.n_usable_trials(ichan) < op.min_trials_for_good_channel
        resp.usable_chan(ichan) = false; 
        continue; % skip stats analysis if channel had too few good trials
    else
        resp.usable_chan(ichan) = true; 
    end 

    % Test for general task responsivity
    % Above-baseline prod activity
    mean_stim_trial_resp = mean(resp.stim{ichan}(good_trials,:),2); % averaged within-trial response from 3 stim periods
    [~, resp.p_stim(ichan)] = ttest(mean_stim_trial_resp, zeroes_vec); 

    % Above-baseline preparatory activity
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials), zeroes_vec); 
    
    % Above-baseline prod activity
    mean_prod_trial_resp = mean(resp.prod{ichan}(good_trials,:),2); % averaged within-trial response from 3 prod periods
    [~, resp.p_prod(ichan)] = ttest(mean_prod_trial_resp, zeroes_vec); 

    % Test for nonrandom responsiveness across baseline vs. stim vs. prep
    resp.p_rspv(ichan) = anova1([zeroes_vec, mean_stim_trial_resp, resp.prep{ichan}(good_trials), mean_prod_trial_resp], [], 'off'); 
    resp.rspv(ichan) = resp.p_rspv(ichan) < op.responsivity_alpha; 

    % Rank-order selectivity during production
    resp.p_rank(ichan) = anova1(resp.prod{ichan}(good_trials,:),[],'off');
    resp.prod_syl_mean(ichan,:) = nanmean(resp.prod{ichan}(good_trials,:)); % mean response in each prod position
    [~, resp.prod_rank_pref(ichan)] = max(resp.prod_syl_mean(ichan,:)); % which is the prod rank with highest response

    % Selectivity in prep epoch for each specific phonemic feature in any of the 3 upcoming positions
    % Compare trials in which this consonant was present to those in which it wasn't present
    for icons = 1:n_unqcons
        this_cons_trials = good_trials & trials.cons_present(:,icons);
        [~, resp.p_prep_cons_pref(ichan, icons)] = ttest2(resp.prep{ichan}(this_cons_trials), resp.prep{ichan}(~this_cons_trials));
    end

    % Compare trials in which this vowel was present to those in which it wasn't present
    for ivow = 1:n_unqvow
        this_vow_trials = good_trials & trials.vow_present(:,ivow);
        [~, resp.p_prep_vow_pref(ichan, ivow)] = ttest2(resp.prep{ichan}(this_vow_trials), resp.prep{ichan}(~this_vow_trials));
    end

    % Compare trials in which this syl was present to those in which it wasn't present
    for isyl = 1:n_unqsyl
        this_syl_trials = good_trials & trials.syl_present(:,isyl);
        [~, resp.p_prep_syl_pref(ichan, isyl)] = ttest2(resp.prep{ichan}(this_syl_trials), resp.prep{ichan}(~this_syl_trials));
    end

    % Prep tuning for the unique, unordered constituent cons/vow/syl.... use unique combos of phonemic elements, lumping together repeated elements within a trial
    resp.p_prep_cons_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.cons_constit(good_trials),'off');
    resp.p_prep_vow_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.vow_constit(good_trials),'off');
    resp.p_prep_syl_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.syl_constit(good_trials),'off');

    % Prod tuning regardless of position - concatenate responses to all 3 positions
    resp_stim_concat = reshape(resp.stim{ichan}(good_trials,:), 3*n_good_trials, 1); % stack stim responses
    resp_prod_concat = reshape(resp.prod{ichan}(good_trials,:), 3*n_good_trials, 1); % stack prod responses
    for thisphonunit = phonunits
        thisphonunit = thisphonunit{:};
        labels_concat = reshape(trials{good_trials,thisphonunit}, 3*n_good_trials, 1); % stack trial labels
        resp{ichan, ['p_stim_',thisphonunit,'_allpos']} = anova1(resp_stim_concat, labels_concat,'off');
        resp{ichan, ['p_prod_',thisphonunit,'_allpos']} = anova1(resp_prod_concat, labels_concat,'off');
    end

    % Tuning to phonemic features in only one of the three positions
    for ipos = 1:3
        syl_in_this_pos = triplet_tablevar(trials,{'syl',ipos},good_trials);
        cons_in_this_pos = triplet_tablevar(trials,{'cons',ipos},good_trials);
        vow_in_this_pos = triplet_tablevar(trials,{'vow',ipos},good_trials);

        % Stim-consonant using the response from only when the syllable in question is being heard
        resp.p_stim_cons(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), cons_in_this_pos,'off');

        % Stim-vowel using the response from only when the syllable in question is being produced
        resp.p_stim_vow(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), vow_in_this_pos,'off');
        
        % Stim-syl using the response from only when the syllable in question is being played
        resp.p_stim_syl(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos),syl_in_this_pos,'off');

        % Prep-consonant with position
        resp.p_prep_cons(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), cons_in_this_pos,'off');

        % Prep-vowel with position
        resp.p_prep_vow(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), vow_in_this_pos,'off');

        % Prep-syl with position
        resp.p_prep_syl(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), syl_in_this_pos,'off');

        % Prod-consonant using the response from only when the syllable in question is being produced
        resp.p_prod_cons(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), cons_in_this_pos,'off');

        % Prod-vowel using the response from only when the syllable in question is being produced
        resp.p_prod_vow(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), vow_in_this_pos,'off');

        % Prod-syl using the response from only when the syllable in question is being produced
        resp.p_prod_syl(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos),syl_in_this_pos,'off');
    end

    % Transition selectivity - code by DT
    % Phonotactic Probability Tuning
    for itrans = 1:2
        valid_trials_for_phonotactic = good_trials & ~isnan(trials.PhonotacticProbabilities(:, itrans));
        valid_responses_phonotactic = resp.trans{ichan}(valid_trials_for_phonotactic);
        valid_phonotactic_probabilities = trials.PhonotacticProbabilities(valid_trials_for_phonotactic, itrans);

        % Perform ANOVA if there are valid data for phonotactic probabilities
        if ~isempty(valid_responses_phonotactic) && ~isempty(valid_phonotactic_probabilities)
            resp.p_phonotactic_prob(ichan, itrans) = anova1(valid_responses_phonotactic, valid_phonotactic_probabilities, 'off');
        else
            disp(['Insufficient data for phonotactic probability ANOVA at ichan=', num2str(ichan), ', itrans=', num2str(itrans)]);
        end
    end

    % Transition ID Tuning
    for itrans = 1:2
        % Check if transition_id is a table variable and access it correctly
        if istable(trials) && any(strcmp(trials.Properties.VariableNames, 'transition_id'))
            valid_trials_for_transition = good_trials & ~cellfun(@isempty, trials.transition_id(:, itrans));
            valid_responses_transition = resp.trans{ichan}(valid_trials_for_transition);
            valid_transition_IDs = trials.transition_id(valid_trials_for_transition, itrans);

            % Perform ANOVA if there are valid data for transition IDs, treating them as categorical
            if ~isempty(valid_responses_transition) && ~isempty(valid_transition_IDs)
                resp.p_trans_id(ichan, itrans) = anova1(valid_responses_transition, valid_transition_IDs, 'off');
            else
                disp(['Insufficient data for transition ID ANOVA at ichan=', num2str(ichan), ', itrans=', num2str(itrans)]);
            end
        else
            error('The "transition_id" variable is not found in the "trials" table.');
        end 
    end

    % Electrode info - first matching electrode table channel
    resp.elc_info_row(ichan) = find(strcmp(elc_info.electrode , resp.chan{ichan}), 1);
end

resp = resp(resp.usable_chan,:); % remove channels with few/no usable trials
elc_info_copy = renamevars(elc_info(resp.elc_info_row,:),'electrode','chan');
resp = join(resp, elc_info_copy(:,info_vars_to_copy)); % add elc_info to resp
resp = removevars(resp,{'elc_info_row','usable_chan'}); 
resp.sub = cellstr(repmat(op.sub, height(resp), 1));
resp = movevars(resp,{'sub','chan','fs_anatomy','MOREL_label_1','DISTAL_label_1'},'Before',1);
