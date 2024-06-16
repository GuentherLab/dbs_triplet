% clear
setpaths_dbs_triplet(); 
set_project_specific_variables()

%% parameters 

% analysis parameters
%%%% for baseline window, use the period from -op.base_win_sec(1) to -op.base_win_sec(2) before syl #1 stim onset
%%% baseline should end at least a few 100ms before stim onset in order to not include anticipatory activity in baseline
%%% nb: in Triplet task, time between voice offset and the subsequent trial's stim onset is likely ~2.2sec, so baseline must be shorter than this duration
field_default('op','base_win_sec', [1, 0.3]); 
field_default('op','post_speech_win_sec',0.5); % time to include after 3rd-syllable voice offset in response timecourse
field_default('op','min_trials_for_good_channel', 4); 
field_default('op','responsivity_alpha', 0.05);  % consider electrodes responsive if they have above-baseline responses during one response epoch at this level

% if following option is true, divide responses by baseline after subtracting baseline
%%% this will turn response values into more interpretable units
field_default('op','divide_response_by_baseline',1); 

field_default('op','rereference_method','CTAR')

% for responses during syl 1 production, start the analyzed 'speech period' this early in seconds to capture pre-sound muscle activation
% also end the prep period this early
%%% extending the window in this way for syl 2 and syl 3 might be trickier, because there is often very little time between the preceding offset and the syl2/3 onset
prod_syl1_window_extend_start = 0;  

%% load subject trial tables containing the phonotactic probabilities
% Define the path and filename
filePath = 'Z:\DBS\Analysis\triplet_results_am\archive';
fileName = 'resp_all_subjects_20230203'; % Replace this with the actual filename if different

% Load the 'subs' variable from the file
load(fullfile(filePath, fileName), 'subs');

% Initialize arrays to store num_error_phoneme and phonotactic_probabilities
num_error_phoneme = [];
phonotactic_probabilities = [];

% Loop through each subject in the 'subs' variable
for i = 1:height(subs)
    % Get the trial table for the current subject
    trial_table = subs{i}.trial_table;
    
    % Extract num_error_phoneme and phonotactic_probabilities
    num_error_phoneme = [num_error_phoneme; trial_table.num_error_phoneme];
    phonotactic_probabilities = [phonotactic_probabilities; trial_table.phonotactic_probabilities];
end

%% load and organize data
% add the following variables to the electrodes response table... use 'electrode'/'chan' as key variable
info_vars_to_copy = {'chan','type','connector','port','strip','comment','target','side','nat_x','nat_y','nat_z',...
    'leadDBS_x','leadDBS_y','leadDBS_z','tkRAS_x','tkRAS_y','tkRAS_z','mni_linear_x','mni_linear_y','mni_linear_z',...
    'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z','fs_anatomy','fs_mni152_cvs_x','fs_mni152_cvs_y','fs_mni152_cvs_z',...
    'DISTAL_label_1','DISTAL_weight_1','DISTAL_label_2','DISTAL_weight_2','DISTAL_label_3','DISTAL_weight_3',...
    'MOREL_label_1','MOREL_weight_1','MOREL_label_2','MOREL_weight_2','MOREL_label_3','MOREL_weight_3',...
    'HCPMMP1_label_1','HCPMMP1_weight_1','HCPMMP1_label_2','HCPMMP1_weight_2'};

% specify phonemic features so that they appear in the same order across all subjects
phonunits = {'cons','vow', 'syl'}; nphonunits = length(phonunits); 
unqcons = {'gh','s','t','v'};
unqvow = {'ah','ee','oo'};
unqsyl = {'ghah','ghee','ghoo','sah','see','soo','tah','tee','too','vah','vee','voo'};

loadfile = [FT_FILE_PREFIX op.resp_signal '_trial_ar-',op.art_crit, '_ref-',op.rereference_method, op.denoise_string, '.mat']; 
load(loadfile, 'D_wavpow_trial');
load_triplet_stim_beh_timing()

%% get responses in predefined epochs
%%% 'base' = average during pre-stim baseline

% remove trials from stim trials table that do not have a corresponding fieldtrip event
%%% for the purposes of finding a fully overlapping ft trial, set trialtable's trial start as beginning of baseline and trial end as end of production plus buffer
cfg = [];
cfg.trials = trials_stim_trip; 
    cfg.trials.starts = cfg.trials.starts - op.base_win_sec(1); % adjust trial start to be beginning of baseline period
    for itrial = 1:height(trials_stim_trip)
        matchrow = trials_prod_trip.session_id==trials_stim_trip.session_id(itrial) & trials_prod_trip.trial_id==trials_stim_trip.trial_id(itrial);
        if any(matchrow) % if this trial has stim timing and prod timing, end trial at prod end plus buffer
            cfg.trials.ends(itrial) = trials_prod_trip.ends(matchrow) + op.post_speech_win_sec; 
        elseif ~any(matchrow) % if this trial has stim timing but no prod timing, add buffer to end of stim timing
            cfg.trials.ends(itrial) = cfg.trials.ends(itrial) + 1 + op.post_speech_win_sec; 
        end
    end
cfg.plot_times = 0;
[trials, trials_ft]  = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D_wavpow_trial);

% rename vars
trials.syl = [trials.stim1, trials.stim2, trials.stim3]; 
trials = removevars(trials,{'stim1','stim2','stim3'});

% vars for table construction
nchans = length(D_wavpow_trial.label);
nans_ch = nan(nchans,1); 
nans_ch2 = nan(nchans,2); 
nans_ch3 = nan(nchans,3); 
logch = false(nchans,1);
ntrials_stim = height(trials); 
nans_tr = nan(ntrials_stim,1); 
nans_tr2 = nan(ntrials_stim,2); 
nans_tr3 = nan(ntrials_stim,3); 
cel_tr = cell(ntrials_stim,1); 

% info about our trial timing analysis window
trials = [trials, table(cel_tr, nans_tr3,    nans_tr3,            nans_tr3,    nans_tr3,       nans_tr2, nans_tr2,     false(ntrials_stim,1),          nans_tr, cel_tr, cel_tr, cel_tr,...
     'VariableNames', {'times', 'stim_syl_on', 'stim_syl_off', 'prod_syl_on', 'prod_syl_off', 'trans_on', 'trans_off',     'has_speech_timing', 'ft_trial_idx','cons_constit', 'vow_constit','syl_constit' })];

% stim timing does not mark our trial boundaries, so give more precise names
%%% 'id' needs to be changed - it specifically refers to row of the trialtable that has data for stim timing; not the same row numbers as other trial tables
trials = renamevars(trials,{'duration','id'}, {'dur_stim','stimtrial_id'}); 

n_unqcons = length(unqcons);
    trials.cons_present = false(ntrials_stim, n_unqcons);
n_unqvow = length(unqvow);
    trials.vow_present = false(ntrials_stim, n_unqvow);
n_unqsyl = length(unqsyl);
    trials.syl_present = false(ntrials_stim, n_unqsyl);

% table containing responses during epochs for each chan
cel = repmat({nans_tr},nchans,1); % 1 value per trial per chan
cel2 = repmat({nans_tr2},nchans,1); % 2 values per trial per chan
cel3 = repmat({nans_tr3},nchans,1); % 3 values per trial per chan
nancons = nan(nchans, n_unqcons); 
nanvow = nan(nchans, n_unqvow); 
nansyl = nan(nchans, n_unqsyl); 
resp = table(   D_wavpow_trial.label, cel,   repmat({cel_tr},nchans,1),  cel3,    cel,    cel3,  cel2,    nans_ch, nans_ch, nans_ch,  logch,  nans_ch, nans_ch,  nans_ch3,    nans_ch3,    nans_ch3,     nans_ch3,    nans_ch3,    nans_ch3,     nancons,          nanvow,            nansyl,           nans_ch3,      nans_ch3,    nans_ch3,    nans_ch3,         nans_ch,           logch, ...
  'VariableNames', {'chan', 'base', 'timecourse',             'stim', 'prep', 'prod', 'trans', 'p_stim','p_prep', 'p_prod', 'rspv','p_rspv','p_rank','p_stim_cons','p_stim_vow','p_stim_syl','p_prep_cons','p_prep_vow','p_prep_syl','p_prep_cons_pref','p_prep_vow_pref','p_prep_syl_pref', 'p_prod_cons','p_prod_vow','p_prod_syl','prod_syl_mean','n_usable_trials', 'usable_chan'}); 

% extract epoch-related responses
%%%% trials.times{itrial} use global time coordinates
%%%% ....... start at a fixed baseline window before stim onset
%%%% ....... end at a fixed time buffer after speech offset
for itrial = 1:ntrials_stim % itrial is absolute index across sessions; does not equal "trial_id" from loaded tables
    isess = trials.session_id(itrial); 
    trial_id_in_sess = trials.trial_id(itrial); % session-relative trial number

    % check whether there is speech production timing data for this trial
    %%%...... if so, add speech timing data; if not, mark it for removal
    speechtrial_match = trials_prod_trip.session_id == isess & trials_prod_trip.trial_id == trial_id_in_sess;
    if any(speechtrial_match)
        trials.has_speech_timing(itrial) = true; 
        trials.ends(itrial) = trials.prod_trip.ends(speechtrial_match) + op.post_speech_win_sec;
        trials.duration(itrial) = trials.ends(itrial) - trials.starts(itrial); 
    elseif ~any(speechtrial_match)
        trials.has_speech_timing(itrial) = false; 
        
        % fill in blanks to maintain variable class consistency in cells
        trials.cons(itrial,1:3) = {'','',''};
        trials.vow(itrial,1:3) = {'','',''};
        trials.syl(itrial,1:3) = {'','',''};

        continue
    end

    % get syllable-specific timing
    ft_idx = trials.ft_idx(itrial); % get the trial in fieldtrip struct that corresponds to the current trial
    for isyl = 1:3
                % stim timing
        trials_stim_syl_row = find(trials_stim_syl.trial_id == trial_id_in_sess  &  trials_stim_syl.syl_id == isyl & trials_stim_syl.session_id == isess); 
        trials.stim_syl_on(itrial,isyl) = trials_stim_syl.starts(trials_stim_syl_row);  % stim syllable start time
        trials.stim_syl_off(itrial,isyl) = trials.stim_syl.ends(trials_stim_syl_row); % stim syllable end time  
        
           % prod timing
        trials_prod_syl_row = find(trials_prod_syl.trial_id == trial_id_in_sess  &  trials_prod_syl.syl_id == isyl & trials_prod_syl.session_id == isess); 
        trials.prod_syl_on(itrial,isyl) = trials.prod_syl.starts(trials_prod_syl_row);  % prod syllable start time
        trials.prod_syl_off(itrial,isyl) = trials.prod_syl.ends(trials_prod_syl_row); % prod syllable end time
    end

    % get indices within the trial-specific set of timepoints of D_D_wavpow_trial.time{ft_idx} that match our specified trial window
    match_time_inds = D_wavpow_trial.time{ft_idx} > trials.starts(itrial) & D_wavpow_trial.time{ft_idx} < trials.ends(itrial); 
    trials.times{itrial} = D_wavpow_trial.time{ft_idx}(match_time_inds); % times in this redefined trial window... still using global time coordinates

    % get trial-relative baseline time indices; window time-locked to first stim onset
    base_inds = D_wavpow_trial.time{ft_idx} > [trials.stim_syl_on(itrial,1) - op.base_win_sec(1)] & D_wavpow_trial.time{ft_idx} < [trials.stim_syl_on(itrial,1) - op.base_win_sec(2)]; 

    % baseline activity and timecourse
    for ichan = 1:nchans
        % use mean rather than nanmean, so that trials which had artifacts marked with NaNs will be excluded
        resp.base{ichan}(itrial) = mean( D_wavpow_trial.trial{ft_idx}(ichan, base_inds), 'includenan' ); % mean HG during baseline
        % get baseline-normalized trial timecourse
       resp.timecourse{ichan}{itrial} =  D_wavpow_trial.trial{ft_idx}(ichan, match_time_inds) - resp.base{ichan}(itrial); 

        if op.divide_response_by_baseline
            resp.timecourse{ichan}{itrial} = resp.timecourse{ichan}{itrial} ./ resp.base{ichan}(itrial); 
        end

    end
    
    % syllable-specific responses
    trials_phon_rowmatch = trials_phon.trial_id == trials.trial_id(itrial) & ...
                           trials_phon.session_id == trials.session_id(itrial); 
    for isyl = 1:3
        stim_syl_inds = D_wavpow_trial.time{ft_idx} > trials.stim_syl_on(itrial,isyl) & D_wavpow_trial.time{ft_idx} < trials.stim_syl_off(itrial,isyl); % time idx when syl 1 2 3 were played
        if isyl == 1 % for syl 1, start the analyzed 'speech window' early to account for pre-speech muscle activity
            prod_syl_inds = D_wavpow_trial.time{ft_idx} > [trials.prod_syl_on(itrial,isyl)-prod_syl1_window_extend_start] & D_wavpow_trial.time{ft_idx} < trials.prod_syl_off(itrial,isyl); % time idx when syl 1 was spoken
        elseif ismember(isyl, [2 3])
            prod_syl_inds = D_wavpow_trial.time{ft_idx} > trials.prod_syl_on(itrial,isyl) & D_wavpow_trial.time{ft_idx} < trials.prod_syl_off(itrial,isyl); % time idx when syl 2 3 were spoken
        end
        for ichan = 1:nchans
            resp.stim{ichan}(itrial,isyl) = mean( D_wavpow_trial.trial{ft_idx}(ichan, stim_syl_inds) ) - resp.base{ichan}(itrial); 
            resp.prod{ichan}(itrial,isyl) = mean( D_wavpow_trial.trial{ft_idx}(ichan, prod_syl_inds) ) - resp.base{ichan}(itrial); 
            
            if op.divide_response_by_baseline
                resp.stim{ichan}(itrial,isyl) = resp.stim{ichan}(itrial,isyl) ./  resp.base{ichan}(itrial);
                resp.prod{ichan}(itrial,isyl) = resp.prod{ichan}(itrial,isyl) ./  resp.base{ichan}(itrial);
            end

        end

        % phoneme info
        phon_row_cons = trials_phon_rowmatch & strcmp(trials_phon.type,'consonant') & trials_phon.syl_id==isyl;
            trials.cons{itrial,isyl} = trials_phon.stim{phon_row_cons};
        phon_row_vow = trials_phon_rowmatch & strcmp(trials_phon.type,'vowel') & trials_phon.syl_id==isyl;
            trials.vow{itrial,isyl} = trials_phon.stim{phon_row_vow};

    end
    
    % preparatory responses
    %%%% prep period inds: after stim ends and before first syllable prod onset, adjusted by prod_syl1_window_extend_start
    prep_inds = D_wavpow_trial.time{ft_idx} > trials.stim_syl_off(itrial,3) & D_wavpow_trial.time{ft_idx} < [trials.prod_syl_on(itrial,1) - prod_syl1_window_extend_start]; 
    for ichan = 1:nchans
        resp.prep{ichan}(itrial) = mean( D_wavpow_trial.trial{ft_idx}(ichan, prep_inds) , 'includenan' ) - resp.base{ichan}(itrial);
        
        if op.divide_response_by_baseline
            resp.prep{ichan}(itrial) = resp.prep{ichan}(itrial) ./  resp.base{ichan}(itrial);
        end

    end

    % inter-syllable transition period responses
    for itrans = 1:2
        % 'transition' periods start/end halfway through the syllable
        trials.trans_on(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,1:2) + trials.prod_syl_off(itrial,1:2)]; % avg start/end
        trials.trans_off(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,2:3) + trials.prod_syl_off(itrial,2:3)]; % avg start/end
        trans_inds = D_wavpow_trial.time{ft_idx} > trials.trans_on(itrial,itrans) & D_wavpow_trial.time{ft_idx} < trials.trans_off(itrial,itrans); 
       for ichan = 1:nchans
           resp.trans{ichan}(itrial,itrans) = mean( D_wavpow_trial.trial{ft_idx}(ichan, trans_inds) ) - resp.base{ichan}(itrial);

            if op.divide_response_by_baseline
                resp.trans{ichan}(itrial,itrans) = resp.trans{ichan}(itrial,itrans) ./  resp.base{ichan}(itrial);
            end

       end
    end


end

%% this section determines phonemic transitions between syllables (and phonotactic probabilities) at each trial.... code by Daphne Toglia (DT)
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

%% extract additional trial-specific stim information 
trials(:,{'has_speech_timing'}) = [];
% % % % % % % % % % % % % % % % % % % % % trials(:,{'t_stim_off','t_stim_on','has_speech_timing'}) = [];
ntrials = height(trials); 
get_unq_trialstim = @(x,trials,itrial) cell2mat(unique(trials{itrial,x}));

% mark with logicals which consonants/vowels/syllables occur on which trials
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

    % get the set of unique constituent cons/vow/syl in each trial for later encoding analysis
    trials.cons_constit{itrial} = get_unq_trialstim('cons',trials,itrial);
    trials.vow_constit{itrial} = get_unq_trialstim('vow',trials,itrial);
    trials.syl_constit{itrial} = get_unq_trialstim('syl',trials,itrial);
end


%% test for response types 
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
    

    %%%%%%%%%% % test for general task responsivity
    % above-baseline prod activity
    mean_stim_trial_resp = mean(resp.stim{ichan}(good_trials,:),2); % averaged within-trial response from 3 stim periods
    [~, resp.p_stim(ichan)] = ttest(mean_stim_trial_resp, zeroes_vec) ; 

    % above-baseline preparatory activity
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials), zeroes_vec) ; 
    
    % above-baseline prod activity
    mean_prod_trial_resp = mean(resp.prod{ichan}(good_trials,:),2); % averaged within-trial response from 3 prod periods
    [~, resp.p_prod(ichan)] = ttest(mean_prod_trial_resp, zeroes_vec) ; 

    % test for nonrandom responsiveness across baseline vs. stim vs. prep
    resp.p_rspv(ichan) = anova1([zeroes_vec, mean_stim_trial_resp, resp.prep{ichan}(good_trials), mean_prod_trial_resp],    [],    'off'); 
    resp.rspv(ichan) = resp.p_rspv(ichan) < op.responsivity_alpha; 

    %%%%%%%%% rank-order selectivity during production
    resp.p_rank(ichan) = anova1(resp.prod{ichan}(good_trials,:),[],'off');
    resp.prod_syl_mean(ichan,:) = nanmean(resp.prod{ichan}(good_trials,:)); % mean response in each prod position
    [~, resp.prod_rank_pref(ichan)] = max(resp.prod_syl_mean(ichan,:)); % which is the prod rank with highest respose
    

    %%%%%% selectivity in prep epoch for each specific phonemic feature in any of the 3 upcoming positions
     % compare trials in which this consonant was present to those in which it wasn't present
    for icons = 1:n_unqcons
        this_cons_trials = good_trials & trials.cons_present(:,icons);
        [~, resp.p_prep_cons_pref(ichan, icons)] = ttest2(resp.prep{ichan}(this_cons_trials), resp.prep{ichan}(~this_cons_trials));
    end

    % compare trials in which this vowel was present to those in which it wasn't present
    for ivow = 1:n_unqvow
        this_vow_trials = good_trials & trials.vow_present(:,ivow);
        [~, resp.p_prep_vow_pref(ichan, ivow)] = ttest2(resp.prep{ichan}(this_vow_trials), resp.prep{ichan}(~this_vow_trials));
    end

    % compare trials in which this syl was present to those in which it wasn't present
    for isyl = 1:n_unqsyl
        this_syl_trials = good_trials & trials.syl_present(:,isyl);
        [~, resp.p_prep_syl_pref(ichan, isyl)] = ttest2(resp.prep{ichan}(this_syl_trials), resp.prep{ichan}(~this_syl_trials));
    end

    %%%%% prep tuning for the unique, unordered constituent cons/vow/syl.... use unique combos of phonemic elements, lumping together repeated elements within a trial
    resp.p_prep_cons_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.cons_constit(good_trials),'off');
    resp.p_prep_vow_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.vow_constit(good_trials),'off');
    resp.p_prep_syl_constit(ichan) = anova1(resp.prep{ichan}(good_trials), trials.syl_constit(good_trials),'off');

    %%%%%%%% prod tuning regardless of position - concatenate responses to all 3 positions
    resp_stim_concat = reshape(resp.stim{ichan}(good_trials,:), 3*n_good_trials, 1); % stack stim responses
    resp_prod_concat = reshape(resp.prod{ichan}(good_trials,:), 3*n_good_trials, 1); % stack prod responses
    for thisphonunit = phonunits
        thisphonunit = thisphonunit{:}; %#ok<FXSET> 
        labels_concat = reshape(trials{good_trials,thisphonunit}, 3*n_good_trials, 1); % stack trial labels
        resp{ichan, ['p_stim_',thisphonunit,'_allpos']} = anova1(resp_stim_concat, labels_concat,'off');
        resp{ichan, ['p_prod_',thisphonunit,'_allpos']} = anova1(resp_prod_concat, labels_concat,'off');
    end


    %%%%%%%%%%%%%%% tuning to phonemic features in only one of the three positions
     for ipos = 1:3
        syl_in_this_pos = triplet_tablevar(trials,{'syl',ipos},good_trials);
        cons_in_this_pos = triplet_tablevar(trials,{'cons',ipos},good_trials);
        vow_in_this_pos = triplet_tablevar(trials,{'vow',ipos},good_trials);

         % stim-consonant using the response from only when the syllable in question is being heard
         resp.p_stim_cons(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), cons_in_this_pos,'off');

         % stim-vowel using the response from only when the syllable in question is being produced
         resp.p_stim_vow(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), vow_in_this_pos,'off');
        
        % stim-syl using the response from only when the syllable in question is being played
         resp.p_stim_syl(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos),syl_in_this_pos,'off');

        % prep-consonant with position
        resp.p_prep_cons(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), cons_in_this_pos,'off');

        % prep-vowel with position
        resp.p_prep_vow(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), vow_in_this_pos,'off');

        % prep-syl with position
        resp.p_prep_syl(ichan, ipos) = anova1(resp.prep{ichan}(good_trials), syl_in_this_pos,'off');

        % prod-consonant using the response from only when the syllable in question is being produced
         resp.p_prod_cons(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), cons_in_this_pos,'off');

        % prod-vowel using the response from only when the syllable in question is being produced
         resp.p_prod_vow(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), vow_in_this_pos,'off');

         % prod-syl using the response from only when the syllable in question is being produced
         resp.p_prod_syl(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos),syl_in_this_pos,'off');

     end

    %% transition selectivity - code by DT
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
%%    

    %%% electrode info - first matching electrode table channel
    %  for macro electrode locations, only copied their first listed location from the electrodes table into response table...
    %  ..... they generally move to deep structures over the course of trials, so this location will become less accurate for later trials
    resp.elc_info_row(ichan) = find(strcmp(elc_info.electrode , resp.chan{ichan}), 1);

end
    
resp = resp(resp.usable_chan,:); % remove channels with few/no usuable trials
elc_info_copy = renamevars(elc_info(resp.elc_info_row,:),'electrode','chan');
resp = join(resp, elc_info_copy(:,info_vars_to_copy)); % add elc_info to resp
resp = removevars(resp,{'elc_info_row','usable_chan'}); 
resp.sub = cellstr(repmat(op.sub, height(resp), 1));
resp = movevars(resp,{'sub','chan','fs_anatomy','MOREL_label_1','DISTAL_label_1'},'Before',1);
