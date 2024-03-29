% look for elcs with particualr response profiles

% % % % % data-loading parameters
vardefault('SUBJECT','DBS3012');

DATE=datestr(now,'yyyymmdd');

setpaths_dbs_triplet()
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_PREPROCESSED = [PATH_SUBJECT filesep 'Preprocessed Data'];
PATH_FIELDTRIP = [PATH_PREPROCESSED filesep 'FieldTrip']; 
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];
ARTIFACT_CRIT = 'E'; 
SAMPLE_RATE = 100; % downsample rate in hz for high gamma traces

% analysis parameters
%%%% for baseline window, use the period from -base_win_sec(1) to -base_win_sec(2) before syl #1 stim onset
%%% baseline should end at least a few 100ms before stim onset in order to not include anticipatory activity in baseline
%%% nb: in Triplet task, time between voice offset and the subsequent trial's stim onset is likely ~2.2sec, so baseline must be shorter than this duration
base_win_sec = [1, 0.3]; 
post_speech_win_sec = 0.5; % time to include after 3rd-syllable voice offset in HG timecourse
min_trials_for_good_channel = 4; 

% add the following variables to the electrodes response table... use 'electrode'/'chan' as key variable
info_vars_to_copy = {'chan','type','connector','port','strip','comment','target','side','nat_x','nat_y','nat_z',...
    'leadDBS_x','leadDBS_y','leadDBS_z','tkRAS_x','tkRAS_y','tkRAS_z','mni_linear_x','mni_linear_y','mni_linear_z',...
    'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z','fs_anatomy','fs_mni152_cvs_x','fs_mni152_cvs_y','fs_mni152_cvs_z',...
	'DISTAL_label_1','DISTAL_weight_1','DISTAL_label_2','DISTAL_weight_2','DISTAL_label_3','DISTAL_weight_3',...
    'MOREL_label_1','MOREL_weight_1','MOREL_label_2','MOREL_weight_2','MOREL_label_3','MOREL_weight_3',...
    'HCPMMP1_label_1','HCPMMP1_weight_1','HCPMMP1_label_2','HCPMMP1_weight_2'};

% specify phonemic features so that they appear in the same order across all subjects
unqcons = {'gh','s','t','v'};
unqvow = {'ah','ee','oo'};
unqsyl = {'ghah','ghee','ghoo','sah','see','soo','tah','tee','too','vah','vee','voo'};

%% load data 
cd(PATH_SYNC)

load([PATH_FIELDTRIP filesep SUBJECT '_ft_hg_trial_ref_criteria_' ARTIFACT_CRIT '_denoised.mat']);

trials_stim_trip = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_triplet.txt';]); % stim timing info
trials_stim_syl = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_syllable.txt';]); % stim timing info
trials_trip = readtable([PATH_ANNOT filesep SUBJECT '_produced_triplet.txt';]); % speech timing info
trials_syl = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable.txt';]); % speech timing info
trials_phon = readtable([PATH_ANNOT filesep SUBJECT '_produced_phoneme.txt';]); % speech timing info

elc_info = readtable([PATH_ANNOT filesep SUBJECT '_electrode.txt';]); 

% some subjects don't have _stimulus_syllable.txt (3001 and 3002)- for these, maybe try to derive it from _stimulus_triplet.txt and expected durations? 

%%
% Part 1

% Define your CVC combinations array
vc_combinations = ["vahtoosee", "veetahsoo", "sooveetah", "tooveeghah", "gheeghahvah", ...
    "soovahghee", "seetahvoo", "seegheetoo", "tahsooghee", "tahseeghoo"];

% Loop through each CVC combination
for i = 1:length(vc_combinations)
    vc = char(vc_combinations(i)); % Convert to character array

    % Split the CVC combination into syllables
    syllable1 = vc(1:3); % First 3 characters
    syllable2 = vc(4:6); % Middle 3 characters
    syllable3 = vc(7:end); % Last 3 characters

    % Handle the 'gh' case by including three letters if present in the first syllable
    if contains(syllable1, "gh")
        transition1 = strcat(syllable1(end-2:end), syllable2(1)); % Adjusted for 'gh'
    else
        transition1 = strcat(syllable1(end-1:end), syllable2(1)); % Standard transition
    end

    % Handle the 'gh' case by including three letters if present in the second syllable
    if contains(syllable2, "gh")
        transition2 = strcat(syllable2(end-2:end), syllable3(1)); % Adjusted for 'gh'
    else
        transition2 = strcat(syllable2(end-1:end), syllable3(1)); % Standard transition
    end

    % Store the transitions
    transitions1{i} = transition1;
    transitions2{i} = transition2;
end

% Filtering out the same transitions before part 2
uniqueTransitions1 = unique(transitions1);
uniqueTransitions2 = unique(transitions2);

%Combining the transitions into one variable, filtering repeated transitions
CombinedTransitions = union(uniqueTransitions1,uniqueTransitions2)

%Display combined transitions
disp(CombinedTransitions)

%Modifying variables before part 2, removing 'egha' since it is not needed
% Check if the 7th element exists
if numel(CombinedTransitions) >= 7
    % Create a logical index for all elements except the 7th
    indexToRemove = true(1, numel(CombinedTransitions)); % Initialize all as true
    indexToRemove(7) = false; % Set the 7th element to false (to be removed)

    % Use the logical index to keep all elements except the 7th
    ActualTransitions = CombinedTransitions(indexToRemove);

    % Display ActualTransitions
    disp('Actual Transitions:');
    disp(ActualTransitions);
else
    % If the 7th element does not exist, copy the original array
    ActualTransitions = CombinedTransitions;
    disp('CombinedTransitions does not have a 7th element. No changes made.');
end

%More filtering before part 2, removing 'ghee' now, since it is a cv and not a vc.

% Find the index of 'ghee' in the ActualTransitions array
indexGhee = find(strcmp(ActualTransitions, 'ghee'));

% If 'ghee' is found, create a logical index to exclude it
if ~isempty(indexGhee)
    indexToKeep = true(1, numel(ActualTransitions)); % Initialize all as true
    indexToKeep(indexGhee) = false; % Set the index of 'ghee' to false (to be removed)

    % Use the logical index to keep all elements except 'ghee'
    FinalTransitions = ActualTransitions(indexToKeep);

    % Display FinalTransitions
    disp('Final Transitions after removing ''ghee'':');
    disp(FinalTransitions);
else
    % If 'ghee' is not found, no changes are made
    FinalTransitions = ActualTransitions;
    disp('''ghee'' not found in ActualTransitions. No changes made.');
end


%%
%Part 2: Assinging Phonotactic Probabilities




%% get responses in predefined epochs
% 'base' = average during pre-stim baseline
% all response values except 'base' are baseline-normalized by dividing by that trial's baseline average
ntrials_stim = height(trials); 
nchans = length(D_hg.label);
nans_ch = nan(nchans,1); 
nans_ch2 = nan(nchans,2); 
nans_ch3 = nan(nchans,3); 
logch = false(nchans,1);
nans_tr = nan(ntrials_stim,1); 
nans_tr2 = nan(ntrials_stim,2); 
nans_tr3 = nan(ntrials_stim,3); 
cel_tr = cell(ntrials_stim,1); 

% info about our trial timing analysis window


trials = [trials, table(cel_tr, nans_tr3,    nans_tr3,            nans_tr3,    nans_tr3,       nans_tr2, nans_tr2,     false(ntrials_stim,1),          nans_tr, cel_tr, cel_tr, cel_tr,...
     'VariableNames', {'times', 'stim_syl_on', 'stim_syl_off', 'prod_syl_on', 'prod_syl_off', 'trans_on', 'trans_off',     'has_speech_timing', 'ft_trial_idx','cons_constit', 'vow_constit','syl_constit' })];

% stim timing does not mark our trial boundaries, so give more precise names
%%% 'id' needs to be changed it specifically refers to row of the trialtable that has data for stim timing; not the same row numbers as other trial tables
trials = renamevars(trials,{'starts','ends','duration','id'}, {'t_stim_on','t_stim_off','dur_stim','stimtrial_id'}); 
trials.starts = trials.t_stim_on - base_win_sec(1); % trial starts at beginning of baseline window

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
resp = table(   D_hg.label, cel,   repmat({cel_tr},nchans,1),  cel3,    cel,    cel3,  cel2,    nans_ch,  nans_ch,  nans_ch3,    nans_ch3,    nans_ch3,     nans_ch3,    nans_ch3,    nans_ch3,     nancons,          nanvow,            nansyl,           nans_ch3,      nans_ch3,    nans_ch3,    nans_ch,           logch, ...
  'VariableNames', {'chan', 'base', 'timecourse',             'stim', 'prep', 'prod', 'trans', 'p_prep', 'p_rank','p_stim_cons','p_stim_vow','p_stim_syl','p_prep_cons','p_prep_vow','p_prep_syl','p_prep_cons_pref','p_prep_vow_pref','p_prep_syl_pref', 'p_prod_cons','p_prod_vow','p_prod_syl','n_usable_trials', 'usable_chan'}); 

% extract epoch-related responses
%%%% trials.times{itrial} use global time coordinates
%%%% ....... start at a fixed baseline window before stim onset
%%%% ....... end at a fixed time buffer after speech offset
for itrial = 1:ntrials_stim % itrial is absolute index across sessions; does not equal "trial_id" from loaded tables
    isess = trials.session_id(itrial); 
    trial_id_in_sess = trials.trial_id(itrial); % session-relative trial number

    % check whether there is speech production timing data for this trial
    %%%...... if so, add speech timing data; if not, mark it for removal
    speechtrial_match = trials_trip.session_id == isess & trials_trip.trial_id == trial_id_in_sess;
    if any(speechtrial_match)
        trials.has_speech_timing(itrial) = true; 
        trials.ends(itrial) = trials_trip.ends(speechtrial_match) + post_speech_win_sec;
        trials.duration(itrial) = trials.ends(itrial) - trials.starts(itrial); 
    elseif ~any(speechtrial_match)
        trials.has_speech_timing(itrial) = false; 
        
        % fill in blanks to maintain variable class consistency in cells
        trials.cons(itrial,1:3) = {'','',''};
        trials.vow(itrial,1:3) = {'','',''};
        trials.syl(itrial,1:3) = {'','',''};

        continue
    end

    % get indices within the trial-specific set of timepoints of D_hg.time{ft_idx} that match our specified trial window
    ft_idx = trials.ft_idx(itrial); % get the trial in fieldtrip struct that corresponds to the curret trial
    match_time_inds = D_hg.time{ft_idx} > trials.starts(itrial) & D_hg.time{ft_idx} < trials.ends(itrial); 
    trials.times{itrial} = D_hg.time{ft_idx}(match_time_inds); % times in this redefined trial window... still using global time coordinates
    % get trial-relative baseline time indices; window time-locked to first stim onset
    base_inds = D_hg.time{ft_idx} > trials.starts(itrial) & D_hg.time{ft_idx} < [trials.t_stim_on(itrial) - base_win_sec(2)]; 

    % baseline activity and timecourse
    for ichan = 1:nchans
        % use mean rather than nanmean, so that trials which had artifacts marked with NaNs will be excluded
        resp.base{ichan}(itrial) = mean( D_hg.trial{ft_idx}(ichan, base_inds), 'includenan' ); % mean HG during baseline
        % get baseline-normalized trial timecourse
       resp.timecourse{ichan}{itrial} =  D_hg.trial{ft_idx}(ichan, match_time_inds) - resp.base{ichan}(itrial); 
    end
    
    % syllable-specific responses and timing
    trials_phon_rowmatch = trials_phon.trial_id == trials.trial_id(itrial) & ...
                           trials_phon.session_id == trials.session_id(itrial); 
    for isyl = 1:3
            % prod timing
        trials_prod_syl_row = find(trials_syl.trial_id == trial_id_in_sess  &  trials_syl.syl_id == isyl & trials_syl.session_id == isess); 
        trials.prod_syl_on(itrial,isyl) = trials_syl.starts(trials_prod_syl_row);  % syllable start time
        trials.prod_syl_off(itrial,isyl) = trials_syl.ends(trials_prod_syl_row); % syllable end time
                % stim timing
        trials_stim_syl_row = find(trials_stim_syl.trial_id == trial_id_in_sess  &  trials_stim_syl.syl_id == isyl & trials_stim_syl.session_id == isess); 
        trials.stim_syl_on(itrial,isyl) = trials_stim_syl.starts(trials_stim_syl_row);  % syllable start time
        trials.stim_syl_off(itrial,isyl) = trials_stim_syl.ends(trials_stim_syl_row); % syllable end time   
            % time indices of epochs
        stim_syl_inds = D_hg.time{ft_idx} > trials.stim_syl_on(itrial,isyl) & D_hg.time{ft_idx} < trials.stim_syl_off(itrial,isyl); 
        prod_syl_inds = D_hg.time{ft_idx} > trials.prod_syl_on(itrial,isyl) & D_hg.time{ft_idx} < trials.prod_syl_off(itrial,isyl); 
        for ichan = 1:nchans
            resp.stim{ichan}(itrial,isyl) = mean( D_hg.trial{ft_idx}(ichan, stim_syl_inds) ) ; % times when syl 1 2 3 were played
            %             resp.prod{ichan}(itrial,isyl) = mean( D_hg.trial{ft_idx}(ichan, syl_inds) ) / resp.base{ichan}(itrial);
            resp.prod{ichan}(itrial,isyl) = mean( D_hg.trial{ft_idx}(ichan, prod_syl_inds) ) ; % times when syl 1 2 3 were spoken

        end

        % phoneme info
        phon_row_cons = trials_phon_rowmatch & strcmp(trials_phon.type,'consonant') & trials_phon.syl_id==isyl;
            trials.cons{itrial,isyl} = trials_phon.stim{phon_row_cons};
        phon_row_vow = trials_phon_rowmatch & strcmp(trials_phon.type,'vowel') & trials_phon.syl_id==isyl;
            trials.vow{itrial,isyl} = trials_phon.stim{phon_row_vow};

     
    end
    
    % preparatory responses........  make sure to tabulate syllable timing first
    %%%% prep period inds: after stim ends and before first syllable prod onset
    prep_inds = D_hg.time{ft_idx} > trials.t_stim_off(itrial) & D_hg.time{ft_idx} < trials.prod_syl_on(itrial,1); 
    for ichan = 1:nchans
        resp.prep{ichan}(itrial) = mean( D_hg.trial{ft_idx}(ichan, prep_inds) ) / resp.base{ichan}(itrial);
    end
    
    % transition responses
    for itrans = 1:2
        % 'transition' periods start/end halfway through the syllable
        trials.trans_on(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,1:2) + trials.prod_syl_off(itrial,1:2)]; % avg start/end
        trials.trans_off(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,2:3) + trials.prod_syl_off(itrial,2:3)]; % avg start/end
        trans_inds = D_hg.time{ft_idx} > trials.trans_on(itrial,itrans) & D_hg.time{ft_idx} < trials.trans_off(itrial,itrans); 
       for ichan = 1:nchans
           resp.trans{ichan}(itrial,itrans) = mean( D_hg.trial{ft_idx}(ichan, trans_inds) ) / resp.base{ichan}(itrial);
       end
    end
end

%% extract additional trial-specific stim information 
trials(:,{'t_stim_off','t_stim_on','has_speech_timing'}) = [];
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
    resp.n_usable_trials(ichan) = nnz(good_trials); 
    if resp.n_usable_trials(ichan) < min_trials_for_good_channel
        resp.usable_chan(ichan) = false; 
        continue; % skip stats analysis if channel had too few good trials
    else
        resp.usable_chan(ichan) = true; 
    end 
    
    % preparatory activity............................ need to not use absolute value so that we can have negative values
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials), resp.base{ichan}(good_trials)); 
    
    % rank-order selectivity
    resp.p_rank(ichan) = anova1(resp.prod{ichan}(good_trials,:),[],'off');
    

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

%%%%%%%%%%%%%%%%%% tuning to phonemic features in only one of the three positions
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

    % transition selectivity
    
    

    %%% electrode info - first matching electrode table channel
    %  for macro electrode locations, only copied their first listed location from the electrodes table into response table...
    %  ..... they generally move to deep structures over the course of trials, so this location will become less accurate for later trials
    resp.elc_info_row(ichan) = find(strcmp(elc_info.electrode , resp.chan{ichan}), 1);

end

%%%% next step to implement should be: instead of averaging, use GLM [MANOVA?] or machine learning to predict presence of multiple phon features
resp.p_prod_cons_mean  = geomean(resp.p_prod_cons,2);
resp.p_prod_vow_mean  = geomean(resp.p_prod_vow,2);
resp.p_prod_syl_mean = geomean(resp.p_prod_syl,2);
resp.p_prep_cons_mean  = geomean([resp.p_prep_cons],2);
resp.p_prep_vow_mean  = geomean([resp.p_prep_vow],2);
resp.p_prep_syl_mean  = geomean([resp.p_prep_syl],2);
    
resp = resp(resp.usable_chan,:); % remove channels with few/no usuable trials
elc_info_copy = renamevars(elc_info(resp.elc_info_row,:),'electrode','chan');
resp = join(resp, elc_info_copy(:,info_vars_to_copy)); % add elc_info to resp
resp = removevars(resp,{'elc_info_row','usable_chan'}); 
resp.sub = repmat(SUBJECT, height(resp), 1);

%%DT attempt at permutation test
