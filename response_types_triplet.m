% look for elcs with particualr response profiles

setpaths_dbs_triplet()

%% parameters 
% % % % % data-loading parameters
vardefault('SUBJECT','DBS3012');

DATE=datestr(now,'yyyymmdd');

setpaths_dbs_triplet()
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_PREPROCESSED = [PATH_SUBJECT filesep 'Preprocessed Data'];
PATH_FIELDTRIP = [PATH_PREPROCESSED filesep 'FieldTrip']; 
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];
    stimsylpath = [PATH_ANNOT filesep SUBJECT '_stimulus_syllable.txt']; 
ARTIFACT_CRIT = 'E'; 
SAMPLE_RATE = 100; % downsample rate in hz for high gamma traces

% files with info to be used if stim syl timing info is missing for a subject
%%% created by average_stim_syl_timing.m
syl_onsets_filename = [PATH_STIM_INFO filesep 'stim_syl_onset_timing_stats']; 
syl_dur_filename = [PATH_STIM_INFO filesep 'stim_syl_durations']; 

% analysis parameters
%%%% for baseline window, use the period from -base_win_sec(1) to -base_win_sec(2) before syl #1 stim onset
%%% baseline should end at least a few 100ms before stim onset in order to not include anticipatory activity in baseline
%%% nb: in Triplet task, time between voice offset and the subsequent trial's stim onset is likely ~2.2sec, so baseline must be shorter than this duration
base_win_sec = [1, 0.3]; 
post_speech_win_sec = 0.5; % time to include after 3rd-syllable voice offset in HG timecourse
min_trials_for_good_channel = 4; 

% for responses during syl 1 production, start the analyzed 'speech period' this early in seconds to capture pre-sound muscle activation
% also end the prep period this early
%%% extending the window in this way for syl 2 and syl 3 might be trickier, because there is often very little time between the preceding offset and the syl2/3 onset
prod_syl1_window_extend_start = 0;  

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


%% load and organize data
load([PATH_FIELDTRIP filesep SUBJECT '_ft_hg_trial_ref_criteria_' ARTIFACT_CRIT '_denoised.mat']);

% for some subjects (3030,4061,4072,4077,4078,4080,4084,4085), AM had to make an xlsx version of this table
%%% for an unknown reason, distal and morel labels were being imported as nans....
%%% ..... this was due to some unexpected cell data formatting (usually in the first 9 rows) of the table....
%%% ..... fixed this problem by filling in one cell in these rows of distal/morel with a space...
%%% ...... which made them import correctly as cells/strings
if exist([PATH_ANNOT filesep SUBJECT '_electrode_fixed-distal-morel-labels.xlsx'], 'file')
    elc_info = readtable([PATH_ANNOT filesep SUBJECT '_electrode_fixed-distal-morel-labels.xlsx']); 
else
    elc_info = readtable([PATH_ANNOT filesep SUBJECT '_electrode.txt']); 
end

trials_prod_syl = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable.txt';]); % load syllable timing info
trials_stim_trip = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_triplet.txt';]); % stim timing info
trials_prod_trip = readtable([PATH_ANNOT filesep SUBJECT '_produced_triplet.txt';]); % speech timing info
trials_phon = readtable([PATH_ANNOT filesep SUBJECT '_produced_phoneme.txt';]); % speech timing info

% account for some tables using onset/duration convention vs. starts/ends convention
if ~any(contains(trials_prod_trip.Properties.VariableNames,'starts'))
    trials_prod_trip.starts = trials_prod_trip.onset;
end
if ~any(contains(trials_prod_trip.Properties.VariableNames,'ends'))
    trials_prod_trip.ends = trials_prod_trip.starts + trials_prod_trip.duration; 
end

% remove trials from stim trials table that do not have a corresponding fieldtrip event
%%% for the purposes of finding a fully overlapping ft trial, set trialtable's trial start as beginning of baseline and trial end as end of production plus buffer
cfg = [];
cfg.trials = trials_stim_trip; 
    cfg.trials.starts = cfg.trials.starts - base_win_sec(1); % adjust trial start to be beginning of baseline period
    for itrial = 1:height(trials_stim_trip)
        matchrow = trials_prod_trip.session_id==trials_stim_trip.session_id(itrial) & trials_prod_trip.trial_id==trials_stim_trip.trial_id(itrial);
        if any(matchrow) % if this trial has stim timing and prod timing, end trial at prod end plus buffer
            cfg.trials.ends(itrial) = trials_prod_trip.ends(matchrow) + post_speech_win_sec; 
        elseif ~any(matchrow) % if this trial has stim timing but no prod timing, add buffer to end of stim timing
            cfg.trials.ends(itrial) = cfg.trials.ends(itrial) + 1 + post_speech_win_sec; 
        end
    end
cfg.plot_times = 0;
[trials, trials_ft]  = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D_hg);

% rename vars
trials.syl = [trials.stim1, trials.stim2, trials.stim3]; 
trials = removevars(trials,{'stim1','stim2','stim3'});

% vars for table construction
%%% note that these variables are taken from stim trials table, which may be larger than 'trials' variable....
%%% ... because some trials listed in trials_stim_trip may have been cut from 'trials' due to not matching up with fieldtrip data
ntrials_stim = height(trials_stim_trip); 
nans_tr_stim = nan(ntrials_stim,1); 
cel_tr_stim = cell(ntrials_stim,1); 

%%% check whether subject is missing _stimulus_syllable.txt
%%% if it is missing, use stim timing estimates averaged from other subjects
if exist(stimsylpath, 'file') % subject has stim syl timing
    trials_stim_syl = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_syllable.txt';]); % stim timing info
elseif ~exist(stimsylpath, 'file') % subject doesn't have stim syl timing
    stim_syl_durations = readtable(syl_dur_filename, ReadRowNames=true); % fixed durations from other subjects
    stats_stim_trip = readtable(syl_onsets_filename, ReadRowNames=true); % average syl onset times from other subjects

    % construct an estimated stim syl timing table
    
    trials_stim_syl = table(nans_tr_stim, nans_tr_stim, nans_tr_stim,   nans_tr_stim,   nans_tr_stim,    nans_tr_stim, cel_tr_stim, 'VariableNames',...
                           {'starts',     'ends',       'duration',     'session_id',   'trial_id',      'syl_id',     'stim'}); 
        trials_stim_syl = repmat(trials_stim_syl,3,1); % triple in height because we have 3 syls per trial
        trials_stim_syl.session_id = repelem(trials_stim_trip.session_id, 3, 1); 
        trials_stim_syl.trial_id = repelem(trials_stim_trip.trial_id, 3, 1); 
        trials_stim_syl.syl_id = repmat([1:3]', ntrials_stim, 1); 
        sylcat = [trials_stim_trip.stim1, trials_stim_trip.stim2, trials_stim_trip.stim3;]';
            trials_stim_syl.stim = sylcat(:);

    for itrial_stim_trip = 1:ntrials_stim % this counter refers to full trials (triplets), not individual stim
        triptabstats_row = string(trials_stim_trip.stim{itrial_stim_trip}) == string(stats_stim_trip.stim); 
        isess = trials_stim_trip.session_id(itrial_stim_trip); 
        trial_id_in_sess = trials_stim_trip.trial_id(itrial_stim_trip); % session-relative trial number
        for isyl = 1:3
            trials_stim_syl_row = find(trials_stim_syl.session_id == isess & trials_stim_syl.trial_id == trial_id_in_sess  &  trials_stim_syl.syl_id == isyl); 
            syldurtab_row = string(trials_stim_syl.stim{itrial_stim_trip}) == string(stim_syl_durations.stim); 
            this_syl_dur = stim_syl_durations.duration(syldurtab_row);
            trials_stim_syl.duration(trials_stim_syl_row) = this_syl_dur;
            switch isyl
                case 1
                    trials_stim_syl.starts(trials_stim_syl_row) = trials_stim_trip.starts(itrial_stim_trip); % syl 1 starts at beginning of triplet
                    trials_stim_syl.ends(trials_stim_syl_row)   = trials_stim_trip.starts(itrial_stim_trip) + this_syl_dur; 
                case 2              % find timepoints relative to triplet beginning
                    ons1_to_ons2 = stats_stim_trip.mean_syl_ons2ons_1(triptabstats_row); % time between syl 1 onset and syl 2 onset (estimated)
                    trials_stim_syl.starts(trials_stim_syl_row) = trials_stim_trip.starts(itrial_stim_trip) + ons1_to_ons2; 
                    trials_stim_syl.ends(trials_stim_syl_row)   = trials_stim_trip.starts(itrial_stim_trip) + ons1_to_ons2 + this_syl_dur; 
                case 3              % find timepoints relative to triplet ending
                    trials_stim_syl.starts(trials_stim_syl_row) = trials_stim_trip.ends(itrial_stim_trip) - this_syl_dur;
                    trials_stim_syl.ends(trials_stim_syl_row)   = trials_stim_trip.ends(itrial_stim_trip); % syl 3 ends at ending of triplet
            end
        end
    end

end

%% get responses in predefined epochs
%%% 'base' = average during pre-stim baseline

% vars for table construction
nchans = length(D_hg.label);
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
    speechtrial_match = trials_prod_trip.session_id == isess & trials_prod_trip.trial_id == trial_id_in_sess;
    if any(speechtrial_match)
        trials.has_speech_timing(itrial) = true; 
        trials.ends(itrial) = trials_prod_trip.ends(speechtrial_match) + post_speech_win_sec;
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
        trials.stim_syl_off(itrial,isyl) = trials_stim_syl.ends(trials_stim_syl_row); % stim syllable end time  
        
           % prod timing
        trials_prod_syl_row = find(trials_prod_syl.trial_id == trial_id_in_sess  &  trials_prod_syl.syl_id == isyl & trials_prod_syl.session_id == isess); 
        trials.prod_syl_on(itrial,isyl) = trials_prod_syl.starts(trials_prod_syl_row);  % prod syllable start time
        trials.prod_syl_off(itrial,isyl) = trials_prod_syl.ends(trials_prod_syl_row); % prod syllable end time
    end

    % get indices within the trial-specific set of timepoints of D_hg.time{ft_idx} that match our specified trial window
    match_time_inds = D_hg.time{ft_idx} > trials.starts(itrial) & D_hg.time{ft_idx} < trials.ends(itrial); 
    trials.times{itrial} = D_hg.time{ft_idx}(match_time_inds); % times in this redefined trial window... still using global time coordinates

    % get trial-relative baseline time indices; window time-locked to first stim onset
    base_inds = D_hg.time{ft_idx} > [trials.stim_syl_on(itrial,1) - base_win_sec(1)] & D_hg.time{ft_idx} < [trials.stim_syl_on(itrial,1) - base_win_sec(2)]; 

    % baseline activity and timecourse
    for ichan = 1:nchans
        % use mean rather than nanmean, so that trials which had artifacts marked with NaNs will be excluded
        resp.base{ichan}(itrial) = mean( D_hg.trial{ft_idx}(ichan, base_inds), 'includenan' ); % mean HG during baseline
        % get baseline-normalized trial timecourse
       resp.timecourse{ichan}{itrial} =  D_hg.trial{ft_idx}(ichan, match_time_inds) - resp.base{ichan}(itrial); 
    end
    
    % syllable-specific responses
    trials_phon_rowmatch = trials_phon.trial_id == trials.trial_id(itrial) & ...
                           trials_phon.session_id == trials.session_id(itrial); 
    for isyl = 1:3
        stim_syl_inds = D_hg.time{ft_idx} > trials.stim_syl_on(itrial,isyl) & D_hg.time{ft_idx} < trials.stim_syl_off(itrial,isyl); % time idx when syl 1 2 3 were played
        if isyl == 1 % for syl 1, start the analyzed 'speech window' early to account for pre-speech muscle activity
            prod_syl_inds = D_hg.time{ft_idx} > [trials.prod_syl_on(itrial,isyl)-prod_syl1_window_extend_start] & D_hg.time{ft_idx} < trials.prod_syl_off(itrial,isyl); % time idx when syl 1 was spoken
        elseif ismember(isyl, [2 3])
            prod_syl_inds = D_hg.time{ft_idx} > trials.prod_syl_on(itrial,isyl) & D_hg.time{ft_idx} < trials.prod_syl_off(itrial,isyl); % time idx when syl 2 3 were spoken
        end
        for ichan = 1:nchans
            resp.stim{ichan}(itrial,isyl) = mean( D_hg.trial{ft_idx}(ichan, stim_syl_inds) ) - resp.base{ichan}(itrial); 
            resp.prod{ichan}(itrial,isyl) = mean( D_hg.trial{ft_idx}(ichan, prod_syl_inds) ) - resp.base{ichan}(itrial); 
        end

        % phoneme info
        phon_row_cons = trials_phon_rowmatch & strcmp(trials_phon.type,'consonant') & trials_phon.syl_id==isyl;
            trials.cons{itrial,isyl} = trials_phon.stim{phon_row_cons};
        phon_row_vow = trials_phon_rowmatch & strcmp(trials_phon.type,'vowel') & trials_phon.syl_id==isyl;
            trials.vow{itrial,isyl} = trials_phon.stim{phon_row_vow};

    end
    
    % preparatory responses
    %%%% prep period inds: after stim ends and before first syllable prod onset, adjusted by prod_syl1_window_extend_start
    prep_inds = D_hg.time{ft_idx} > trials.stim_syl_off(itrial,3) & D_hg.time{ft_idx} < [trials.prod_syl_on(itrial,1) - prod_syl1_window_extend_start]; 
    for ichan = 1:nchans
        resp.prep{ichan}(itrial) = mean( D_hg.trial{ft_idx}(ichan, prep_inds) , 'includenan' ) - resp.base{ichan}(itrial);
    end

    % inter-syllable transition period responses
    for itrans = 1:2
        % 'transition' periods start/end halfway through the syllable
        trials.trans_on(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,1:2) + trials.prod_syl_off(itrial,1:2)]; % avg start/end
        trials.trans_off(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,2:3) + trials.prod_syl_off(itrial,2:3)]; % avg start/end
        trans_inds = D_hg.time{ft_idx} > trials.trans_on(itrial,itrans) & D_hg.time{ft_idx} < trials.trans_off(itrial,itrans); 
       for ichan = 1:nchans
           resp.trans{ichan}(itrial,itrans) = mean( D_hg.trial{ft_idx}(ichan, trans_inds) ) - resp.base{ichan}(itrial);
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
    if resp.n_usable_trials(ichan) < min_trials_for_good_channel
        resp.usable_chan(ichan) = false; 
        continue; % skip stats analysis if channel had too few good trials
    else
        resp.usable_chan(ichan) = true; 
    end 
    
    % preparatory activity
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials), zeroes_vec) ; 
    
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
resp.sub = repmat(SUBJECT, height(resp), 1);
resp = movevars(resp,{'sub','chan','fs_anatomy','MOREL_label_1','DISTAL_label_1'},'Before',1);