% look for elcs with response types
%
% updated by AM 2022/9/23

if ~exist('D_hg','var')
    load('Z:\DBS\DBS3012\Preprocessed Data\FieldTrip\DBS3012_ft_hg_trial_criteria_E.mat')
end

% Loading packages
ft_defaults
bml_defaults
format long

% Loading parameters
SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_PREPROCESSED = [PATH_SUBJECT filesep 'Preprocessed Data'];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];
PATH_FIELDTRIP = [PATH_SUBJECT '/Preprocessed Data/FieldTrip'];
ARTIFACT_CRIT = 'E'; 
SAMPLE_RATE = 100; % downsample rate in hz for high gamma traces

% analysis parameters
%%%% for baseline window, use the period from -base_win_sec(1) to -base_win_sec(2) before syl #1 stim onset
%%% baseline should end at least a few 100ms before stim onset in order to not include anticipatory activity in baseline
%%% nb: in Triplet task, time between voice offset and the subsequent trial's stim onset is likely ~2.2sec, so baseline must be shorter than this duration
base_win_sec = [1, 0.3]; 
post_speech_win_sec = 0.5; % time to include after 3rd-syllable voice offset in HG timecourse


%% load data 
cd(PATH_SYNC)

% load([PATH_FIELDTRIP filesep SUBJECT '_ft_hg_trial_criteria_' ARTIFACT_CRIT '.mat']);

trials_stim_trip = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_triplet.txt';]); % stim timing info
trials_stim_syl = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_syllable.txt';]); % stim timing info
trials_trip = readtable([PATH_ANNOT filesep SUBJECT '_produced_triplet.txt';]); % speech timing info
trials_syl = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable.txt';]); % speech timing info
trials_phon = readtable([PATH_ANNOT filesep SUBJECT '_produced_phoneme.txt';]); % speech timing info

%% get responses in predefined epochs
% 'base' = average durng pre-stim baseline
% all response values except 'base' are baseline-normalized by dividing by that trial's baseline average
ntrials = height(trials_trip);
nchans = length(D_hg.label);
nans_ch = nan(nchans,1); 
nans_tr = nan(ntrials,1); 
nans_tr2 = nan(ntrials,2); 
nans_tr3 = nan(ntrials,3); 
cel_tr = cell(ntrials,1); 

% info about our trial timing analysis window
trials = trials_stim_trip; 
trials = [trials, table(cel_tr, nans_tr3,    nans_tr3,            nans_tr3,    nans_tr3,         nans_tr2, nans_tr2,...
     'VariableNames', {'times', 'stim_syl_on', 'stim_syl_off', 'prod_syl_on', 'prod_syl_off',  'trans_on', 'trans_off' })];
trials = renamevars(trials,{'starts','ends','duration'}, {'t_stim_on','t_stim_off','dur_stim'}); % stim timing does not mark our trial boundaries
trials.starts = trials.t_stim_on - base_win_sec(1); % trial starts at beginning of baseline window
trials.ends = trials_trip.ends + post_speech_win_sec; % trial ends at fixed time after voice offset
trials.duration = trials.ends - trials.starts; 

% table containing responses during epochs for each chan
cel = repmat({nans_tr},nchans,1); % 1 value per trial per chan
cel2 = repmat({nans_tr2},nchans,1); % 2 values per trial per chan
cel3 = repmat({nans_tr3},nchans,1); % 3 values per trial per chan
resp = table(   D_hg.label, cel,   repmat({cel_tr},nchans,1),  cel,    cel,    cel3,  cel2,    nans_ch, nans_ch, nans_ch, nans_ch,  ....
  'VariableNames', {'chan', 'base', 'timecourse',             'stim', 'prep', 'syl', 'trans', 'p_prep', 'p_rank', 'p_syl', 'p_trans'}); 

% extract epoch-related responses
%%%% trials.times{itrial} use global time coordinates
%%%% ....... start at a fixed baseline window before stim onset
%%%% ....... end at a fixed time buffer after speech offset
for itrial = 1:ntrials % itrial is absolute index across sessions; does not equal "trial_id" from loaded tables
    isess = trials.session_id(itrial); 
    trial_id_in_sess = trials.trial_id(itrial); % session-relative trial number
    % get indices within the trial-specific set of timepoints of D_hg.time{itrial} that match our specified trial window
    match_time_inds = D_hg.time{itrial} > trials.starts(itrial) & D_hg.time{itrial} < trials.ends(itrial); 
    trials.times{itrial} = D_hg.time{itrial}(match_time_inds); % times in this redefined trial window... still using global time coordinates
    % get trial-relative baseline time indices; window time-locked to first stim onset
    base_inds = D_hg.time{itrial} > trials.starts(itrial) & D_hg.time{itrial} < [trials.t_stim_on(itrial) - base_win_sec(2)]; 

    % baseline activity and timecourse
    for ichan = 1:nchans
        % use mean rather than nanmean, so that trials which had artifacts marked with NaNs will be excluded
        resp.base{ichan}(itrial) = mean( D_hg.trial{itrial}(ichan, base_inds), 'includenan' ); % mean HG during baseline
        % get baseline-normalized trial timecourse
       resp.timecourse{ichan}{itrial} =  D_hg.trial{itrial}(ichan, match_time_inds) - resp.base{ichan}(itrial); 
    end
    
    % syllable-specific responses and timing
    for isyl = 1:3
        trials_syl_row = find(trials_syl.trial_id == trial_id_in_sess  &  trials_syl.syl_id == isyl & trials_syl.session_id == isess); 
        trials.prod_syl_on(itrial,isyl) = trials_syl.starts(trials_syl_row);  % syllable start time
        trials.prod_syl_off(itrial,isyl) = trials_syl.ends(trials_syl_row); % syllable end time
        syl_inds = D_hg.time{itrial} > trials.prod_syl_on(itrial,isyl) & D_hg.time{itrial} < trials.prod_syl_off(itrial,isyl); 
        for ichan = 1:nchans
            resp.syl{ichan}(itrial,isyl) = mean( D_hg.trial{itrial}(ichan, syl_inds) ) / resp.base{ichan}(itrial);
        end

        % stim timing
        trials.stim_syl_on(itrial,isyl) = trials_stim_syl.starts(trials_syl_row);  % syllable start time
        trials.stim_syl_off(itrial,isyl) = trials_stim_syl.ends(trials_syl_row); % syllable end time        
    end
    
    % preparatory responses........  make sure to tabulate syllable timing first
    %%%% prep period inds: after stim ends and before first syllable prod onset
    prep_inds = D_hg.time{itrial} > trials.t_stim_off(itrial) & D_hg.time{itrial} < trials.prod_syl_on(itrial,1); 
    for ichan = 1:nchans
        resp.prep{ichan}(itrial) = mean( D_hg.trial{itrial}(ichan, prep_inds) ) / resp.base{ichan}(itrial);
    end
    
    % transition responses
    for itrans = 1:2
        % 'transition' periods start/end halfway through the syllable
        trials.trans_on(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,1:2) + trials.prod_syl_off(itrial,1:2)]; % avg start/end
        trials.trans_off(itrial,:) = 0.5 * [trials.prod_syl_on(itrial,2:3) + trials.prod_syl_off(itrial,2:3)]; % avg start/end
        trans_inds = D_hg.time{itrial} > trials.trans_on(itrial,itrans) & D_hg.time{itrial} < trials.trans_off(itrial,itrans); 
       for ichan = 1:nchans
           resp.trans{ichan}(itrial,itrans) = mean( D_hg.trial{itrial}(ichan, trans_inds) ) / resp.base{ichan}(itrial);
       end
    end
end
trials(:,{'t_stim_off','t_stim_on'}) = [];

%% test for response types 
for ichan = 1:nchans
    good_trials = ~isnan(resp.base{ichan}); % non-artifactual trials for this channel
    if nnz(good_trials) == 0; continue; end % skip stats analysis if channel had no good trials
    
    % preparatory activity............................ need to not use absolute value so that we can have negative values
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials)); 
    
    % rank-order selectivity
    resp.p_rank(ichan) = anova1(resp.syl{ichan}(good_trials,:),[],'off');
    
    % syllable selectivity
    
    % transtiion selectivity
    
    
    
    
end
    
