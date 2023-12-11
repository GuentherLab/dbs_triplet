% look for elcs with particualr response profiles
%
% updated by AM 2022/9/23



% Loading packages
ft_defaults
bml_defaults
format long

% % % % % Loading parameters
% SUBJECT='DBS3012';
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

% add the following variables to the electrodes response table... use 'electrode'/'chan' as key variable
info_vars_to_copy = {'chan','type','connector','port','strip','comment','target','side','nat_x','nat_y','nat_z',...
    'leadDBS_x','leadDBS_y','leadDBS_z','tkRAS_x','tkRAS_y','tkRAS_z','mni_linear_x','mni_linear_y','mni_linear_z',...
    'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z','fs_anatomy','fs_mni152_cvs_x','fs_mni152_cvs_y','fs_mni152_cvs_z',...
	'DISTAL_label_1','DISTAL_weight_1','DISTAL_label_2','DISTAL_weight_2','DISTAL_label_3','DISTAL_weight_3',...
    'MOREL_label_1','MOREL_weight_1','MOREL_label_2','MOREL_weight_2','MOREL_label_3','MOREL_weight_3',...
    'HCPMMP1_label_1','HCPMMP1_weight_1','HCPMMP1_label_2','HCPMMP1_weight_2'};

%% load data 
cd(PATH_SYNC)

% if ~exist('D_hg','var')
%     load([PATH_FIELDTRIP, filesep, SUBJECT, '_ft_hg_trial_ref_criteria_E_denoised.mat'])
% end

load([PATH_FIELDTRIP filesep SUBJECT '_ft_hg_trial_ref_criteria_' ARTIFACT_CRIT '_denoised.mat']);

trials_stim_trip = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_triplet.txt';]); % stim timing info
trials_stim_syl = readtable([PATH_ANNOT filesep SUBJECT '_stimulus_syllable.txt';]); % stim timing info
trials_trip = readtable([PATH_ANNOT filesep SUBJECT '_produced_triplet.txt';]); % speech timing info
trials_syl = readtable([PATH_ANNOT filesep SUBJECT '_produced_syllable.txt';]); % speech timing info
trials_phon = readtable([PATH_ANNOT filesep SUBJECT '_produced_phoneme.txt';]); % speech timing info

elc_info = readtable([PATH_ANNOT filesep SUBJECT '_electrode.txt';]); 

% some subjects don't have _stimulus_syllable.txt - for these, maybe try to derive it from _stimulus_triplet.txt and expected durations? 


% account for some tables using onset/duration convention vs. starts/ends convention
if ~any(contains(trials_trip.Properties.VariableNames,'starts'))
    trials_trip.starts = trials_trip.onset;
end
if ~any(contains(trials_trip.Properties.VariableNames,'ends'))
    trials_trip.ends = trials_trip.starts + trials_trip.duration; 
end

%% get responses in predefined epochs
% 'base' = average during pre-stim baseline
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

unqsyl = unique([trials.stim1(:); trials.stim2(:); trials.stim3(:)]);
n_unqsyl = length(unqsyl);
trials.prod_syl_present = false(ntrials, n_unqsyl);

% table containing responses during epochs for each chan
cel = repmat({nans_tr},nchans,1); % 1 value per trial per chan
cel2 = repmat({nans_tr2},nchans,1); % 2 values per trial per chan
cel3 = repmat({nans_tr3},nchans,1); % 3 values per trial per chan
nansyl = nan(nchans, n_unqsyl); 
resp = table(   D_hg.label, cel,   repmat({cel_tr},nchans,1),  cel3,    cel,    cel3,  cel2,    nans_ch, nans_ch,  nansyl,        nansyl, ....
  'VariableNames', {'chan', 'base', 'timecourse',             'stim', 'prep', 'prod', 'trans', 'p_prep', 'p_rank','p_prep_syl', 'p_prod_syl'}); 

%% PROBABLY NEED TO REWRITE THE TIME-MATCHING LINES TO LOOK FOR MATCHING TIMES ACROSS ALL TRIALS.....
%%% .... due to missing trials in the fieldtrip structures, trialtable indices may not match to fieldtrip trial indices

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
    trials_phon_rowmatch = trials_phon.trial_id == trials.trial_id(itrial) & ...
                           trials_phon.session_id == trials.session_id(itrial); 
    for isyl = 1:3
            % prod timing
        trials_syl_row = find(trials_syl.trial_id == trial_id_in_sess  &  trials_syl.syl_id == isyl & trials_syl.session_id == isess); 
        trials.prod_syl_on(itrial,isyl) = trials_syl.starts(trials_syl_row);  % syllable start time
        trials.prod_syl_off(itrial,isyl) = trials_syl.ends(trials_syl_row); % syllable end time
                % stim timing
        trials.stim_syl_on(itrial,isyl) = trials_stim_syl.starts(trials_syl_row);  % syllable start time
        trials.stim_syl_off(itrial,isyl) = trials_stim_syl.ends(trials_syl_row); % syllable end time   
            % time indices of epochs
        stim_syl_inds = D_hg.time{itrial} > trials.stim_syl_on(itrial,isyl) & D_hg.time{itrial} < trials.stim_syl_off(itrial,isyl); 
        prod_syl_inds = D_hg.time{itrial} > trials.prod_syl_on(itrial,isyl) & D_hg.time{itrial} < trials.prod_syl_off(itrial,isyl); 
        for ichan = 1:nchans
            resp.stim{ichan}(itrial,isyl) = mean( D_hg.trial{itrial}(ichan, stim_syl_inds) ) ; % times when syl 1 2 3 were played
            %             resp.prod{ichan}(itrial,isyl) = mean( D_hg.trial{itrial}(ichan, syl_inds) ) / resp.base{ichan}(itrial);
            resp.prod{ichan}(itrial,isyl) = mean( D_hg.trial{itrial}(ichan, prod_syl_inds) ) ; % times when syl 1 2 3 were spoken

        end

        % phoneme info
        phon_row_cons = trials_phon_rowmatch & strcmp(trials_phon.type,'consonant') & trials_phon.syl_id==isyl;
        trials.cons{itrial,isyl} = trials_phon.stim{phon_row_cons};
        phon_row_vow = trials_phon_rowmatch & strcmp(trials_phon.type,'vowel') & trials_phon.syl_id==isyl;
        trials.vow{itrial,isyl} = trials_phon.stim{phon_row_vow};

     
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

for itrial = 1:ntrials
    for isyl = 1:n_unqsyl
        this_syl = unqsyl{isyl};
        if any(strcmp(this_syl, {trials.stim1{itrial}, trials.stim2{itrial}, trials.stim3{itrial}}))
            trials.prod_syl_present(itrial,isyl) = true;
        end
    end
end


%% test for response types 
resp.elc_info_row = nan(nchans,1); 
resp.nodata = false(nchans,1); 
for ichan = 1:nchans
    good_trials = ~isnan(resp.base{ichan}); % non-artifactual trials for this channel
    if nnz(good_trials) == 0
       resp.nodata(ichan) = 1; % if no data in this channel, remove it later
        continue; 
    end % skip stats analysis if channel had no good trials
    
    % preparatory activity............................ need to not use absolute value so that we can have negative values
    [~, resp.p_prep(ichan)] = ttest(resp.prep{ichan}(good_trials), ones(size(resp.prep{ichan}(good_trials)))); 
    
    % rank-order selectivity
    resp.p_rank(ichan) = anova1(resp.prod{ichan}(good_trials,:),[],'off');
    
    % syllable prep selectivity
    for isyl = 1:n_unqsyl
        resp.p_prep_syl(ichan,isyl) = anova1(resp.prep{ichan}(good_trials), trials.prod_syl_present(good_trials,isyl), "off"); 
    end

    % syllable prod selectivity - full production period
    for isyl = 1:n_unqsyl
        resp.p_prod_syl(ichan,isyl) = anova1(mean(resp.prod{ichan}(good_trials), 2), trials.prod_syl_present(good_trials,isyl), "off"); 
    end

    % prep-syl with position
    [p_prep_syl1, p_prep_anovatab] = anova1(resp.prep{ichan}(good_trials),trials.stim1(good_trials),'off');
    resp.p_prep_syl1(ichan) = p_prep_syl1;

    [p_prep_syl2, p_prep_anovatab2] = anova1(resp.prep{ichan}(good_trials),trials.stim2(good_trials),'off');
    resp.p_prep_syl2(ichan) = p_prep_syl2;

    [p_prep_syl3, p_prep_anovatab3] = anova1(resp.prep{ichan}(good_trials),trials.stim3(good_trials),'off');
    resp.p_prep_syl3(ichan) = p_prep_syl3;

     % prod-syl using the response from only when the syllable in question is being produced
     for ipos = 1:3
         resp.p_prod_syl_position(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), trials{good_trials,['stim',num2str(ipos)]},'off');
     end

       % stim-syl using the response from only when the syllable in question is being played
     for ipos = 1:3
         resp.p_stim_syl_position(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), trials{good_trials,['stim',num2str(ipos)]},'off');
     end

     % prod-consonant using the response from only when the syllable in question is being produced
     for ipos = 1:3
         resp.p_prod_cons_position(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), trials.cons(good_trials, ipos),'off');
     end

      % prod-vowel using the response from only when the syllable in question is being produced
     for ipos = 1:3
         resp.p_prod_vow_position(ichan,ipos) = anova1(resp.prod{ichan}(good_trials,ipos), trials.vow(good_trials, ipos),'off');
     end

          % stim-consonant using the response from only when the syllable in question is being produced
     for ipos = 1:3
         resp.p_stim_cons_position(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), trials.cons(good_trials, ipos),'off');
     end

      % stim-vowel using the response from only when the syllable in question is being produced
     for ipos = 1:3
         resp.p_stim_vow_position(ichan,ipos) = anova1(resp.stim{ichan}(good_trials,ipos), trials.vow(good_trials, ipos),'off');
     end

    % prep-syl - syllable-selective activity during prep period
    for isyl = 1:n_unqsyl
        this_syl_trials = good_trials & trials.prod_syl_present(:,isyl);
        p_prep_syl(ichan, isyl) = ttest2(resp.prep{ichan}(this_syl_trials), resp.prep{ichan}(~this_syl_trials));
    end

    % transition selectivity
    
    

    %%% electrode info - first matching electrode table channel
    %  for macro electrode locations, only copied their first listed location from the electrodes table into response table...
    %  ..... they generally move to deep structures over the course of trials, so this location will become less accurate for later trials
    resp.elc_info_row(ichan) = find(strcmp(elc_info.electrode , resp.chan{ichan}), 1);

end
    
resp = resp(~resp.nodata,:); % remove empty channels
elc_info_copy = renamevars(elc_info(resp.elc_info_row,:),'electrode','chan');
resp = join(resp, elc_info_copy(:,info_vars_to_copy)); % add elc_info to resp
resp = removevars(resp,{'elc_info_row','nodata'}); 
resp.sub = repmat(SUBJECT, height(resp), 1);