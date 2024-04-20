 %%%% load triplet subject-specific tables for stim and speech timing

setpaths_dbs_triplet()
field_default('op','base_win_sec',[1, 0.3]);


% for some subjects (3030,4061,4072,4077,4078,4080,4084,4085), AM had to make an xlsx version of this table
%%% for an unknown reason, distal and morel labels were being imported as nans....
%%% ..... this was due to some unexpected cell data formatting (usually in the first 9 rows) of the table....
%%% ..... fixed this problem by filling in one cell in these rows of distal/morel with a space...
%%% ...... which made them import correctly as cells/strings
if exist([PATH_ANNOT filesep op.sub '_electrode_fixed-distal-morel-labels.xlsx'], 'file')
    elc_info = readtable([PATH_ANNOT filesep op.sub '_electrode_fixed-distal-morel-labels.xlsx']); 
else
    elc_info = readtable([PATH_ANNOT filesep op.sub '_electrode.txt']); 
end

trials_prod_syl = readtable([PATH_ANNOT filesep op.sub '_produced_syllable.txt';]); % load syllable timing info
trials_stim_trip = readtable([PATH_ANNOT filesep op.sub '_stimulus_triplet.txt';]); % stim timing info
trials_prod_trip = readtable([PATH_ANNOT filesep op.sub '_produced_triplet.txt';]); % speech timing info
trials_phon = readtable([PATH_ANNOT filesep op.sub '_produced_phoneme.txt';]); % speech timing info

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
    trials_stim_syl = readtable([PATH_ANNOT filesep op.sub '_stimulus_syllable.txt';]); % stim timing info
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