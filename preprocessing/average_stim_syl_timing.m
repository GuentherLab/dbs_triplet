%%% get the average stimulus inter-syllable time by averaging from across 3000-series subjects
%%% .... these values will be used for estimating stim timing for subjects that do not have _stimulus_syllable.txt tables available
%   for each unique triplet, this script will determine the average syl onset-to-onset timing
%   compute latency between [syl 1 and syl 2] and between [syl 1 and syl 3] onsets
%
%  also save a record of the stim syl durations; these are always identical (for a given syl) across _stimulus_syllable.txt....
%  ... so just take them from the first subject

clear
setpaths_dbs_triplet()
syl_onsets_filename = [PATH_STIM_INFO filesep 'stim_syl_onset_timing_stats.txt']; 
syl_dur_filename = [PATH_STIM_INFO filesep 'stim_syl_durations.txt']; 

subs = readtable([PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']); % TURBO path
nsubs = height(subs);
has_stimsyl = false(nsubs,1);

%% syllable onset-to-onset timing
% find the subjects that do have stimsyl table
for isub = 1:nsubs
    thissub = subs.subject{isub};
    stimsylpath = [PATH_DATA filesep thissub filesep 'Preprocessed Data' filesep 'Sync' filesep 'annot' filesep thissub '_stimulus_syllable.txt'];
    has_stimsyl(isub) = exist(stimsylpath, 'file');
end

subs_w_stimsyl = subs(has_stimsyl,:);
nsubs_w_stimsyl = height(subs_w_stimsyl);
trials_stim_trip_cat = table; 
for isub = 1:nsubs_w_stimsyl
    thissub = subs_w_stimsyl.subject{isub};
    stimsylpath = [PATH_DATA filesep thissub filesep 'Preprocessed Data' filesep 'Sync' filesep 'annot' filesep thissub '_stimulus_syllable.txt'];
    stimtrippath = [PATH_DATA filesep thissub filesep 'Preprocessed Data' filesep 'Sync' filesep 'annot' filesep thissub '_stimulus_triplet.txt'];
    trials_stim_syl = readtable(stimsylpath);
    trials_stim_trip = readtable(stimtrippath);
    ntrials = height(trials_stim_trip);
    trials_stim_trip.syl_ons2ons = nan(ntrials,2); % time between onsets of stim syl 1-2 and syl 1-3
    for itrial = 1:ntrials
       matchrow = find(trials_stim_syl.session_id == trials_stim_trip.session_id(itrial) & ...
                       trials_stim_syl.trial_id == trials_stim_trip.trial_id(itrial) & ...
                       trials_stim_syl.syl_id == 1); % find syl 1 of this trial's triplet in the syl timing table
        trials_stim_trip.syl_ons2ons(itrial,1) = trials_stim_syl.starts(matchrow+1) - trials_stim_syl.starts(matchrow); 
        trials_stim_trip.syl_ons2ons(itrial,2) = trials_stim_syl.starts(matchrow+2) - trials_stim_syl.starts(matchrow);
    end

    trials_stim_trip_cat = [trials_stim_trip_cat; trials_stim_trip]; % stack subjects into 1 table
end

stats_stim_trip = grpstats(trials_stim_trip_cat,"stim",["mean","std","min","max","range"],"DataVars","syl_ons2ons");
stats_stim_trip = sortrows(stats_stim_trip, "stim"); % alphabetize

% range within this table indicates the maximum amount that stim timing might be off if we use these estimates...
% .... (assuming the actual timing we're trying to estimate falls within this table's min and max)
mean_range = mean(stats_stim_trip.range_syl_ons2ons)
std_range = std(stats_stim_trip.range_syl_ons2ons)
max_range = max(stats_stim_trip.range_syl_ons2ons)

mkdir(PATH_STIM_INFO)
writetable(stats_stim_trip, syl_onsets_filename, "FileType","text", "WriteRowNames",true, "Delimiter","tab")

%% syllable duration
% stim syl durations are identical for a given syl, so just get durations from first subject
stimsylpath = [PATH_DATA filesep subs_w_stimsyl.subject{1} filesep 'Preprocessed Data' filesep 'Sync' filesep 'annot' filesep subs_w_stimsyl.subject{1} '_stimulus_syllable.txt'];
trials_stim_syl = readtable(stimsylpath);

% can check the table 'durstats_stim_syl' to confirm that the listed durations do not vary within a syl ID
durstats_stim_syl = grpstats(trials_stim_syl,"stim",["mean","std","min","max","range"],"DataVars","duration");
stim_syl_durations = sortrows(durstats_stim_syl(:,{'stim','mean_duration'}), 'stim');
stim_syl_durations = renamevars(stim_syl_durations, "mean_duration", "duration"); 
writetable(stim_syl_durations, syl_dur_filename, "FileType","text", "WriteRowNames",true, "Delimiter","tab")



