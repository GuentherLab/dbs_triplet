% artifact removal for a single subject - protocol 09
% Andrew Meier

clear

%% Defining artifact rejection parameters
CRITERIA = 'AM_1'; % high thresh = 3SDs, low thresh = 2SDs

HIGH_PASS_FILTER = 'yes'; %should a high pass filter be applied
HIGH_PASS_FILTER_FREQ = 1; %cutoff frequency of high pass filter

do_bsfilter = 'yes'; 
line_noise_harm_freqs=[60 120 180 240]; % for notch filters for 60hz harmonics

ECOG_ENVELOPE_BIN_SIZE_SECONDS = 1; %envelope bin size in seconds
ECOG_THRESHOLD_STD_FACTORS = [2, 3]; %factors to determine detection thresholds 
ECOG_CONSOLIDATION_TIME_TOLERANCE = 3; %min time allowed between adjacent artifacts
ELECTRODE_COVERAGE_THRESHOLD = 0.5; %max allowed fraction of time with artifacts
CONNECTOR_THRESHOLD = [4, 8]; %detection threshold for faulty electrodes in connector

%loading packages 
ft_defaults
bml_defaults
format long

%defining paths
SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
cd(PATH_SYNC)

session= bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
% if doing manual rejection 
% artifact_manual = bml_annot_read(['annot/' SUBJECT '_artifact_manual.txt']);

% artifact removal for a single subject - protocol 09
% Andrew Meier

clear

%% Defining artifact rejection parameters
CRITERIA = 'AM_1'; % high thresh = 3SDs, low thresh = 2SDs

HIGH_PASS_FILTER = 'yes'; %should a high pass filter be applied
HIGH_PASS_FILTER_FREQ = 1; %cutoff frequency of high pass filter

do_bsfilter = 'yes'; 
line_noise_harm_freqs=[60 120 180 240]; % for notch filters for 60hz harmonics

ECOG_ENVELOPE_BIN_SIZE_SECONDS = 1; %envelope bin size in seconds
ECOG_THRESHOLD_STD_FACTORS = [2, 3]; %factors to determine detection thresholds 
ECOG_CONSOLIDATION_TIME_TOLERANCE = 3; %min time allowed between adjacent artifacts
ELECTRODE_COVERAGE_THRESHOLD = 0.5; %max allowed fraction of time with artifacts
CONNECTOR_THRESHOLD = [4, 8]; %detection threshold for faulty electrodes in connector

%loading packages 
ft_defaults
bml_defaults
format long

%defining paths
SUBJECT='DBS3012';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
cd(PATH_SYNC)

session= bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
% if doing manual rejection 
% artifact_manual = bml_annot_read(['annot/' SUBJECT '_artifact_manual.txt']);

load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_raw_session.mat']);
nTrials = numel(D.trial);

%mask the NaN values you have with the value of 0 
cfg=[];
cfg.value = 0;
cfg.remask_nan = true;
D_mask = bml_mask(cfg,D);


%high pass filtering


cfg=[];
cfg.hpfilter=HIGH_PASS_FILTER;	 % setting a high pass filter
cfg.hpfreq=HIGH_PASS_FILTER_FREQ; % 1Hz cutoff frequency	
cfg.hpfilttype='but';    		 % using butterworth filter
cfg.hpfiltord=5;         		 % of order 5
cfg.hpfiltdir='twopass'; 		 % two passes not to introduce time shift
cfg.bsfilter=do_bsfilter;
cfg.bsfreq= [line_noise_harm_freqs-1; line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'}; %include channels to filter
D_hpf = ft_preprocessing(cfg,D_mask);

cfg=[];
cfg.channel = 'ecog_*'; % add in dbs leads once ecog preprocessing is working
D_hpf_ecog = ft_selectdata(cfg,D_hpf);


%% get trial timepoints
% adapted from protocol P10
coding = bml_annot_read(['annot/' SUBJECT '_coding.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);


epoch = table(); %creating empty table
epoch.syl1_onset = coding.syl1_onset; %copying onset of syllable 1 from coding table
epoch.trial_id = coding.trial_id; 
epoch.session_id = coding.session_id; 
epoch.starts = coding.syl1_onset - 2; %creating starts
epoch.ends = coding.syl1_onset + 3; %creating ends
epoch = bml_annot_table(epoch); %creating id and duration columns
epoch = epoch(~ismissing(epoch.starts),:); %removing rows with NaNs

%% get high gamma for artifact rejection


hg_lowerfreq = 62;
hg_upperfreq = 150; 
nfreqs = 1; 
TF_FOI = round(logspace(log10(hg_lowerfreq),log10(hg_upperfreq),4),nfreqs,'significant'); % frequencies of interest to extract
TF_RATE = 20; %Hz Sampling rate for wavelet transform
cfg = []; 
cfg.foi = TF_FOI; % only get the freq bands within high gamma
cfg.dt = 1/TF_RATE;
% cfg.toilim = [epoch.starts, epoch.ends]; %%% analyze data within trial boundaries
cfg.toilim = [epoch.starts(1), epoch.ends(end)]; %%% analyze data within trial boundaries
D_hg_ecog = bml_freqanalysis_power_wavelet(cfg, D_hpf_ecog);
% D_hg_ecog.trial = squeeze(mean(D_hg_ecog.powspctrm,3)); % average across HG freq bands



nTrials = numel(D.trial);

%mask the NaN values you have with the value of 0 
cfg=[];
cfg.value = 0;
cfg.remask_nan = true;
D_mask = bml_mask(cfg,D);


%high pass filtering


cfg=[];
cfg.hpfilter=HIGH_PASS_FILTER;	 % setting a high pass filter
cfg.hpfreq=HIGH_PASS_FILTER_FREQ; % 1Hz cutoff frequency	
cfg.hpfilttype='but';    		 % using butterworth filter
cfg.hpfiltord=5;         		 % of order 5
cfg.hpfiltdir='twopass'; 		 % two passes not to introduce time shift
cfg.bsfilter=do_bsfilter;
cfg.bsfreq= [line_noise_harm_freqs-1; line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'}; %include channels to filter
D_hpf = ft_preprocessing(cfg,D_mask);

cfg=[];
cfg.channel = 'ecog_*'; % add in dbs leads once ecog preprocessing is working
D_hpf_ecog = ft_selectdata(cfg,D_hpf);


%% get trial timepoints
% adapted from protocol P10
coding = bml_annot_read(['annot/' SUBJECT '_coding.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);


epoch = table(); %creating empty table
epoch.syl1_onset = coding.syl1_onset; %copying onset of syllable 1 from coding table
epoch.trial_id = coding.trial_id; 
epoch.session_id = coding.session_id; 
epoch.starts = coding.syl1_onset - 2; %creating starts
epoch.ends = coding.syl1_onset + 3; %creating ends
epoch = bml_annot_table(epoch); %creating id and duration columns
epoch = epoch(~ismissing(epoch.starts),:); %removing rows with NaNs

%% get high gamma for artifact rejection


hg_lowerfreq = 62;
hg_upperfreq = 150; 
nfreqs = 1; 
TF_FOI = round(logspace(log10(hg_lowerfreq),log10(hg_upperfreq),4),nfreqs,'significant'); % frequencies of interest to extract
TF_RATE = 20; %Hz Sampling rate for wavelet transform
cfg = []; 
cfg.foi = TF_FOI; % only get the freq bands within high gamma
cfg.dt = 1/TF_RATE;
% cfg.toilim = [epoch.starts, epoch.ends]; %%% analyze data within trial boundaries
cfg.toilim = [epoch.starts(1), epoch.ends(end)]; %%% analyze data within trial boundaries
D_hg_ecog = bml_freqanalysis_power_wavelet(cfg, D_hpf_ecog);
% D_hg_ecog.trial = squeeze(mean(D_hg_ecog.powspctrm,3)); % average across HG freq bands


