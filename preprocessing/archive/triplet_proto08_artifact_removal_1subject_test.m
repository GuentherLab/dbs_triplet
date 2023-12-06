% artifact removal for a single subject - protocol 09
% Andrew Meier

clear

load(['Z:\DBS\DBS3012' filesep 'Preprocessed Data' filesep 'FieldTrip' filesep 'DBS3012_ft_raw_session.mat']);

cfg=[];
cfg.channel = 'ecog_*'; % add in dbs leads once ecog preprocessing is working
D_ecog = ft_selectdata(cfg,D);

cfg = [];
cfg.foi = 100; % high gamma
cfg.dt = 1/20; % 20 hz sampling rate
cfg.toilim = [D.time{1}(1), D.time{5}(end)]; % [first timepoint of first session, last timepoint of last session]
timefreq = bml_freqanalysis_power_wavelet(cfg, D_ecog);

% compare time points
D_ecog.time{1}(1)
D_ecog.time{5}(1)
timefreq.time(1)



