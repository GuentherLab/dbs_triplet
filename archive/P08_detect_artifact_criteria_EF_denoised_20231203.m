% Detect artifacts, create an artifact annotation table
%%% Triplet dataset version
%%% .... uses power-averaging function from dbs_seq_analysis

% function P08_detect_artifact_criteria_E(SUBJECT, param)

%% load packages
ft_defaults
bml_defaults
format long

%% Defining paths, loading parameters
SUBJECT='DBS3012';

%%% choose an artifact criterion version
ARTIFACT_CRIT = 'E'; % 70-250hz high gamma; identifier for the criteria implemented in this script
% ARTIFACT_CRIT = 'F'; % beta ; identifier for the criteria implemented in this script
% ARTIFACT_CRIT = 'G'; %identifier for the criteria implemented in this script

% highpass and line-noise filtered have already been applied to vibration-denoised data; do not apply again
high_pass_filter = 'no'; %
    high_pass_filter_freq = 1; %cutoff frequency of high pass filter
do_bsfilter = 'no'; 
    line_noise_harm_freqs=[60 120 180 240]; % for notch filters for 60hz harmonics
sample_rate = 100; % downsample rate in hz for high gamma traces

%%%%%%%%%%%%%%% PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%
PATH_DATASET = 'Z:\DBS';
PATH_SUB = [PATH_DATASET filesep SUBJECT];  
PATH_GROUPANALYSIS = [PATH_DATASET filesep 'triplet_results_am']; 
%
PATH_PREPROCESSED = [PATH_SUB filesep 'Preprocessed Data']; % derivatives
PATH_SYNC = [PATH_PREPROCESSED filesep 'Sync'];
PATH_ANALYSIS = [PATH_PREPROCESSED filesep 'analysis'];
PATH_FIELDTRIP = [PATH_PREPROCESSED filesep 'fieldtrip'];
PATH_ANNOT = [PATH_SYNC filesep 'annot'];
%
PATH_RAW = [PATH_SUB filesep 'Raw'];
PATH_AUDIO = [PATH_RAW filesep 'audio']; 
%
PATH_ART_PROTOCOL = [PATH_DATASET filesep 'Batch' filesep 'P08_artifact_criteria_', ARTIFACT_CRIT]; % group folder where artifact protocol and results are stored
PATH_ARTPROT_FIGURES = [PATH_ART_PROTOCOL filesep 'figures']; 
FILEPATH_ARTPROT_ANNOT = [PATH_ANNOT filesep SUBJECT '_artifact_criteria_' ARTIFACT_CRIT '.txt'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(PATH_SYNC)

session= readtable(['annot' filesep SUBJECT '_session.txt']);
electrodes = readtable(['annot' filesep SUBJECT '_electrode.txt']);

%% Loading FieldTrip data 
load([PATH_FIELDTRIP filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'],'D','loaded_epoch');
ntrials = numel(D.trial);

%remasking nans with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan=true;
D_trials=bml_mask(cfg,D);

%% loading electrode type band table
if ~exist('el_band','var')
  param = readtable([PATH_ART_PROTOCOL, filesep, 'artifact_', ARTIFACT_CRIT , '_params.txt'],'FileType','text');
  param_default = param(param.subject == "default",:);
  param_subject = param(strcmp(param.subject,SUBJECT),:);
  if ~isempty(param_subject)
    param = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%% Applying High Pass and Notch Filters
%   skip for vibration-denoised data
cfg=[];
cfg.hpfilter=high_pass_filter;
cfg.hpfreq=high_pass_filter_freq;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.bsfilter=do_bsfilter;
cfg.bsfreq= [line_noise_harm_freqs-1; line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
D_hpf_trials = ft_preprocessing(cfg,D_trials);

%% Artifact rejection - ECoG channels 
% iterating over bands and electrode types
artifact = table();
hfig = figure();
for idx = 1:height(param)
  
  fprintf('doing %s %s \n',SUBJECT,param.name{idx});

    cfg = [];
    cfg.out_freq = sample_rate; 
    cfg.suppress_output = 1; 
    cfg.param = param(idx,:); 
    
    D_avgpow_eltype_trials = multifreq_avg_power(cfg, D_hpf_trials);

    cfg.runs = session; 
    D_avgpow_eltype = concat_trials_within_run(cfg, D_avgpow_eltype_trials); 
    
    if isempty(D_avgpow_eltype_trials)
        %channel type not available
        continue
    end

    % concatenate trials within each run so that we can detect artifacts from the continuous (wavelet transformed) data
    D_avgpow_eltype = concat_trials_within_run(cfg, D_avgpow_eltype_trials); 

  ENVELOPE_BIN_SIZE_SECONDS = param.env_bin_size(idx); %envelope bin size in seconds
  THRESHOLD_STD_FACTORS = [param.th_factor_std_low(idx), param.th_factor_std_high(idx)]; %factors to determine detection thresholds 
  THRESHOLD_FIX = [param.th_fix_min(idx), param.th_fix_max(idx)]; %fix thresholds to filter data before applying robust estimates
  CONSOLIDATION_TIME_TOLERANCE = param.th_consolidation(idx); %min time allowed between consecutive artifacts
  ELECTRODE_COVERAGE_THRESHOLD = param.th_frac_coverage(idx); %max allowed fraction of time with artifacts
  CONNECTOR_THRESHOLD = [param.th_conn_low(idx), param.th_conn_high(idx)]; %detection threshold for number of electrodes in a connector  

  cfg=[];
  cfg.freq = 1/ENVELOPE_BIN_SIZE_SECONDS;
  D_avgpow_eltype_env = bml_envelope_binabs(cfg,D_avgpow_eltype);

  %calculating log10 transform (envelopes have log normal distributions)
  D_hg_eltype_env_log10 = bml_apply(@(x) param.env_mult_factor(idx) .* log10(x),D_avgpow_eltype_env);
  
  cfg=[];
  cfg.remask_inf=true;
  cfg.value=NaN;
  D_hg_eltype_env_log10 = bml_mask(cfg,D_hg_eltype_env_log10);
  
  %calculating distribution robust statistics. 
  THRESHOLD = nan(ntrials,2);
  max_v=nan(1,ntrials);
  min_v=nan(1,ntrials);
  for itrial=1:ntrials
    v = reshape(D_hg_eltype_env_log10.trial{itrial},1,[]);
    v1 = v((v>THRESHOLD_FIX(1)) & (v<THRESHOLD_FIX(2)));
    m = median(v1);
    std = bml_robust_std(v1);
    if ~isempty(v1)
      max_v(itrial)=max(v);
      min_v(itrial)=min(v);
      THRESHOLD(itrial,:) = m + THRESHOLD_STD_FACTORS.*std;
    end
  end

  %plotting histogram to assess threshold levels
  if ~isempty(v1)
      clf(hfig); set(hfig,'Position',[0 0 600 600]);
      for itrial=1:ntrials
        subplot(ceil(ntrials/2),2,itrial)
        hold on;
        h=histogram(D_hg_eltype_env_log10.trial{itrial},linspace(min(min_v),max(max_v),61),...
          'FaceAlpha',0.1,'EdgeAlpha',1);
        maxBinCount = max(h.BinCounts);
        plot([THRESHOLD(itrial,1),THRESHOLD(itrial,1)],[0,maxBinCount .* 1.1]);
        plot([THRESHOLD(itrial,2),THRESHOLD(itrial,2)],[0,maxBinCount .* 1.1]);
        %set(gca,'YScale','log')
        title(['session ' num2str(itrial)]);
      end
      %saveas(hfig,[PATH_FIGURES filesep SUBJECT '_' pname '_artifact_env_log10_hist.png'])
  elseif isempty(v1)
      warning(['For electrode type ''', param.electrode_type{idx}, ''' (sub ', SUBJECT, '), no timepoints found between low threshold (', ...
          num2str(THRESHOLD_FIX(1)), ') and high threshold (', num2str(THRESHOLD_FIX(2)), ')'])
  end

  %detecting segments of time for each channel above threshold
  artifact_eltype_1 = table();
  for itrial=1:ntrials
    cfg=[];
    cfg.threshold = THRESHOLD(itrial,:);
    cfg.trials = itrial;
    artifact_eltype_1 = bml_annot_rowbind(artifact_eltype_1, bml_annot_detect(cfg,D_hg_eltype_env_log10));
  end

  if isempty(artifact_eltype_1)
    continue
  end
  
  %making figure with random snippets of detected artifacts
  cfg=[];
  cfg.n = 1;
  cfg.groupby  = 'label';
  artifact_eltype_1_sample = bml_annot_sample(cfg, artifact_eltype_1);
  artifact_eltype_1_sample = bml_annot_extend(artifact_eltype_1_sample,2,2);

  cfg=[];
  cfg.n = 60;
  artifact_eltype_1_sample = bml_annot_sample(cfg, artifact_eltype_1_sample);
  
  cfg=[];
  cfg.epoch = artifact_eltype_1_sample;
  [D_hpf_eltype_sample, epoch_hpf_eltype_sample] = bml_redefinetrial(cfg,D_hpf_eltype);
 
  D_p = D_hpf_eltype_sample;
  E_p = epoch_hpf_eltype_sample;
  nx=10; ny=floor(numel(D_p.trial)/nx);
  if ny==0
    ny=1; nx=numel(D_p.trial);
  end
  clf(hfig); set(hfig,'Position',[0 0 nx*200 ny*200]);
  for itrial=1:ny
      for j=1:nx
          pidx = (itrial-1)*nx+j;
          l = E_p.label(pidx);
          l_idx = bml_getidx(l,D_p.label);
          subplot(ny,nx,pidx);
          plot(D_p.time{pidx},D_p.trial{pidx}(l_idx,:));
          title(E_p.label(pidx));
      end
  end
  %saveas(hfig,[PATH_FIGURES filesep SUBJECT '_' pname '_artifact_snippets.png'])

  
  %consolidating annotations with CONSOLIDATION_TIME_TOLERANCE margin of overlap
  cfg=[];
  cfg.criterion = @(x) (x.starts(end) - max(x.ends(1:(end-1))) < CONSOLIDATION_TIME_TOLERANCE);
  cfg.groupby = 'label';
  artifact_eltype_2 = bml_annot_consolidate(cfg,artifact_eltype_1);

%   %creating ft_raw from annotations for visualization
%   cfg=[];
%   cfg.template = D_hpf_eltype_env_log10;
%   cfg.annot_label_colname='label';
%   artifact_eltype_3_raw = bml_annot2raw(cfg,artifact_eltype_2);
% 
%   %raster plot of artifacts for session 1
%   f=figure();
%   bml_plot_raster(artifact_eltype_3_raw)

  %check if excluded segments are correct
%   cfg=[];
%   cfg.label_colname = 'label';
%   cfg.annot = artifact_eltype_2;
%   cfg.value = NaN;
%   D_hpf_eltype_mask = bml_mask(cfg, D_hpf_eltype);
% 
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';
%   cfg.continuous = 'yes';
%   ft_databrowser(cfg,D_hpf_eltype_mask);

  %% rejecting faulty channels 

  %decide which artifacts to include. Usually just ECoG artifacts
  %artifact = bml_annot_rowbind(artifact_ecog_3,artifact_macro_3,artifact_dbs_3);
  artifact_1 = artifact_eltype_2;


  cfg = [];
  cfg.groupby = 'label';
  artifact_1_session_cvg = bml_annot_coverage(cfg,artifact_1,session);

  %histogram(artifact_1_session_cvg.coverage,linspace(0,1,51))

  %if a channel in a session has more than COVERAGE_THRESHOLD of the time with
  %artifacts, the entire channel gets rejected for that session

  artifact_1_session_cvg_sel = artifact_1_session_cvg(artifact_1_session_cvg.coverage >= ELECTRODE_COVERAGE_THRESHOLD,:);
  artifact_2 = bml_annot_rowbind(artifact_1,artifact_1_session_cvg_sel);
  cfg=[];
  cfg.groupby = 'label';
  artifact_2 = bml_annot_consolidate(cfg,artifact_2);

%   %creating ft_raw from annotations for visualization
%   cfg=[];
%   cfg.template = D_hpf_env;
%   cfg.annot_label_colname='label';
%   artifact2_raw = bml_annot2raw(cfg,artifact_2);
% 
%   %raster plot of artifacts for session 1
%   f=figure();
%   bml_plot_raster(artifact2_raw)

  %% checking coverage per connector group
  %if several channels of the same connector group have an artifact, reject
  %the entire connector group

  %adding connector information to artifac annotation table
  electrodes.conn_label = strcat({'conn'},num2str(electrodes.connector));
  artifact_2.conn_label = bml_map(artifact_2.label, electrodes.name, electrodes.conn_label);

% 	%calculating absolute value envelope at 1Hz (1s chunks)
%   cfg=[];
%   cfg.freq=ENVELOPE_BIN_SIZE_SECONDS;
%   D_hpf_env = bml_envelope_binabs(cfg,D_hpf);

  %for each connector and bin, count number of faulty channels
  cfg=[];
  cfg.roi = bml_raw2annot(D_avgpow_eltype_env);
  cfg.annot_label_colname = 'conn_label';
  connector_artifact_2_cvg_raw = bml_annot2raw(cfg,artifact_2);

%   f=figure();
%   cfg.colorbar = true;
%   bml_plot_raster(cfg,connector_artifact_2_cvg_raw)

  %detecting faulty connectors
  cfg=[];
  cfg.threshold = CONNECTOR_THRESHOLD;
  connector_artifact_3 = bml_annot_detect(cfg,connector_artifact_2_cvg_raw);

  if ~isempty(connector_artifact_3)
    %for each period a connector is faulty, create table with all channels
    %corresponding to that connctor
    cfg=[];
    cfg.groupby_x='conn_label'; %grouping variable in electrode table
    cfg.groupby_y='label'; %corresponding grouping variable in connector_artifact_3
    artifact_4=bml_annot_intersect(cfg,electrodes,connector_artifact_3);

    if ~isempty(artifact_4)
      %combining with previously detected artifacts
      artifact_4.label = artifact_4.name;
      artifact_5 = bml_annot_rowbind(artifact_2, artifact_4);
      cfg=[];
      cfg.groupby = 'label';
      artifact_5 = bml_annot_consolidate(cfg,artifact_5);
    else
      artifact_5 = artifact_2;
    end
  else
    artifact_5 = artifact_2;
  end

  %final raster plot for artifacts
  cfg=[];
  cfg.template = D_avgpow_eltype_env;
  cfg.annot_label_colname = 'label';
  artifact_5_raw = bml_annot2raw(cfg,artifact_5);

  clf(hfig); set(hfig,'Position',[0 0 600 600]);
  cfg.trial_name='session';
  bml_plot_raster(cfg,artifact_5_raw)
  saveas(hfig,[PATH_FIGURES filesep SUBJECT '_' pname '_artifact_mask.png'])

  artifact_5.pname = repmat({pname},height(artifact_5),1);
  
%% saving  artifact annotation table
  artifact = bml_annot_rowbind(artifact,...
    artifact_5(:,{'id','starts','ends','duration','label','conn_label','pname'}));
  
end




%archiving 
if isfile(FILEPATH_ARTPROT_ANNOT)
  copyfile(FILEPATH_ARTPROT_ANNOT,...
        [PATH_ANNOT filesep 'archive' filesep SUBJECT '_artifact-criteria-' ARTIFACT_CRIT '_' datestr(now,'yyyymmdd_HHMM') '.txt'])
end

writetable(artifact,FILEPATH_ARTPROT_ANNOT, 'delimiter','\t', 'FileType','text');


