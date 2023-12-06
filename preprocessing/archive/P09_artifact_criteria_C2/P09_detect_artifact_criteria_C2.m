function P09_detect_artifact_criteria_C2(SUBJECT, param)

%% Detect artifacts
% creates an artifact annotation table

% CRITERIA C paramaters valus
CRITERIA = 'C2'; %identifier for the criteria implemented in this script

HIGH_PASS_FILTER = 'yes'; %should a hisgh pass filter be applied
HIGH_PASS_FILTER_FREQ = 1; %cutoff frequency of high pass filter

%   Parameters now defined in artifact_C_params
%   ENVELOPE_BIN_SIZE_SECONDS = 1; %envelope bin size in seconds
%   ECOG_THRESHOLD_STD_FACTORS = [1.5, 3]; %factors to determine detection thresholds 
%   ECOG_CONSOLIDATION_TIME_TOLERANCE = 3; %min time allowed between consecutive artifacts
%   ELECTRODE_COVERAGE_THRESHOLD = 0.5; %max allowed fraction of time with artifacts
%   CONNECTOR_THRESHOLD = [4, 8]; %detection threshold for number of electrodes in a connector

%% load packages
%addpath('E:\MATLAB\fieldtrip')
%addpath('E:\MATLAB\NPMK')
%addpath('E:\MATLAB\bml')

ft_defaults
bml_defaults
format long

%% Definig paths
%SUBJECT='DBS3022';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='Z:\DBS';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_PROTOCOL = 'C:\Users\amsmeier\Documents\MATLAB\P09_artifact_criteria_C2';

cd(PATH_SYNC)

session= bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);



%% Loading FieldTrip data 
load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_raw_session.mat'],'D','loaded_epoch');
nTrials = numel(D.trial);

%remasking nans with zeros
cfg=[];
cfg.value=0;
cfg.remask_nan=true;
D=bml_mask(cfg,D);

%% working in protocol folder
cd(PATH_PROTOCOL)

%% loading electrode type band table
if ~exist('el_band','var')
  param = readtable('artifact_C2_params.txt');
  param_default = param(param.subject == "default",:);
  param_subject = param(strcmp(param.subject,SUBJECT),:);
  if ~isempty(param_subject)
    param = bml_annot_rowbind(param_default(~ismember(param_default.name,param_subject.name),:),param_subject);
  end
end

%% Applying High Pass Filter
cfg=[];
cfg.hpfilter=HIGH_PASS_FILTER;
cfg.hpfreq=HIGH_PASS_FILTER_FREQ;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.channel={'ecog_*','macro_*','micro_*','dbs_*'};
D_hpf = ft_preprocessing(cfg,D);

%% Artifact rejection - ECoG channels 
% iterating over bands and electrode types
artifact = table();
f = figure();
for idx = 1:height(param)
  
  fprintf('doing %s %s \n',SUBJECT,param.name{idx});
  
  el_type = strip(param.electrode_type{idx});
  wav_freq = param.wav_freq(idx);
	wav_width = param.wav_width(idx);
  env_mult_factor =  param.env_mult_factor(idx);
  pname = strip(param.name{idx});
  
  ENVELOPE_BIN_SIZE_SECONDS = param.env_bin_size(idx); %envelope bin size in seconds
  THRESHOLD_STD_FACTORS = [param.th_factor_std_low(idx), param.th_factor_std_high(idx)]; %factors to determine detection thresholds 
  THRESHOLD_FIX = [param.th_fix_min(idx), param.th_fix_max(idx)]; %fix thresholds to filter data before applying robust estimates
  CONSOLIDATION_TIME_TOLERANCE = param.th_consolidation(idx); %min time allowed between consecutive artifacts
  ELECTRODE_COVERAGE_THRESHOLD = param.th_frac_coverage(idx); %max allowed fraction of time with artifacts
  CONNECTOR_THRESHOLD = [param.th_conn_low(idx), param.th_conn_high(idx)]; %detection threshold for number of electrodes in a connector  
  
  %selecting ECoG channels for artifact rejection
  cfg=[];
  cfg.channel = [el_type,'_*'];
  D_hpf_eltype = ft_selectdata(cfg,D_hpf);

  if isempty(D_hpf_eltype.label)
    %channel type not available
    continue
  end
  
%   %viasually inspect signal
%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';4
%   cfg.continuous = 'yes';
%   cfg.channel = {'ecog_11*', 'audio_*'};
%   ft_databrowser(cfg,D);

  %calculating absolute value envelope at 1Hz (1s chunks)
  cfg=[];
  cfg.out_freq = 100;
  cfg.wav_freq = wav_freq;
  cfg.wav_width = wav_width;
  D_hpf_eltype_env = bml_envelope_wavpow(cfg,D_hpf_eltype);

%   cfg=[];
%   cfg.viewmode = 'vertical';
%   cfg.blocksize = 30;
%   cfg.ylim = 'maxmin';
%   cfg.continuous = 'yes';
%   ft_databrowser(cfg,D_hpf_ecog_env);

  cfg=[];
  cfg.freq = 1/ENVELOPE_BIN_SIZE_SECONDS;
  D_hpf_eltype_env = bml_envelope_binabs(cfg,D_hpf_eltype_env);

  %calculating log10 transform (envelopes have log normal distributions)
  D_hpf_eltype_env_log10 = bml_apply(@(x) env_mult_factor .* log10(x),D_hpf_eltype_env);
  
  cfg=[];
  cfg.remask_inf=true;
  cfg.value=NaN;
  D_hpf_eltype_env_log10 = bml_mask(cfg,D_hpf_eltype_env_log10);
  
  %calculating distribution robust statistics. 
  THRESHOLD = nan(nTrials,2);
  max_v=nan(1,nTrials);
  min_v=nan(1,nTrials);
  for i=1:nTrials
    v = reshape(D_hpf_eltype_env_log10.trial{i},1,[]);
    v1 = v((v>THRESHOLD_FIX(1)) & (v<THRESHOLD_FIX(2)));
    m = median(v1);
    std = bml_robust_std(v1);
    if ~isempty(v1)
      max_v(i)=max(v);
      min_v(i)=min(v);
      THRESHOLD(i,:) = m + THRESHOLD_STD_FACTORS.*std;
    end
  end

  %plotting histogram to asses threshold levels
  clf(f); set(f,'Position',[0 0 600 600]);
  for i=1:nTrials
    subplot(ceil(nTrials/2),2,i)
    hold on;
    h=histogram(D_hpf_eltype_env_log10.trial{i},linspace(min(min_v),max(max_v),61),...
      'FaceAlpha',0.1,'EdgeAlpha',1);
    maxBinCount = max(h.BinCounts);
    plot([THRESHOLD(i,1),THRESHOLD(i,1)],[0,maxBinCount .* 1.1]);
    plot([THRESHOLD(i,2),THRESHOLD(i,2)],[0,maxBinCount .* 1.1]);
    %set(gca,'YScale','log')
    title(['session ' num2str(i)]);
  end
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_env_log10_hist.png'])

  %detecting segments of time for each channel above threshold
  artifact_eltype_1 = table();
  for i=1:nTrials
    cfg=[];
    cfg.threshold = THRESHOLD(i,:);
    cfg.trials = i;
    artifact_eltype_1 = bml_annot_rowbind(artifact_eltype_1, bml_annot_detect(cfg,D_hpf_eltype_env_log10));
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
  clf(f); set(f,'Position',[0 0 nx*200 ny*200]);
  for i=1:ny
      for j=1:nx
          pidx = (i-1)*nx+j;
          l = E_p.label(pidx);
          l_idx = bml_getidx(l,D_p.label);
          subplot(ny,nx,pidx);
          plot(D_p.time{pidx},D_p.trial{pidx}(l_idx,:));
          title(E_p.label(pidx));
      end
  end
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_snippets.png'])

  
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

  %% cheking coverage per connector group
  %if several channels of the same connector group have an artifact, reject
  %the entire connector group

  %adding connector information to artifac annotation table
  electrode.conn_label = strcat({'conn'},num2str(electrode.connector));
  artifact_2.conn_label = bml_map(artifact_2.label, electrode.electrode, electrode.conn_label);

% 	%calculating absolute value envelope at 1Hz (1s chunks)
%   cfg=[];
%   cfg.freq=ENVELOPE_BIN_SIZE_SECONDS;
%   D_hpf_env = bml_envelope_binabs(cfg,D_hpf);

  %for each connector and bin, count number of faulty channels
  cfg=[];
  cfg.roi = bml_raw2annot(D_hpf_eltype_env);
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
    artifact_4=bml_annot_intersect(cfg,electrode,connector_artifact_3);

    if ~isempty(artifact_4)
      %combining with previously detected artifacts
      artifact_4.label = artifact_4.electrode;
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
  cfg.template = D_hpf_eltype_env;
  cfg.annot_label_colname = 'label';
  artifact_5_raw = bml_annot2raw(cfg,artifact_5);

  clf(f); set(f,'Position',[0 0 600 600]);
  cfg.trial_name='session';
  bml_plot_raster(cfg,artifact_5_raw)
  saveas(f,['figures/' SUBJECT '_' pname '_artifact_mask.png'])

  artifact_5.pname = repmat({pname},height(artifact_5),1);
  
%% saving  artifact annotation table
  artifact = bml_annot_rowbind(artifact,...
    artifact_5(:,{'id','starts','ends','duration','label','conn_label','pname'}));
  
end

cd(PATH_SYNC)
artifact_annot_path = ['annot/' SUBJECT '_artifact_criteria_' CRITERIA '.txt'];


% AM commented out archiving because permissions for sync and annot folders are not available

% % % % %archiving 
% % % % if isfile(artifact_annot_path)
% % % %   copyfile(artifact_annot_path,...
% % % %         [PATH_SYNC filesep 'archive' filesep SUBJECT '_artifact_criteria_' CRITERIA '_' datestr(now,'yyyymmdd_HHMM') '.txt'])
% % % % end





% AM changed write-to folder because permissions for sync and annot folders are not available

% bml_annot_write(artifact,['annot/' SUBJECT '_artifact_criteria_' CRITERIA '.txt']);
bml_annot_write(artifact,[PATH_PROTOCOL, filesep, 'annot', filesep, SUBJECT '_artifact_criteria_' CRITERIA '.txt']);

