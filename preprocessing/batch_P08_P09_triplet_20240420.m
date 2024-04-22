
%%%% run artifact detection, rereferencing, and computation of signal of interest (e.g. high gamma) for triplet subjects

%%% if we want to include time tolerance data for each subject, it is available in : 
%%% ...... Z:\DBS\DBS_subject_lists\time-tolerance.tsv

clear


%% params
op.art_crit = 'E'; op.resp_signal = 'hg';
% op.art_crit = 'F'; op.resp_signal = 'beta';

op.denoised = 1; % work with vibration-denoised data

setpaths_dbs_triplet() % need to run after setting art crit

% op.rereference_method = 'none';
op.rereference_method = 'CTAR';

op.out_freq = 100; % freq of wavpow output files

% % % % PROTOCOL_TABLE = 'P08_Subjects_3000.txt';
PROTOCOL_TABLE = [PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt'];  % created by generate_triplet_subject_list.m

subject_table = readtable(PROTOCOL_TABLE);         

sub_inds_to_run = 1:height(subject_table);
% sub_inds_to_run = [1 2];



%% subject loop

for isub = sub_inds_to_run
  op.sub = subject_table.subject{isub};
  set_project_specific_variables(); % set paths etc. based on data collection site, load timing and electrode data
  fprintf('\n* Preprocessing subject %s...',op.sub)
    cd(PATH_SYNC)
    prodtrip = readtable([PATH_ANNOT filesep op.sub '_produced_triplet.txt']);
    if ~any("onset" == string(prodtrip.Properties.VariableNames))
        prodtrip.onset = prodtrip.starts;
    end

    close all force
    
    P08A09_wavpow_from_denoised(op); %%% compute and save wav power

    P08A09_detect_artifact_denoised(op);

    P09_redefine_trial_common_avg_ref_denoised(op);

    P09_wavpow_from_denoised_rereferenced(op); %%% check that this produces identical outputs as wavpow version, then remove/archive it


end
diary('off')
