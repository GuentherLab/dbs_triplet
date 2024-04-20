 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas
  % load resp_all_subjects first
setpaths_dbs_triplet()

% close all

 %% params
vardefault('show_barplot',1);

analyze_responsive_elcs_only = 1; 
newfig = 1; 
warn_about_unassigned_elcs = 0; 

%%% define anatomical regions composed of smaller areas
regiondef = {   'mfg',  {'rostralmiddlefrontal' , 'caudalmiddlefrontal'};... middle frontal gyrus... maybe also inf front sulcus
                'ifg',  {'parsopercularis', 'parsorbitalis',  'parstriangularis'};... % inferior frontal gyrus
                'smc',  {'postcentral', 'precentral'};...                   % sensorimotor cortex
                'suptemp', {'superiortemporal', 'bankssts' }; ... % superior temporal
                % 'thal', {'Left-Thalamus-Proper' , 'Right-Thalamus-Proper' };... % thalamus.... use MOREL labels instead
                % 'wm',   {'Left-Cerebral-White-Matter', 'Right-Cerebral-White-Matter'}; ... % white matter.... use MOREL labels instead
                % 'ventdc', {'Left-VentralDC', 'Right-VentralDC'}... %% ? not sure what VentralDC is.... use MOREL labels instead
                 % 'sfg'   {'superiorfrontal'};... % superior frontal gyrus..... should be excluded - only has 2 electrodes
                 'stn', {'STh_L', 'STh_L'};...
                 'thal', {'VApc_L','VLa_L' , 'VLpv_L', 'VM_L' 'VM_R', 'VPM_L'};...
%                  'gp', {'GPe_L','GPi_sensorimotor_L'}; % <20 electrodes, so maybe not worth including
                };

% param = 'p_rank'; 
% param = 'p_prep';

% param = {'p_stim_cons',1};
% param = {'p_stim_cons',2};
% param = {'p_stim_cons',3};
% param = {'p_stim_syl',1};
% param = {'p_stim_syl',2};
% param = {'p_stim_syl',3};
% param = {'p_stim_vow',1};
% param = {'p_stim_vow',2};
% param = {'p_stim_vow',3};

% param = {'p_prep_cons',1};
param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};

% param = 'p_stim_cons_allpos';
% param = 'p_stim_vow_allpos';
% param = 'p_stim_syl_allpos';

% param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; 
% param = 'p_prep_syl_constit';

% param = 'p_prod_cons_allpos';
% param = 'p_prod_vow_allpos';
% param = 'p_prod_syl_allpos';


alpha = 0.05; 

bar_face_color = [0.5 0.5 0.5]; 


[paramvals, param_name, full_param_string] = triplet_tablevar(resp, param); 

compare_areal_tuning() % in ieeg_funcs_am



