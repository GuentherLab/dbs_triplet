%%%% set paths for AM triplet analysis depending on computer

[~,compname] = system('hostname'); compname = string(deblank(compname));

 switch compname
     case {'MSI','677-GUE-WL-0010'} % AM personal computer, work laptop
         PATH_DATA='D:\triplet'; %%% may not have a copy on all machines.... 600GB 
         PATH_RESULTS = ['C:\triplet_results_am']; %%% use the SSD for faster load times
         PATH_CODE = 'C:\docs\code'; % AM laptop top directory for all code repos
         PATH_TRIPLET_CODE = [PATH_CODE filesep 'dbs_triplet']; 
         PATH_IEEG_FT_FUNCS_AM = [PATH_CODE filesep 'ieeg_ft_funcs_am']; % ieeg processing code shared across AM projects
         PATH_BML = [PATH_CODE filesep 'bml']; 
         PATH_FIELDTRIP_CODE = [PATH_CODE filesep 'fieldtrip']; 
         PATH_LEADDBS = [PATH_CODE filesep ]; 
         PATH_AVERAGE_MNI = [PATH_RESULTS filesep 'atlases' filesep 'CortexLowRes_15000V_MNI_ICBM_2009b_NLIN_ASYM.mat']; 
         PATH_SUBCORT_ATLAS_VIM = [PATH_RESULTS filesep 'atlases' filesep 'atlas_index_DISTAL_Ewert2017.mat']; % too large to keep in github; put atlases in Results
         PATH_STN_ATLAS = [PATH_RESULTS filesep 'atlases' filesep 'atlas_index_subcort_Ewert_v2.1.7.mat']; 
     case 'NSSBML01' % TURBO - BML server computer
         PATH_DATA='Z:\DBS';
         PATH_RESULTS = [PATH_DATA filesep 'Analysis\triplet_results_am'];
         PATH_TRIPLET_CODE = 'C:\Users\amsmeier\dbs_triplet'; 
             % PATH_ARTIFACT = [PATH_DATA filesep 'Batch\P08_artifact_criteria_E']; %%%% periodically should copy from triplet repo to this folder
         PATH_IEEG_FT_FUNCS_AM = 'C:\Users\amsmeier\ieeg_ft_funcs_am'; % ieeg processing code shared across AM projects
         PATH_BML = 'C:\Program Files\Brain-Modulation-Lab\bml'; 
         PATH_FIELDTRIP_CODE = 'Y:\Users\lbullock\MATLAB_external_libs_Turbo20230907\fieldtrip'; 
         PATH_LEADDBS = 'C:\Program Files\LeadDBS';
         PATH_AVERAGE_MNI = 'Z:/DBS/DBS_subject_lists/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexLowRes_15000V.mat';
         PATH_SUBCORT_ATLAS_VIM = 'C:\Program Files\LeadDBS_Classic\leaddbs\templates\space\MNI_ICBM_2009b_NLIN_ASYM\atlases\DISTAL (Ewert 2017)/atlas_index.mat';
         PATH_STN_ATLAS = 'Z:\Resources\STN-Atlas\atlas_index.mat';
     otherwise 
         error('computer name not recognized; please add computer to setpaths_dbs_triplet.m')
 end

addpath([PATH_IEEG_FT_FUNCS_AM, filesep, 'util']); % add first so we can use vardefault and field_default

vardefault('op',struct);
field_default('op','art_crit','E'); % E = high gamma, F = beta

 % common paths

PATH_ARTIFACT = [PATH_TRIPLET_CODE filesep 'P08_artifact_criteria_' op.art_crit]; % keep in repo to sync across devices
PATH_STIM_INFO = [PATH_TRIPLET_CODE filesep 'stim_info']; 

paths_to_add = {PATH_DATA;... % derivatives and (if on server) sourcedata
                PATH_RESULTS;... % outputs of post-derivatives analyses by AM
                PATH_ARTIFACT;... % artifact param definition; also may contain subject artifact tables, though these are also stored in sub annot folders
                PATH_STIM_INFO;... % information about average stim timing; used to fill in info for subjects in which it is missing
                PATH_TRIPLET_CODE;... % code by AM for triplet analysis
                    [PATH_TRIPLET_CODE filesep 'preprocessing'];...
                    [PATH_TRIPLET_CODE filesep 'daphne_scripts'];  ... % 2024 rotation project code
                    [PATH_TRIPLET_CODE filesep 'util_triplet']
                PATH_IEEG_FT_FUNCS_AM;.... 
                PATH_BML;... % Brain Modulation Lab repo
                PATH_FIELDTRIP_CODE;...
                PATH_LEADDBS;...
    };
addpath(paths_to_add{:});

ft_defaults()
bml_defaults()

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
% set(0, 'DefaultAxesTickLabelInterpreter', 'none')

format long

 clearvars compname paths_to_add PATH_CODE PATH_BML