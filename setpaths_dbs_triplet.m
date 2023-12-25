%%%% set paths for AM triplet analysis depending on computer

[~,compname] = system('hostname'); compname = string(deblank(compname));


 switch compname
     case 'MSI' % AM personal computer
         PATH_DATA='D:\triplet'; 
         PATH_RESULTS = [PATH_DATA filesep 'triplet_results_am']; 
         PATH_ARTIFACT = [PATH_DATA filesep 'P08_artifact_criteria_E']; 
         PATH_CODE = 'C:\docs\code'; % AM laptop top directory for all code repos
         PATH_TRIPLET_CODE = [PATH_CODE filesep 'dbs_triplet']; 
         PATH_BML = [PATH_CODE filesep 'bml']; 
         PATH_FIELDTRIP_CODE = [PATH_CODE filesep 'fieldtrip-20210616']; 
     case 'NSSBML01' % TURBO - BML server computer
         PATH_DATA='Z:\DBS';
         PATH_RESULTS = [PATH DATA filesep '\Analysis\triplet_results_am'];
         PATH_ARTIFACT = [PATH_DATA filesep 'Batch\P08_artifact_criteria_E']; 
         PATH_TRIPLET_CODE = 'C:\Users\amsmeier\dbs_triplet'; 
         PATH_BML = 'C:\Program Files\Brain-Modulation-Lab\bml'; 
         PATH_FIELDTRIP_CODE = 'Y:\Users\lbullock\MATLAB_external_libs_Turbo20230907\fieldtrip'; 
     otherwise 
         error('computer name not recognized; please add computer to setpaths_dbs_triplet.m')
 end

paths_to_add = {PATH_DATA;... % derivatives and (if on server) sourcedata
                PATH_RESULTS;... % outputs of post-derivatives analyses by AM
                PATH_ARTIFACT;... % artifact param definition; also may contain subject artifact tables, though these are also stored in sub annot folders
                PATH_TRIPLET_CODE;... % code by AM for triplet analysis
                    [PATH_TRIPLET_CODE filesep 'preprocessing'];...
                    [PATH_TRIPLET_CODE filesep 'util'];  
                PATH_BML;... % Brain Modulation Lab repo
                PATH_FIELDTRIP_CODE;...
    };
addpath(paths_to_add{:});

ft_defaults()
bml_defaults()

set(0, 'DefaultTextInterpreter', 'none')
set(0, 'DefaultLegendInterpreter', 'none')
% set(0, 'DefaultAxesTickLabelInterpreter', 'none')

format long

 clearvars compname paths_to_add PATH_CODE PATH_FIELDTRIP_CODE PATH_BML