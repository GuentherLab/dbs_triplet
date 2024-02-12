
% andrew meier

clear

PATH_PROTOCOL = 'Z:\DBS\Batch\P08_artifact_criteria_E';
PROTOCOL_FUNCTION = 'P08_detect_artifact_criteria_E';
PROTOCOL_TABLE = 'P08_Subjects_3000.txt';
ARTIFACT_CRIT = 'E'; 
exe_daytime = datestr(now,'yyyymmdd_HHMM');
addpath(PATH_PROTOCOL);
diary([PATH_PROTOCOL filesep 'batch_' PROTOCOL_FUNCTION '_' exe_daytime '.log'])

SKIP_OK = false; %Should previous protocols run by this script successfully be skipped
FORCE = true; %Archive all previous versions of the script and run current 
              %overrides any manual modification

PATH_DATA = 'Z:\DBS';
cd(PATH_PROTOCOL)

subject_table = readtable(PROTOCOL_TABLE);         
fprintf('=== Running protocol %s ===\n',PROTOCOL_FUNCTION)
if FORCE; fprintf('Forced run, overwritting any manual change.\n'); end
if SKIP_OK; fprintf('Skipping previously successfully executed protocols.\n'); end

sub_inds_to_run = 28:height(subject_table);
% sub_inds_to_run = [4];

%% subject loop
for isub = sub_inds_to_run
  SUBJECT = subject_table.subject{isub};
  timetol = subject_table.timetol(isub);
  
  PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
  PATH_PREPROCESSED=[PATH_SUBJECT filesep 'Preprocessed Data'];
  PATH_SYNC=[PATH_PREPROCESSED filesep 'Sync'];
  PATH_ARCHIVE=[PATH_SYNC filesep 'archive'];
  fprintf('\n* Subject %s: ',SUBJECT)
  
    rereferenced_ft_file = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip'...
        filesep SUBJECT '_ft_raw_filt_trial_denoised_ref_criteria_' ARTIFACT_CRIT '.mat'];

    % name of the file to save trialwise referenced HG responses into
    rereferenced_highgamma_savename = [PATH_SUBJECT '/Preprocessed Data/FieldTrip' filesep SUBJECT...
        '_ft_hg_trial_ref_criteria_' ARTIFACT_CRIT '_denoised.mat'];

	new_proto_fname = [SUBJECT '_' PROTOCOL_FUNCTION '_' exe_daytime]; 
  old_proto = dir([PATH_SYNC filesep SUBJECT '_' PROTOCOL_FUNCTION '*.m']);
  
       
    %running protocol
        PATH_DATA='Z:\DBS';
        PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
        PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
        PATH_ANNOT = [PATH_SYNC '/annot']; 
        cd(PATH_SYNC)
        prodtrip = readtable([PATH_ANNOT filesep SUBJECT '_produced_triplet.txt']);
        if ~any("onset" == string(prodtrip.Properties.VariableNames))
            prodtrip.onset = prodtrip.starts;
        end

%         if ~any(isnan(prodtrip.onset))
            SUBJECT

            close all force
            
            P08A09_highgamma_from_denoised(SUBJECT,timetol); 
    
            P08A09_detect_artifact_criteria_E(SUBJECT);
    
            P09_redefine_trial_common_average_reference_denoised(SUBJECT);
    
            P09_highgamma_from_denoised_rereferenced(SUBJECT,rereferenced_ft_file,ARTIFACT_CRIT,rereferenced_highgamma_savename);

%         end

end
diary('off')
