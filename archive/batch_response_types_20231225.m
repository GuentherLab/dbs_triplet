% get triplet electrode response profiles for a list of subjects 

%%%% this version includes extra checks on whether the trialwise rerefenced denoised high gammma fieldtrip file exists

setpaths_dbs_triplet()

% subject_list = readtable('C:\Users\amsmeier\Documents\MATLAB\P08_Subjects_to_analyze.txt'); % TURBO path
subject_list = readtable([PATH_ARTIFACT filesep 'P08_Subjects_3000.txt']);
nsubs = height(subject_list);

subinds_to_run = [3:21, 23:nsubs];





no_hg_trial_denoised = {};



for isub = subinds_to_run
    SUBJECT = subject_list.subject{isub}
        clearvars -except PATH_DATA SUBJECT subject_list isub nsubs PATH_RESULTS subinds_to_run no_hg_trial_denoised





    PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
    PATH_PREPROCESSED = [PATH_SUBJECT filesep 'Preprocessed Data'];
    PATH_FIELDTRIP = [PATH_PREPROCESSED filesep 'FieldTrip']; 
    PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];
    ARTIFACT_CRIT = 'E'; 
    if ~exist([PATH_FIELDTRIP filesep SUBJECT '_ft_hg_trial_ref_criteria_' ARTIFACT_CRIT '_denoised.mat'])
        no_hg_trial_denoised = [no_hg_trial_denoised; SUBJECT]
        continue
    end














        response_types()
        savefile = [PATH_RESULTS, filesep, SUBJECT '_responses'];
        save(savefile, 'trials', 'resp')
end

compile_resp_subjects()