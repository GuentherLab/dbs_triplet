 %%% for each 3000 subject, check that the number of trials listed in the trial_epoch table.....
 %%%%% .... matches that in the vibration-denoised fieldtrip file

 PATH_PROTOCOL = 'Z:\DBS\Batch\P08_artifact_criteria_E';
PATH_DATA = 'Z:\DBS';
save_filepath = 'Z:\DBS\Analysis\triplet_results_am\archive\trialnumber_mismatches_in_denoised.txt';

subject_table = readtable(PROTOCOL_TABLE);
nsubs = height(subject_table);  
subject_table.ntrials_epoch_table = nan(nsubs,1);
subject_table.ntrials_fieldtrip_denoised = nan(nsubs,1);


for isub = 1:nsubs
    try 
  SUBJECT = subject_table.subject{isub}
  PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
  PATH_PREPROCESSED=[PATH_SUBJECT filesep 'Preprocessed Data'];
  PATH_FIELDTRIP = [PATH_PREPROCESSED filesep 'FieldTrip']; 
  PATH_SYNC=[PATH_PREPROCESSED filesep 'Sync'];
  PATH_ANNOT = [PATH_SYNC filesep 'annot']; 
  PATH_ARCHIVE=[PATH_SYNC filesep 'archive'];

  clear D trials
  load([PATH_FIELDTRIP filesep SUBJECT '_ft_raw_filt_trial_denoised.mat'],'D')
  trials = readtable([PATH_ANNOT filesep SUBJECT '_trial_epoch.txt']);

  subject_table.ntrials_epoch_table(isub) = height(trials);
  subject_table.ntrials_fieldtrip_denoised(isub) = length(D.trial); 
    end

end    

subject_table.ntrials_missing_in_fieldtrip = subject_table.ntrials_epoch_table - subject_table.ntrials_fieldtrip_denoised;
subject_table.trialnumber_mismatch = subject_table.ntrials_missing_in_fieldtrip ~= 0; 

writeable(subject_table, save_filepath)