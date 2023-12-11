 %%% for each 3000 subject, check that it has an annot table with times of each stim syllable/phoneme onset and offset

 PATH_PROTOCOL = 'Z:\DBS\Batch\P08_artifact_criteria_E';
 PROTOCOL_TABLE = [PATH_PROTOCOL filesep 'P08_Subjects_3000.txt']; 
PATH_DATA = 'Z:\DBS';
save_filepath = ['Z:\DBS\Analysis\triplet_results_am\archive\subjects_missing_stimsyl_table_', datestr(now,'yyyy-mm-dd'), '.txt'];

subject_table = readtable(PROTOCOL_TABLE);
nsubs = height(subject_table);  
subject_table.has_stimsyl_table = nan(nsubs,1);
subject_table.has_stimphon_table = nan(nsubs,1);


for isub = 1:nsubs
  SUBJECT = subject_table.subject{isub};
  PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
  PATH_PREPROCESSED=[PATH_SUBJECT filesep 'Preprocessed Data'];
  PATH_SYNC=[PATH_PREPROCESSED filesep 'Sync'];
  PATH_ANNOT = [PATH_SYNC filesep 'annot']; 

  if exist([PATH_ANNOT, filesep, SUBJECT, '_stimulus_syllable.txt'],'file')
      subject_table.has_stimsyl_table(isub) = 1;
  elseif ~exist([PATH_ANNOT, filesep, SUBJECT, '_stimulus_syllable.txt'],'file')
    subject_table.has_stimsyl_table(isub) = 0;
  end

  if exist([PATH_ANNOT, filesep, SUBJECT, '_stimulus_phoneme.txt'],'file')
      subject_table.has_stimphon_table(isub) = 1;
  elseif ~exist([PATH_ANNOT, filesep, SUBJECT, '_stimulus_phoneme.txt'],'file')
    subject_table.has_stimphon_table(isub) = 0;
  end

end    

writetable(subject_table, save_filepath, 'delimiter', '\t')