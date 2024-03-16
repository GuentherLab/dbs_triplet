%%%% generate and save list of subjects to run

clear 
setpaths_dbs_triplet()

% location in which save the new subject list
subj_list_savename = [PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']; 

subs_to_exclude = {...
    'DBS4060';... % missing stim timing info
    'DBS4062';... % missing stim timing info
    'DBS4069';... % "Missed first trial, very short session, couldn't stay awake"... only 18 trials
    'DBS4071';... % got stuck during preprocessing, maybe artifact detection... might be fixable
    'DBS4079';... % got stuck during preprocessing, maybe artifact detection... might be fixable
    'DBS4082';... % zero trials
};

pitt_subject_list_file = 'Z:\DBS\DBS_subject_lists\REDCap-20191112.xlsx';

subs_pitt_all = readtable(pitt_subject_list_file);
subs_triplet = subs_pitt_all(string(subs_pitt_all.Task)=="U01 speech",:);
subs_triplet = renamevars(subs_triplet, 'ID', 'subject');
subs_triplet(contains(subs_triplet.subject, subs_to_exclude),:) = []; % remove specific subjects

writetable(subs_triplet, subj_list_savename, 'Delimiter', 'tab')
