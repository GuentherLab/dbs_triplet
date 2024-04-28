%%%% generate and save list of subjects to run

clear 
setpaths_dbs_triplet()

% location in which save the new subject list
subj_list_savename = {[PATH_TRIPLET_CODE filesep 'P08_artifact_criteria_E' filesep 'P08_Subjects_to_analyze.txt'];
                        [PATH_TRIPLET_CODE filesep 'P08_artifact_criteria_F' filesep 'P08_Subjects_to_analyze.txt']}; 

subs_to_exclude = {...
%     'DBS3004';... % got stuck during P08 preprocessing, maybe artifact detection... might be fixable; was not a problem with HG preprocessing scripts used before 2024/4/20
%     'DBS3006';... % got stuck during P08 preprocessing, maybe artifact detection... might be fixable; was not a problem with HG preprocessing scripts used before 2024/4/20
%     'DBS3011';... % got stuck during P08 preprocessing, maybe artifact detection [block 2]... might be fixable; was not a problem with HG preprocessing scripts used before 2024/4/20
    'DBS3025';... % causes error in P09_highgamma_from_denoised_rereferenced.m during a specific trial; also had hearing issues - see DBS3000-Sessions and intraop notes
    'DBS4060';... % missing stim timing info
    'DBS4062';... % missing stim timing info
    'DBS4068';... % electrodes table is missing mni info
    'DBS4069';... % "Missed first trial, very short session, couldn't stay awake"... only 18 trials
    'DBS4071';... % got stuck during P08 preprocessing, maybe artifact detection... might be fixable
    'DBS4079';... % got stuck during P08 preprocessing, maybe artifact detection... might be fixable
    'DBS4081';... % didn't do triplet task because couldn't hear stim; see intraop notes
    'DBS4082';... % aborted task because patient starting vomiting
};

pitt_subject_list_file = 'Z:\DBS\DBS_subject_lists\REDCap-20191112.xlsx';

subs_pitt_all = readtable(pitt_subject_list_file);
subs_triplet = subs_pitt_all(string(subs_pitt_all.Task)=="U01 speech",:);
subs_triplet = renamevars(subs_triplet, 'ID', 'subject');
subs_triplet(contains(subs_triplet.subject, subs_to_exclude),:) = []; % remove specific subjects

npaths = length(subj_list_savename);
for ipath = 1:npaths
    writetable(subs_triplet, subj_list_savename{ipath}, 'Delimiter', 'tab')
end
