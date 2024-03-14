%%%% generate and save list of subjects to run

% location in which save the new subject list
subj_list_savename = [PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']; 

pitt_subject_list_file = 'Z:\DBS\DBS_subject_lists\REDCap-20191112.xlsx';

subs_pitt_all = readtable(pitt_subject_list_file);
subs_triplet = subs_pitt_all(string(subs_pitt_all.Task)=="U01 speech",:);
subs_triplet = renamevars(subs_triplet, 'ID', 'subject');

writetable(subs_triplet, subj_list_savename, 'Delimiter', 'tab')
