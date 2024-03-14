%%%% generate and save list of subjects to run

pitt_subject_list_file = 'Z:\DBS\DBS_subject_lists\REDCap-20191112.xlsx';
subs_pitt_all = readtable(pitt_subject_list_file);
subs_triplet = subs_pitt_all(string(subs_pitt_all.Task)=="U01 speech",:);


