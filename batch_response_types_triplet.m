% get triplet electrode response profiles for a list of subjects 

setpaths_dbs_triplet()

subject_list = readtable([PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']); % TURBO path
% % % % % subject_list = readtable([PATH_ARTIFACT filesep 'P08_Subjects_3000.txt']);
nsubs = height(subject_list);

subinds_to_run = [3:21, 23:28]; % 3000-series; subs causing errors
subinds_to_run = [29:nsubs]; % 4000-series subs

for isub = subinds_to_run
    SUBJECT = subject_list.subject{isub}
        clearvars -except SUBJECT subject_list isub nsubs PATH_RESULTS subinds_to_run

        response_types()
        savefile = [PATH_RESULTS, filesep, SUBJECT '_responses'];
        save(savefile, 'trials', 'resp')
end

compile_resp_subjects()