 
PATH_RESULTS = 'Z:\DBS\Analysis\triplet_results_am'; 

% subject_list = readtable('C:\Users\amsmeier\Documents\MATLAB\P08_Subjects_to_analyze.txt'); 
subject_list = readtable('Z:\DBS\Batch\P08_artifact_criteria_E\P08_Subjects_3000.txt');
nsubs = height(subject_list);

subinds_to_run = [3:21, 23:nsubs];

for isub = subinds_to_run
    SUBJECT = subject_list.subject{isub}
%     try 
        clearvars -except SUBJECT subject_list isub nsubs PATH_RESULTS subinds_to_run
        response_types()
        savefile = [PATH_RESULTS, filesep, SUBJECT '_responses'];
        save(savefile, 'trials', 'resp')
%     end
end

compile_resp_subjects()