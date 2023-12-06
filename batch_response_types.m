 
PATH_ANALYSIS = 'Z:\DBS\Analysis\triplet_analysis_am'; 

subject_list = readtable('C:\Users\amsmeier\Documents\MATLAB\P09_Subjects_to_analyze.txt'); 
nsubs = height(subject_list);

for isub = 1:nsubs
    SUBJECT = subject_list.subject{isub}
    try 
        clearvars -except SUBJECT subject_list isub nsubs PATH_ANALYSIS
        response_types()
        savefile = [PATH_ANALYSIS, filesep, SUBJECT '_responses'];
        save(savefile, 'trials', 'resp')
    end
end
