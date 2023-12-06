PATH_DATASET = 'Z:\DBS';
PATH_ANALYSIS = [PATH_DATASET '\Analysis\triplet_analysis_am']; 

electrodes_vars_to_add = {'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'};

subs = readtable('C:\Users\amsmeier\Documents\MATLAB\P09_Subjects_to_analyze.txt'); 
nsubs = height(subs);

resp_temp = table;

for isub = 1:nsubs
    SUBJECT = subs.sub{isub}
    PATH_SYNC = [PATH_DATASET filesep SUBJECT filesep 'Preprocessed Data\Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];
    electrodes_file = [PATH_ANNOT filesep SUBJECT '_electrode.txt'];
    try 
        loadfile = ['Z:\DBS\Analysis\triplet_analysis_am\' SUBJECT '_responses'];
        load(loadfile)
        electrodes_sub = readtable(electrodes_file);

        for ielc = 1:height(resp)
            rowmatch = find(strcmp(electrodes_sub.electrode, resp.chan{ielc}), 1); % sometimes matches multiple rows of electrode.txt ......
            resp{ielc,electrodes_vars_to_add} = electrodes_sub{rowmatch, electrodes_vars_to_add};
        end

        resp_temp = [resp_temp; resp];
        subs.trials{isub} = trials; 
    end
end

resp = resp_temp; clear resp_temp
save([PATH_ANALYSIS, filesep, 'resp_all_subjects'],'resp','subs')
