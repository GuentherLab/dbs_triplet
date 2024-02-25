%%
%Original Code
% get triplet electrode response profiles for a list of subjects 

setpaths_dbs_triplet()

% subject_list = readtable('C:\Users\amsmeier\Documents\MATLAB\P08_Subjects_to_analyze.txt'); % TURBO path
subject_list = readtable([PATH_ARTIFACT filesep 'P08_Subjects_3000.txt']);
nsubs = height(subject_list);

subinds_to_run = [3:21, 23:nsubs];

for isub = subinds_to_run
    SUBJECT = subject_list.subject{isub}
        clearvars -except SUBJECT subject_list isub nsubs PATH_RESULTS subinds_to_run

        response_types()
        savefile = [PATH_RESULTS, filesep, SUBJECT '_responses'];
        save(savefile, 'trials', 'resp')
end

compile_resp_subjects()














%% DT Addition
% Initialize a summary table or structure
summary_results = struct();

% Loop through each subject file
for isub = subinds_to_run
    subject_id = subject_list.subject{isub};
    loadfile = [PATH_RESULTS filesep subject_id '_responses.mat'];
    if isfile(loadfile)
        load(loadfile, 'resp');
        
        % Extract the ANOVA results
        summary_results.(subject_id).p_phonotactic_prob = resp.p_phonotactic_prob;
        summary_results.(subject_id).p_trans_id = resp.p_trans_id;
    else
        warning('File for %s does not exist.', subject_id);
    end
end

% Now 'summary_results' contains the ANOVA results for all subjects
