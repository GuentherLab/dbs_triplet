% get triplet electrode response profiles for a list of subjects 
%%%% optionally, also put all electrodes from analyzed subjects into a single table
%%%% ........ and put trial table from all subjects into a single variable

setpaths_dbs_triplet()

%% params
subject_list_filename = [PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']; % created by generate_triplet_subject_list.m
% subject_list_filename = [PATH_ARTIFACT filesep 'P08_Subjects_3000.txt'];

compiled_responses_filepath = [PATH_RESULTS, filesep, 'resp_all_subjects']; 




subs = readtable(subject_list_filename);
nsubs = height(subs);

subinds_to_run = [1:nsubs]; % all subjects
% subinds_to_run = []; % 
% subinds_to_run = [28:nsubs]; % 4000-series subs

%% get responses types for individual subjects

for isub = subinds_to_run
    clearvars -except subs isub nsubs subinds_to_run compiled_responses_filepath subject_list_filename
    
    op.sub = subs.subject{isub};
    fprintf([op.sub, '\n'])
       
    response_types_triplet()
    savefile = [PATH_RESULTS, filesep, op.sub '_responses'];
    save(savefile, 'trials', 'resp')
end



%% compile response tables
fprintf(['Compiling response tables in %s \n'], compiled_responses_filepath);
subs = readtable(subject_list_filename);
nsubs = height(subs);
resp_temp = table;

for isub = subinds_to_run
    op.sub = subs.subject{isub};
    PATH_SYNC = [PATH_DATA filesep op.sub filesep 'Preprocessed Data\Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];

%     try 
        loadfile = [PATH_RESULTS filesep op.sub '_responses.mat'];
        load(loadfile)

        %%% the problem addressed here might be fixable by looking at electrods.txt, but 'port' isn't an important variable for our analysis 
        if isnumeric(resp.port) 
            resp.port = num2cell(resp.port);
        end

        resp_temp = [resp_temp; resp];
        subs.trials{isub} = trials; 
%     end
end

resp = resp_temp; clear resp_temp
save(compiled_responses_filepath,'resp','subs',  '-v7.3')