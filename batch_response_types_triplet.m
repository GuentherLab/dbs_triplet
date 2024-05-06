% get triplet electrode response profiles for a list of subjects 
%%%% optionally, also put all electrodes from analyzed subjects into a single table
%%%% ........ and put trial table from all subjects into a single variable

clear

%% params
% op.art_crit = 'E'; op.resp_signal = 'hg';
op.art_crit = 'F'; op.resp_signal = 'beta';

op.denoised = 1; % work with vibration-denoised data

op.rereference_method = 'none';
% op.rereference_method = 'CTAR';

op.out_freq = 100; % freq of wavpow output files

%%%%%%%%%%%%%%%%%%%%%%%
setpaths_dbs_triplet(); set_project_specific_variables(); % need to run after setting art crit
%%%%%%%%%%%%%%%%%%

subject_list_filename = [PATH_ARTIFACT filesep 'P08_Subjects_to_analyze.txt']; % created by generate_triplet_subject_list.m
% subject_list_filename = [PATH_ARTIFACT filesep 'P08_Subjects_3000.txt'];

subs = readtable(subject_list_filename); nsubs = height(subs);

subinds_to_run = [1:nsubs]; % all subjects
% subinds_to_run = []; % 
% subinds_to_run = [28:nsubs]; % 4000-series subs

%% get responses types for individual subjects
analysis_spec_string = [op.resp_signal, '_ar-', op.art_crit, '_ref-',op.rereference_method, op.denoise_string]; 
subs_resp_dir =  [PATH_RESULTS, filesep, 'subs_', analysis_spec_string];
if ~exist("subs_resp_dir", 'dir')
    mkdir(subs_resp_dir)
end

tic
for isub = subinds_to_run
    clearvars -except op subs isub nsubs subinds_to_run subject_list_filename subs_resp_dir analysis_spec_string
    
    op.sub = subs.subject{isub};
    fprintf(['.... Getting response types (', analysis_spec_string, ') for subject: ' op.sub, '\n'])
       
    response_types_triplet()
    savefile = [subs_resp_dir, filesep, op.sub '_responses_' analysis_spec_string];
    save(savefile, 'trials', 'resp', 'op')
end
toc


%% compile response tables
compiled_responses_filepath = [PATH_RESULTS, filesep, 'resp_all-subjects_', analysis_spec_string]; 
fprintf('Compiling response tables in %s \n', compiled_responses_filepath);
subs = readtable(subject_list_filename);
nsubs = height(subs);
resp_temp = table;

tic
for isub = subinds_to_run
    op.sub = subs.subject{isub};
    PATH_SYNC = [PATH_DATA filesep op.sub filesep 'Preprocessed Data\Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];

%     try 
        loadfile = [subs_resp_dir, filesep, op.sub '_responses_' analysis_spec_string];
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
op = rmfield(op,'sub'); % save non-subject-specific version
save(compiled_responses_filepath,'resp','subs','op',  '-v7.3')

toc