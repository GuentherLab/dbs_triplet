
% put all electrodes from analyzed subjects into a single table


subs = readtable([PATH_ARTIFACT filesep 'P08_Subjects_3000.txt']);
nsubs = height(subs);

subinds_to_run = [3:21, 23:nsubs];

resp_temp = table;

for isub = subinds_to_run
    SUBJECT = subs.subject{isub}
    PATH_SYNC = [PATH_DATA filesep SUBJECT filesep 'Preprocessed Data\Sync'];
    PATH_ANNOT = [PATH_SYNC filesep 'annot'];

%     try 
        loadfile = [PATH_RESULTS filesep SUBJECT '_responses.mat'];
        load(loadfile)
% % % % % %             electrodes_file = [PATH_ANNOT filesep SUBJECT '_electrode.txt'];
% % % % % %         electrodes_sub = readtable(electrodes_file);

        if isnumeric(resp.port)
            resp.port = num2cell(resp.port);
        end
        if isnumeric(resp.MOREL_label_3)
            resp.MOREL_label_3 = num2cell(resp.MOREL_label_3);
        end

        resp_temp = [resp_temp; resp];
        subs.trials{isub} = trials; 
%     end
end

resp = resp_temp; clear resp_temp
save([PATH_RESULTS, filesep, 'resp_all_subjects'],'resp','subs',  '-v7.3')
