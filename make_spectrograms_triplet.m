%%% create spectrams for each trial - triplet data
%

setpaths_dbs_triplet()

%%
op.sub = 'DBS3012'; 
op.base_win_sec = [1, 0.3]; 
op.post_speech_win_sec = 0.5; % time to include after 3rd-syllable voice offset


%% get stim and behavioral timing, align fieldtrip timing
set_project_specific_variables()
load_triplet_stim_beh_timing()


%%
% make_spectrograms()