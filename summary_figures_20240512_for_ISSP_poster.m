%%%% plot summary figures for presentation
% make a 2x3 panel plot which shows proportions of electrodes across areas tuned for a parameter...
% .... along w/ an illustrative ctx and subctx example timecourse
% .... top row = beta, bot row = high gamma
%%%% examples are intended to highlight general trends of responsivity/selectivity in different regions and freqbands
% 
% if resp_hg and resp_beta are not yet loaded, they will be loaded
%%% make sure and compare_areal_tuning_triplet and plot_resp_timecourse_triplet aren't clearing/resetting important params

%% selection one of the following option sets for different figures


%%%%%%% other hg subcort examples which show peak before each syllable onset: 
%%%%%%% ...... , DBS4077_dbs_L1;    'DBS4083','dbs_R3';   'DBS3019','macro_c'
%%% beta subcort alternatives: 'DBS4077','dbs_L2C'
% plotop.elc_param = 'p_prep'; 
% plotop.trial_sort_param = [];
% plotop.sub_elc = {  {'DBS3028','ecog_158'},{'DBS3017','macro_c'};... % beta
%                     {'DBS3004','ecog_135'}, {'DBS3019','macro_c'} };  % hg
% plotop.x_ax_hardlims = [-3 2]; 
% 
% 
%%%% ecog hg alternates: DBS3014_ecog_136
% plotop.elc_param = 'p_prod';
% plotop.trial_sort_param = [];
% plotop.sub_elc = { {'DBS3029','ecog_113'},{'DBS3028','macro_m'};...  % beta
%                     {'DBS3021','ecog_161'},{'DBS3018','dbs_L4'} }; % hg
% plotop.x_ax_hardlims = [-3 2]; 

%%%%%%%% ecog beta alternates: DBS3001___ecog_111, DBS3014___ecog_137, DBS4083___ecog_116
%%%%%%%% ecog hg alternates: DBS3027___ecog_119
plotop.elc_param = {'p_prep_cons',1};
plotop.trial_sort_param = {'cons',1};
plotop.sub_elc = {{'DBS3021','ecog_123'}, {'DBS3001','macro_c'};... % beta....
                    {'DBS3021','ecog_133'},{'DBS3011','dbs_L4'} }; % hg
plotop.x_ax_hardlims = [-3 1.7]; 

% plotop.elc_param = ;
% plotop.trial_sort_param = [];
% plotop.sub_elc = { {'DBS',''},{'DBS',''};... % beta
%                     {'DBS',''},{'DBS',''} }; % hg


%% general params



set(0,'DefaultFigureWindowStyle','normal')

plotop.ncols = 3; 
smooth_timecourses = 1; 
    smooth_windowsize = 20; 
%      smooth_method = 'movmean'; 
     smooth_method = 'gaussian';



%%
setpaths_dbs_triplet()

%%% load resp tables for hg and beta first
if ~exist('resp_beta','var')
    load([PATH_RESULTS, filesep, 'resp_all-subjects_beta_ar-F_ref-CTAR_denoised.mat']); resp_beta = resp;
end
if ~exist('resp_hg','var')
    load([PATH_RESULTS, filesep, 'resp_all-subjects_hg_ar-E_ref-CTAR_denoised.mat']); resp_hg = resp;
end

op.newfig = 0; 
global op subs htl

hfig = figure('Position', get(0, 'Screensize')+[20 50 -50 -150]);
htl = tiledlayout(2, plotop.ncols, 'Padding','tight');

%%%%%%% top row - beta
resp = resp_beta; op.art_crit = 'F'; op.resp_signal = 'beta';
plotop.tile_row = 1; 
make_figure(plotop, resp)

%%%% bot row - hg
resp = resp_hg; op.art_crit = 'E';  op.resp_signal = 'hg';
plotop.tile_row = 2; 
make_figure(plotop, resp)

print -dbitmap % copy figure to clipboard

%% function for plotting one row of panels
function make_figure(plotop, resp)

global op subs htl

srt = resp; 
sort_cond = plotop.trial_sort_param;
newfig = 0; 

%%%% areal proportion comparison
tile_index = (plotop.tile_row - 1) * plotop.ncols + 1; % Calculate the tile index
nexttile
param = plotop.elc_param;
compare_areal_tuning_triplet();

%%% example timecourse 1
tile_index = (plotop.tile_row - 1) * plotop.ncols + 2; % Calculate the tile index
hax = nexttile(htl, tile_index);
sub_elc = plotop.sub_elc{ plotop.tile_row, 1};
thissub = sub_elc{1};
thiselc = sub_elc{2};
srt_row = find(strcmp(srt.sub, thissub) & strcmp(srt.chan, thiselc));
plot_resp_timecourse_triplet()

%%% example timecourse 2
tile_index = (plotop.tile_row - 1) * plotop.ncols + 3; % Calculate the tile index
hax = nexttile(htl, tile_index);
sub_elc = plotop.sub_elc{ plotop.tile_row, 2};
thissub = sub_elc{1};
thiselc = sub_elc{2};
srt_row = find(strcmp(srt.sub, thissub) & strcmp(srt.chan, thiselc));
plot_resp_timecourse_triplet()

end
