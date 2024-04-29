%%%% plot summary figures for presentation
% need to load resp tables for hg and beta first

%% option sets for different figures

%%%% 'DBS4083','dbs_R2B' shows peak before each syllable onset
plotop.elc_param = 'p_prep'; 
plotop.trial_sort_param = [];
plotop.sub_elc = {  {'DBS3010','ecog_222'},{'DBS3011','dbs_L4'};... % beta
                    {'DBS3020','ecog_212'}, {'DBS4083','dbs_R2B'} };  % hg


plotop.elc_param = 'p_prod';
plotop.trial_sort_param = [];
plotop_sub_elc = { {'DBS3029','ecog_113'},{'DBS3028','macro_m'};...  % beta
                    {'DBS3021','ecog_161'},{'DBS3018','dbs_L4'} }; % hg




% plotop.elc_param = ;
% plotop.trial_sort_param = [];
% plotop_sub_elc = { {'DBS',''},{'DBS',''};... % beta
%                     {'DBS',''},{'DBS',''} }; % hg


%% general params
smooth_timecourses = 1; 
    smooth_windowsize = 20; 
%      smooth_method = 'movmean'; 
     smooth_method = 'gaussian';

%%
op.newfig = 0; 
global op subs htl

htl = tiledlayout(2, 3);

%%%%%%% top row - beta
resp = resp_beta; op.resp_signal = 'beta';
plotop.tile_rowcol(1) = 2; 


%%%% bot row - hg
resp = resp_hg; op.resp_signal = 'hg';
plotop.tile_row = 2; 
make_figure(plotop, resp)


%%
function make_figure(plotop, resp)

global op subs htl

%%%% areal proportion comparison
nexttile(tilenum(htl, plotop.tile_row, 1);
param = plotop.elc_param;
compare_areal_tuning();

%%% example timecourse 1
sub_elc = plotop.sub_elc{ plotop.tile_rowcol(1), plotop.tile_rowcol(2)};
thissub = sub_elc{1};
thiselc = sub_elc{2};
srt_row = find(strcmp(resp.subject, thissub) & strcmp(resp.chan, thiselc));
plot_resp_timecourse_triplet()

%%% example timecourse 2
sub_elc = plotop.sub_elc{ plotop.tile_rowcol(1), plotop.tile_rowcol(2)};
thissub = sub_elc{1};
thiselc = sub_elc{2};
srt_row = find(strcmp(resp.subject, thissub) & strcmp(resp.chan, thiselc));
plot_resp_timecourse_triplet()
