 %%%% plot the timecourses of electrode's responses for different phonemes/syllables
 % load resp_all_subjects first
 %
 % align all trialwise responses to speech onset of syllable #1
 %
 % % % updated by AM 2022/8/21
 
close all

%% params
srt_row = 35;
smooth_timecourses = 1; 
    smooth_windowsize = 30; 
%      smooth_method = 'movmean'; 
     smooth_method = 'gaussian';
show_error_bars = 0; 

y_ax_hardlims = []; % cut off y axis if it's lesser/greater than this value
% y_ax_hardlims = [-1 4]; % cut off y axis if it's lesser/greater than this value
% y_ax_hardlims = [-4 10]; % cut off y axis if it's lesser/greater than this value

xlimits = [-3 2]; 

plotops.linewidth = 2; 

condval_inds_to_plot = []; % plot all vals
% condval_inds_to_plot = [1 4 7 10]; 
% condval_inds_to_plot = [2 5 8 11]; 
% condval_inds_to_plot = [3 6 9 12]; 
% condval_inds_to_plot = [1:3]; 
% condval_inds_to_plot = [1:6];
% condval_inds_to_plot = [7:12]; 
% condval_inds_to_plot = [1:12]; 

%%% choose the stimulus variable which will be used to sort trials
% sort_cond = []; % average all trials in a single trace
% sort_cond = 'stim_volume'; 
% sort_cond = {'cons',1};
% sort_cond = {'cons',2};
% sort_cond = {'cons',3};
% sort_cond = {'vow',1};
% sort_cond = {'vow',2};
% sort_cond = {'vow',3};
% sort_cond = {'syl',1}; 
% sort_cond = {'syl',2}; 
% sort_cond = {'syl',3}; 
% sort_cond = 'cons_constit';
% sort_cond = 'vow_constit'; 
sort_cond = 'syl_constit';
 
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

plot_timecourses = 1; 
plot_raster = 0; 

xline_color_stim_syl_on = [0.3 0.3 0.8];
xline_color_stim_syl_off = [0.8 0.3 0.3];
xline_color_prod_syl_on = [0.7 0.7 0.7];
xline_color_prod_syl_off = [0.3 0.3 0.3];
xline_style = '--';
xline_width = 0.25; 

%%%%% method for finding time landmarks from trial times
% xline_fn = @mean; 
xline_fn = @median;

yline_zero_width = 0.25; 
yline_zero_color = [0.8 0.8 0.8]; 
yline_zero_style = '-';

%%%% how to find the time length that trials will be cut/padded to be
% trial_time_adj_method = 'median';
% trial_time_adj_method = 'median_plus_sd'; % median plus stdev
trial_time_adj_method = 'max';



%% align responses
channame = srt.chan{srt_row};
thissub = srt.sub(srt_row,:);

subrow = find(subs.subject == string(thissub));
trials = subs.trials{subrow}; 


ntrials = height(trials);
 % align resp_allonses to first syl onset
 nans_tr = nan(ntrials,1); 
 trials_tmp = trials; % temporary copy of trials table
 


%%
% times relative to produced syllable #1 onset
trials_tmp.stim_syl_on_adj = trials_tmp.stim_syl_on - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.stim_syl_off_adj = trials_tmp.stim_syl_off - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.prod_syl_on_adj = trials_tmp.prod_syl_on - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.prod_syl_off_adj = trials_tmp.prod_syl_off - trials_tmp.prod_syl_on(:,1) ; 

%% organize responses by grouping var
% find trial details for the appropriate subject
if isempty(sort_cond) % plot all trials in a single trace
    trial_conds = true(ntrials,1); 
    full_var_string = 'all_trials';
elseif ~isempty(sort_cond)
    [trial_conds, ~, full_var_string] = triplet_tablevar(trials,sort_cond); 
end

 %%
 plot_response_timecourse()

 %%

htitle = title([thissub, '___', channame, '... ', full_var_string], 'Interpreter','none');
%     hleg = legend(resp_grpd.condval{condval_inds_to_plot});

f=get(gca,'Children');

if ~isempty(sort_cond)
    hleg = legend(flipud(f(end-nvals_to_plot+1:end)),resp_grpd.condval{condval_inds_to_plot});
end

% stim syllable onsets
hstim_syl(1) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,1),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off');
hstim_syl(2) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off');
hstim_syl(3) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off');

% stim offset
hstim_syl(3) = xline(mean(trials_tmp.stim_syl_off_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_off, 'LineStyle',xline_style, 'HandleVisibility','off');

% produced syllable onsets
hprod_syl(1) = xline(0, 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off'); 
hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_on_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off');
hprod_syl(3) = xline(xline_fn(trials_tmp.prod_syl_on_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style, 'HandleVisibility','off');

% produced syllable offsets
hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,1),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style, 'HandleVisibility','off');
hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style, 'HandleVisibility','off');
hprod_syl(3) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style, 'HandleVisibility','off');

hyline = yline(0, 'LineWidth',yline_zero_width, 'Color',yline_zero_color, 'LineStyle',yline_zero_style, 'HandleVisibility','off');

xlabel('Time (sec)')
ylabel('HG power (normed)')
