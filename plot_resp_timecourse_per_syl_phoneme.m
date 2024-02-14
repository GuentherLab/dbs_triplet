 %%%% plot the timecourses of electrode's responses for different phonemes/syllables
 % load resp_all_subjects first
 %
 % align all trialwise responses to speech onset of syllable #1
 %
 % % % updated by AM 2022/8/21
 
% close all

%% params
srt_row = 5;
show_error_bars = 0; 

y_ax_hardlims = []; % cut off y axis if it's lesser/greater than this value
% y_ax_hardlims = [-1 4]; % cut off y axis if it's lesser/greater than this value
% y_ax_hardlims = [-4 10]; % cut off y axis if it's lesser/greater than this value

xlimits = [-3 2]; 

plotops.linewidth = 2; 

groupval_inds_to_plot = []; % plot all vals
% groupval_inds_to_plot = [1 4 7 10]; 
% groupval_inds_to_plot = [2 5 8 11]; 
% groupval_inds_to_plot = [3 6 9 12]; 
% groupval_inds_to_plot = [1:3]; 
% groupval_inds_to_plot = [1:6];
% groupval_inds_to_plot = [7:12]; 
% groupval_inds_to_plot = [1:12]; 

%%% choose the stimulus variable which will be used to sort trials
trial_grouping_var = {'cons',1};
% trial_grouping_var = {'cons',2};
% trial_grouping_var = {'cons',3};
% trial_grouping_var = {'vow',1};
% trial_grouping_var = {'vow',2};
% trial_grouping_var = {'vow',3};
% trial_grouping_var = {'syl',1}; 
% trial_grouping_var = {'syl',2}; 
% trial_grouping_var = {'syl',3}; 
% trial_grouping_var = 'cons_constit';
% trial_grouping_var = 'vow_constit'; 
% trial_grouping_var = 'syl_constit';
 
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
 trials_tmp = [trials_tmp, table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
 
 samp_period = 1e-5 * round(1e5 * diff(trials.times{1}(1:2))); % sampling interval
 
 %%% find trial lengths pre- and post-onset
    %%%% trials.times{itrial} use global time coordinates
    %%%% ....... start at a fixed baseline window before stim onset
    %%%% ....... end at a fixed time buffer after speech offset
for itrial = 1:ntrials
    % n timepoints before or at voice onset [time from beginning of pre-stim baseline to voice onset]
    trials_tmp.tpoints_pre_onset(itrial) = nnz(trials_tmp.times{itrial} <= trials_tmp.prod_syl_on(itrial,1)); 
    % n timepoints after voice onset [time from voice onset to end of post-voice time buffer]
    trials_tmp.tpoints_post_onset(itrial) = nnz(trials_tmp.times{itrial} > trials_tmp.prod_syl_on(itrial,1)); 
end
 
% pad or cut each trial to fit a specific size, so that we can align and average trials
switch trial_time_adj_method
    case 'median_plus_sd'
        n_tpoints_pre_fixed = round(median(trials_tmp.tpoints_pre_onset) + std(trials_tmp.tpoints_pre_onset)); 
        n_tpoints_post_fixed = round(median(trials_tmp.tpoints_post_onset) + std(trials_tmp.tpoints_post_onset)); 
    case 'median'
        n_tpoints_pre_fixed = median(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = median(trials_tmp.tpoints_post_onset); 
    case 'max'
        n_tpoints_pre_fixed = max(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = max(trials_tmp.tpoints_post_onset); 
end
tpoints_tot = n_tpoints_pre_fixed + n_tpoints_post_fixed; 

% nan-pad or cut trial windows so that they are all the same duration
%%% pad and cut values must be non-negative
resp_align = struct; 
resp_align.resp = NaN(ntrials, tpoints_tot); % aligned responses for this electrode; rows = trials, columns = timepoints
for itrial = 1:ntrials
   pre_pad = max([0, n_tpoints_pre_fixed - trials_tmp.tpoints_pre_onset(itrial)]); 
   pre_cut = max([0, -n_tpoints_pre_fixed + trials_tmp.tpoints_pre_onset(itrial)]); 
   % inds from resp.timecourse... if pre_cut > 0, some timepoints from this trial will not be used
   pre_inds = 1+pre_cut:trials_tmp.tpoints_pre_onset(itrial); 
   % fill in pre-onset data... fill in electrode responses starting after the padding epoch
   resp_align.resp(itrial, pre_pad+1 : n_tpoints_pre_fixed) = srt.timecourse{srt_row}{itrial}(pre_inds); 

   post_pad = max([0, n_tpoints_post_fixed - trials_tmp.tpoints_post_onset(itrial)]);
   post_cut = max([0, -n_tpoints_post_fixed + trials_tmp.tpoints_post_onset(itrial)]); 
   post_inds = trials_tmp.tpoints_pre_onset(itrial) +  [1 : trials_tmp.tpoints_post_onset(itrial)-post_cut]; % inds from resp.timecourse
   resp_align.resp(itrial, n_tpoints_pre_fixed+1:end-post_pad) = srt.timecourse{srt_row}{itrial}(post_inds); % fill in post-onset data
   
   trials_tmp.trial_onset_adjust(itrial) = samp_period * [pre_pad - pre_cut]; % number of timepoints to add to time landmarks
end
resp_align.mean = mean(resp_align.resp,'omitnan'); % mean response timecourse
resp_align.std = std(resp_align.resp, 'omitnan'); % stdev of response timecourses
resp_align.std_lims = [resp_align.mean + resp_align.std; resp_align.mean - resp_align.std]; 
resp_align.n_nonnan_tpoints = sum(~isnan(resp_align.resp));
resp_align.sem = resp_align.std ./ sqrt(resp_align.n_nonnan_tpoints);
resp_align.sem_lims = [resp_align.mean + resp_align.sem; resp_align.mean - resp_align.sem]; 

% times relative to produced syllable #1 onset
trials_tmp.stim_syl_on_adj = trials_tmp.stim_syl_on - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.stim_syl_off_adj = trials_tmp.stim_syl_off - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.prod_syl_on_adj = trials_tmp.prod_syl_on - trials_tmp.prod_syl_on(:,1) ; 
trials_tmp.prod_syl_off_adj = trials_tmp.prod_syl_off - trials_tmp.prod_syl_on(:,1) ; 

%% organize responses by grouping var
% find trial details for the appropriate subject
[triptabvals, ~, full_var_string] = triplet_tablevar(trials,trial_grouping_var); 
[unq_grouping_vals, ~, grouping_var_ind] = unique( triptabvals );
ngroupvals = length(unq_grouping_vals);
celcol = cell(ngroupvals,1);
resp_grpd = table(unq_grouping_vals,celcol,celcol,'VariableNames',{'groupval','resp','resp_mean'}); 
for igroupval = 1:ngroupvals
    these_trial_inds = grouping_var_ind == igroupval;
    resp_grpd.resp{igroupval} = resp_align.resp(these_trial_inds,:);
    resp_grpd.resp_mean{igroupval} = mean(resp_grpd.resp{igroupval},1,'omitnan');
end

%% plotting
xtime = 0.5 + [linspace(-n_tpoints_pre_fixed, -1, n_tpoints_pre_fixed), linspace(0, n_tpoints_post_fixed-1, n_tpoints_post_fixed)];
xtime = samp_period * xtime; 


if plot_timecourses 

hfig = figure('Color',[1 1 1]); 
hold off

% if grouping val indices not specified, plot them all
if isempty (groupval_inds_to_plot)
    groupval_inds_to_plot = 1:ngroupvals;
end
nvals_to_plot = length(groupval_inds_to_plot); 

% error bars
if show_error_bars
    hfill = fill([xtime, fliplr(xtime)], [resp_align.sem_lims(1,:), fliplr(resp_align.sem_lims(2,:))], [0.8 0.8 0.8]); % standard error
        hfill.LineStyle = 'none'; % no border
    hold on 
end

hplot = plot(xtime,cell2mat(resp_grpd.resp_mean(groupval_inds_to_plot,:))'); 
%         hplot.LineWidth = 1;
hax = gca;
for ival = 1:nvals_to_plot
    hplot(ival).LineWidth = plotops.linewidth;
end

xlim(xlimits)

ylimdefault = ylim;
if ~isempty(y_ax_hardlims)
    ylim([max(y_ax_hardlims(1),ylimdefault(1)), min(y_ax_hardlims(2),ylimdefault(2))])
end

htitle = title([thissub, '___', channame, '... ', full_var_string], 'Interpreter','none');
%     hleg = legend(resp_grpd.groupval{groupval_inds_to_plot});

f=get(gca,'Children');
hleg = legend(flipud(f(end-nvals_to_plot+1:end)),resp_grpd.groupval{groupval_inds_to_plot});

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

end

if plot_raster
    imagesc(resp_align.resp)
%     xlabel('Time (sec)')
    ylabel('Trial')

% % % % % % % % %     hold on
% % % % % % % % %         % stim syllable onsets
% % % % % % % % %     hstim_syl(1) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,1),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style);
% % % % % % % % %     hstim_syl(2) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style);
% % % % % % % % %     hstim_syl(3) = xline(xline_fn(trials_tmp.stim_syl_on_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_on, 'LineStyle',xline_style);
% % % % % % % % %     
% % % % % % % % %     % stim offset
% % % % % % % % %     hstim_syl(3) = xline(mean(trials_tmp.stim_syl_off_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_stim_syl_off, 'LineStyle',xline_style);
% % % % % % % % %     
% % % % % % % % %     % produced syllable onsets
% % % % % % % % %     hprod_syl(1) = xline(0, 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style); 
% % % % % % % % %     hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_on_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style);
% % % % % % % % %     hprod_syl(3) = xline(xline_fn(trials_tmp.prod_syl_on_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_on, 'LineStyle',xline_style);
% % % % % % % % %     
% % % % % % % % %     % produced syllable offsets
% % % % % % % % %     hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,1),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style);
% % % % % % % % %     hprod_syl(2) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,2),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style);
% % % % % % % % %     hprod_syl(3) = xline(xline_fn(trials_tmp.prod_syl_off_adj(:,3),'omitnan'), 'LineWidth',xline_width, 'Color',xline_color_prod_syl_off, 'LineStyle',xline_style);

end