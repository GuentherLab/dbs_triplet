 % remove rows from a trialtable that do not have a corresponding trial in a fieldtrip object
 %
 %  [trialtable_out, trials_ft]  = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D)
%
% inputs: 
%%% cfg.trials = trialtable with rows to be removed
%%% cfg.plot_times = toggle whether figure is created indicated input trial times, fieldtrip times, output trial times
%%% D = fieldtrip_object
%
% outputs: 
%%% trialtable_out = original trialtable with non-fieldtrip rows removed
%%% trials_ft = trialtable created from the fieltrip object timing
%
%   created to cut out of trialtables all rows for which the vibration-denoised trial is missing

function [trialtable_out, trials_ft]  = P08_correct_fieldtrip_trialtable_discrepancies(cfg,D)

% plot colors
cfg.color_ft = [1 0 0]; 
cfg.color_trials_in = [0 0 1]; 
cfg.color_trials_out = [0 0.6 0];

field_default('cfg','plot_times',false)

trials_in = cfg.trials; 
ntrials_in = height(trials_in);
trials_in.ft_idx = nan(ntrials_in,1); 

% create trialtable from 'times' field in fieldtrip object
trials_ft = table(cellfun(@(x)x(1),D.time)', cellfun(@(x)x(end),D.time)', 'VariableNames', {'starts','ends'});

trials_to_keep = false(ntrials_in,1);
for itrial = 1:ntrials_in
    startmatch = trials_in.starts(itrial) >= trials_ft.starts & trials_in.starts(itrial) <= trials_ft.ends;
    endmatch = trials_in.ends(itrial) >= trials_ft.starts & trials_in.ends(itrial) <= trials_ft.ends;
    trialmatch = find(startmatch & endmatch);
    if isempty(trialmatch)
       trials_to_keep(itrial) = false; 
    elseif ~isempty(trialmatch)
       trials_to_keep(itrial) = true; 
       trials_in.ft_idx(itrial) = trialmatch(1); % index of the first trial in the fieldtrip object that contains this trialtable timing
    end
end

trialtable_out = trials_in(trials_to_keep,:); 

n_ft_trials = height(trials_ft); 
trials_ft.session_id = nan(n_ft_trials,1);
trials_ft.session_id(trialtable_out.ft_idx) = trialtable_out.session_id; % copy session id from the epoch trialtable to the fieldtrip trialtable
% for each unclaimed ft trial, assign the either the preceding or next trial's session, depending on which has the closer start time
 for i_ft_trial = 1:n_ft_trials
    if isnan(trials_ft.session_id(i_ft_trial))
       ftstarts = trials_ft.starts;
       ftstarts(i_ft_trial) = inf; % exclude this trial
       ftstarts(isnan(trials_ft.session_id)) = inf; % exlude other trials which do not have session ID assigned
       [~,closest_trial_idx] = min(abs(ftstarts - trials_ft.starts(i_ft_trial)));
        trials_ft.session_id(i_ft_trial) = trials_ft.session_id(closest_trial_idx); 
    end
 end

if cfg.plot_times
    hfig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on

    plot_trialtable(trials_ft, 0, cfg.color_ft)
        han = annotation("textbox",'String','fieldtrip','Position',[.8 .2 .1 .1], 'Color',cfg.color_ft, 'FontWeight','bold');
    plot_trialtable(trials_in, 0.2, cfg.color_trials_in)
        han = annotation("textbox",'String','trials in','Position',[.8 .15 .1 .1], 'Color',cfg.color_trials_in, 'FontWeight','bold');
    plot_trialtable(trialtable_out, 0.4, cfg.color_trials_out) 
        han = annotation("textbox",'String','trials out','Position',[.8 .1 .1 .1], 'Color',cfg.color_trials_out, 'FontWeight','bold');
    xlabel('seconds GTC')

    hold off
end


%% 
function hscat = plot_trialtable(tab, yoffset, color)
    jittermag = 0.0; 
    hold on
    for itrial = 1:height(tab)
        yval = itrial + yoffset + jittermag*rand(1,1); 
        hplot(itrial,1) = plot([tab.starts(itrial), tab.ends(itrial)], [yval yval], 'color', color);
    end
end
%%



end