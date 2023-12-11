 % remove rows from a trialtable that do not have a corresponding trial in a fieldtrip object
 % inputs: 
%%% cfg.trials = trialtable with rows to be removed
%%% D = fieldtrip_object
%
% outputs: 
%%% trialtable_out = original trialtable with non-fieldtrip rows removed
%%% trials_ft = trialtable created from the fieltrip object timing
%
%   created to cut out of trialtables all rows for which the vibration-denoised trial is missing

function [trialtable_out, trials_ft]  = remove_trialrows_missing_from_fieldtrip(cfg,D)

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

if cfg.plot_times
    hfig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on

    plot_trialtable(trials_ft, 0, [1 0 0])
    plot_trialtable(trials_in, 0.2, [0 0 1])
    plot_trialtable(trialtable_out, 0.4, [0 0.6 0])
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