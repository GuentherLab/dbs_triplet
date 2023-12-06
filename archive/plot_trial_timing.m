%%%% plot the onsets and offsets of trials in vibration-denoised data....
% ............. relative to actual stim onset and voice offset
%
%%%% AM 2022/7/23


ft_defaults
bml_defaults

SUBJECT='DBS3012';
PATH_SUBJECT=[PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT filesep 'Preprocessed Data' filesep 'Sync'];
PATH_PROTOCOL = 'C:\Users\amsmeier\Documents\MATLAB\P09_artifact_criteria_E';
PATH_ANNOT = [PATH_SYNC '/annot']; 
trials_to_plot = 1:10;

plot_nonoverlap_trialtimes = 1; % plot redefined trials without overlap

if ~exist('D_hg_trial','var')
    load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_hg_trial_denoised.mat'],'D_hg_trial');
    tt = D_hg_trial.time;
    tr_stim = bml_annot_read([PATH_ANNOT '/' SUBJECT '_stimulus_triplet.txt']);
    tr_prod = bml_annot_read([PATH_ANNOT '/' SUBJECT '_produced_triplet.txt']);
end

trialtimes_no_overlap = table;

close all
hfig = figure('units','normalized','outerposition',[0 0 1 1]);
for itrial = trials_to_plot
    hold on
    col(itrial,:) = rand(1,3); 
    
    % vibration-denoised ephys data trial timing
    hscat(itrial,1) = scatter(tt{itrial}([1,end]),[itrial, itrial],'filled');
    hscat(itrial,1).CData = col(itrial,:);
    
    % stim onset and voice offset trial timing
    hscat(itrial,2) = scatter([tr_stim.starts(itrial), tr_prod.ends(itrial)], [itrial, itrial], 'x');
    hscat(itrial,2).CData = col(itrial,:);
    
    %%%%%% get new trial times that will not overlap
    % trial onsets stay the same
    trialtimes_no_overlap.starts(itrial) = tt{itrial}(1); 
    % for trial offsets, take the latest timepoint of this trial which is....
    % ... before the first timepoint of the subsequent trial
    nonoverlap_timepoints = tt{itrial}(tt{itrial} < tt{itrial+1}(1)); 
    trialtimes_no_overlap.ends(itrial) = max(nonoverlap_timepoints); 
    
    if plot_nonoverlap_trialtimes
        offs = 0.4; 
        hscat(itrial,3) = scatter([trialtimes_no_overlap.starts(itrial), trialtimes_no_overlap.ends(itrial)],...
            [itrial+offs, itrial+offs], 'v');
        hscat(itrial,3).CData = col(itrial,:);
    end
    
end

ylabel('trial')
xlabel('time (sec)')
title([SUBJECT])