%%%% plut muliple triplet timecourses in one figure

% close all
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

rowlist = 1:6; 
% rowlist = 7:12;
% rowlist = 13:18; 
% rowlist = 19:24; 
% rowlist = 25:30;
% rowlist = 31:36;
% rowlist = 1143:1148;

nplotrows = 2; 

%%%%% either plot just timecourses, or timecourses plus brains
plot_brains_on_row2 = 0; 

%%% choose the stimulus variable which will be used to sort trials
% sort_cond = []; % do not sort by trial condition; average all trials
% sort_cond = 'stim_volume'; 
% sort_cond = {'cons',1};
sort_cond = {'cons',2};
% sort_cond = {'cons',3};
% sort_cond = {'vow',1};
% sort_cond = {'vow',2};
% sort_cond = {'vow',3};
% sort_cond = {'syl',1}; 
% sort_cond = {'syl',2}; 
% sort_cond = {'syl',3}; 
% sort_cond = 'cons_constit';
% sort_cond = 'vow_constit'; 
% sort_cond = 'syl_constit';

plotop.x_ax_hardlims = [-3 2]; 


%%
nelcs = length(rowlist);

hfig = figure('Color','w');

if ~plot_brains_on_row2
    for ielc = 1:nelcs
        thisrow = rowlist(ielc); 
        subplot(nplotrows,nelcs/nplotrows,ielc);
        srt_row = thisrow;
        newfig = 0;
        plot_resp_timecourse_triplet()
    
        ylimdefault = ylim;
%         ylim([ylimdefault(1), min(y_axmax,ylimdefault(2))])
    
    end
elseif plot_brains_on_row2
    for ielc = 1:nelcs
        thisrow = rowlist(ielc);
        subplot(2,nelcs,ielc)
        srt_row = thisrow;
        newfig = 0;
        plot_resp_timecourse_triplet()
    
        ylimdefault = ylim;
%         ylim([ylimdefault(1), min(y_axmax,ylimdefault(2))])

        subplot(2, nelcs, nelcs+ielc)
        plot_sorted_resp_mni_on_ctx

    end    
end
 