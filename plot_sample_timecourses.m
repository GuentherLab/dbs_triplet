%%%% plut muliple triplet timecourses in one figure

% close all

% rowlist = 1:6; 
rowlist = 7:12;
% rowlist = 13:18; 
% rowlist = 19:24; 
% rowlist = 25:30;
% rowlist = 31:36;

% y_axmax = 3.5; 

nplotrows = 2; 

%%%%% either plot just timecourses, or timecourses plus brains
plot_brains_on_row2 = 0; 


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
 