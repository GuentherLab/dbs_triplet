rowlist = 1:6; 
% rowlist = 18:21; 
% rowlist = 22:25;

y_axmax = 3.5; 

nplotrows = 2; 

% either plot just timecourses, or timecourses plus brains
plot_brains_on_row2 = 0; 

nelcs = length(rowlist);

if ~plot_brains_on_row2
    for ielc = 1:nelcs
        thisrow = rowlist(ielc); 
        subplot(nplotrows,nelcs/nplotrows,ielc);
        srt_row = thisrow;
        plot_resp_timecourse
    
        ylimdefault = ylim;
        ylim([ylimdefault(1), min(y_axmax,ylimdefault(2))])
    
    end
elseif plot_brains_on_row2
    for ielc = 1:nelcs
        thisrow = rowlist(ielc);
        subplot(2,nelcs,ielc)
        srt_row = thisrow;
        plot_resp_timecourse
    
        ylimdefault = ylim;
        ylim([ylimdefault(1), min(y_axmax,ylimdefault(2))])

        subplot(2, nelcs, nelcs+ielc)
        plot_sorted_resp_mni_on_ctx

    end    
end
 