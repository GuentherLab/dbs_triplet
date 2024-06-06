%%%% brainplot the top electrodes for a particular parameter 

% close all
% clear
clearvars -except resp subs

% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

setpaths_dbs_triplet()

%% set params

 struct_to_plot = 'ctx';
    snap_to_surf = 1; % cortex only - if true, project eletrodes to nearest point on ctx surface
    % shift electrodes so that they aren't covered by the brain surface
    %%% gets applied after snapping to surface
    %%% .... if snapping, offset of -1 should be enough to have points entirely above ctx surface (in L hem)
    x_offset = -1;
%struct_to_plot = 'stn';
% struct_to_plot = 'thal';

%%% pick hemisphere to plot - subcortical only
side = 'L'; 
% side = 'R'; 

inclusion_mode = 'thresh';
    p_thresh = 0.001; 
    % p_thresh = 0.00001; 
    % p_thresh = 0.05 / 3; % bonf corrected 0.05
% inclusion_mode = 'proportion';
    p_proportion = 0.01; 

% param = 'p_rank';
% param = {'p_prod_syl',1};
% param = {'p_prod_syl',2};
% param = {'p_prod_syl',3};
% param = {'p_prod_cons',1};
%  param = {'p_prod_cons',2};
% param = {'p_prod_cons',3};
% param = {'p_prod_vow',1};
% param = {'p_prod_vow',2};
% param = {'p_prod_vow',3};
% % % % % % % % param = 'p_prod_cons_best_anypos';
% % % % % % % % param = 'p_prod_vow_best_anypos';
% % % % % % % param = 'p_prod_syl_best_anypos';
% param = 'p_prep';
% param = 'p_prep_syl_best_anypos';
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};
param = {'p_prep_cons',1};
% param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
% param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; 
% param = 'p_prep_syl_constit';
% param = {'p_trans_id',1}
exclude_if_p_zero = 1; % exclude channels if they have p=0 for the key parameter

also_plot_nonsgnf_elcs = 1; % if true, plot non-significant electrodes alongside significant electrodes

plotcolor = 'r';
plotcolor_nonsgn = [0.6 0.6 0.6]; 

marker_size = 40; % size of electrode marker; scatter 'SizeData' parameter
marker_size_nonsgn = 5; 

view_angle = [-90, 0]; % use [-90, 0] for straight-on lateral left hemisphere



%% load electrode responses and mni coords

% load([PATH_RESULTS filesep 'resp_all_subjects'])
% 
n_elc = height(resp);



%% load brain surfaces
switch struct_to_plot
    case 'ctx'
        %loading cortical reconstructions
        average_mni = load(PATH_AVERAGE_MNI);
    case 'stn'
        subcort_stn = load(PATH_STN_ATLAS);
    case 'thal'
        subcort_vim = load(PATH_SUBCORT_ATLAS_VIM);
        % nii_vimi = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/91.nii.gz');
        % nii_vime = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/94.nii.gz');
        % nii_vimip = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/104.nii.gz');
        % nii_vimep = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/122.nii.gz');
        
        % 
        % subcort_vimi_lh_fv = ea_nii2fv(nii_vimi);
        % subcort_vime_lh_fv = ea_nii2fv(nii_vime);
        % subcort_vimip_lh_fv = ea_nii2fv(nii_vimip);
        % subcort_vimep_lh_fv = ea_nii2fv(nii_vimep);
        
        % % % %loading VL posterior ventral from Morel atlas
        % nii_vlpv = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/MorelAtlasICBM2009b (Jakab 2008)/lh/VLpv.nii.gz');
        % subcort_vlpv_lh_fv = ea_nii2fv(nii_vlpv);
end


%% make brainplot

switch side
    case 'L'
        side_number = 2;
    case 'R'
        side_number = 1; 
end

if exclude_if_p_zero
    excluded_rows = triplet_tablevar(resp, param) == 0; 
elseif ~exclude_if_p_zero
    excluded_rows = false(n_elc,1);
end

switch inclusion_mode
    case 'thresh'
        sgn_rows = triplet_tablevar(resp, param) < p_thresh & ~excluded_rows;
    case 'proportion'
        varvals = triplet_tablevar(resp, param);
        varvals(excluded_rows) = nan; 
        [~, rows_ranked] = sort(varvals);
        sgn_rows = rows_ranked( 1:round(p_proportion * n_elc) ); 
end

hfig = figure;

switch struct_to_plot
    case 'ctx'
        rows_to_plot_sgn = sgn_rows & string(resp.type)=="ecog"; 
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 
        xyz_to_plot_nonsnapped = [elc_to_plot.mni_nonlinear_x, elc_to_plot.mni_nonlinear_y, elc_to_plot.mni_nonlinear_z];

        rows_to_plot_nonsgn = ~sgn_rows & string(resp.type)=="ecog"; 
        elc_to_plot_nonsgn = resp(rows_to_plot_nonsgn,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 
        xyz_to_plot_nonsnapped_nonsgn = [elc_to_plot_nonsgn.mni_nonlinear_x, elc_to_plot_nonsgn.mni_nonlinear_y, elc_to_plot_nonsgn.mni_nonlinear_z];

        if snap_to_surf
            [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped,average_mni.Vertices), [], 2); % find nearest surf points
            xyz_to_plot = average_mni.Vertices(surfpoint_idx,:); 
            [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped_nonsgn,average_mni.Vertices), [], 2); % find nearest surf points
            xyz_to_plot_nonsgn = average_mni.Vertices(surfpoint_idx,:); 
        elseif ~snap_to_surf
            xyz_to_plot = xyz_to_plot_nonsnapped;
            xyz_to_plot_nonsgn = xyz_to_plot_nonsnapped_nonsgn;
        end

        hpatch = patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
            'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
            'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);


    case 'stn'
        rows_to_plot_sgn = sgn_rows & contains(resp.type,{'dbs';'macro'}) & contains(resp.MOREL_label_1,{['STh_',side]});
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 
        xyz_to_plot = [elc_to_plot.mni_nonlinear_x, elc_to_plot.mni_nonlinear_y, elc_to_plot.mni_nonlinear_z];

        rows_to_plot_nonsgn = ~sgn_rows & contains(resp.type,{'dbs';'macro'}) & contains(resp.MOREL_label_1,{['STh_',side]});
        elc_to_plot_nonsgn = resp(rows_to_plot_nonsgn,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 
        xyz_to_plot_nonsgn = [elc_to_plot_nonsgn.mni_nonlinear_x, elc_to_plot_nonsgn.mni_nonlinear_y, elc_to_plot_nonsgn.mni_nonlinear_z];

        hpatch = patch('vertices', subcort_stn.atlases.fv{1,side_number}.vertices, 'faces', subcort_stn.atlases.fv{1,side_number}.faces,...
            'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
            'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);

    case 'thal' %%%%%% need to find an appropriate VIM atlas before using this option
        rows_to_plot_sgn = sgn_rows & contains(resp.type,{'dbs';'macro'});
        rows_to_plot_sgn = rows_to_plot_sgn & contains(resp.MOREL_label_1,{['VApc_',side],['VLa_',side],['VLpv_',side],['VM_',side],['VPM_',side]});
        elc_to_plot = resp(rows_to_plot_sgn,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 
        xyz_to_plot = [elc_to_plot.mni_nonlinear_x, elc_to_plot.mni_nonlinear_y, elc_to_plot.mni_nonlinear_z];

        % hpatch = patch('vertices', subcort_vim.atlases.fv{???????,side_number}.vertices, 'faces', subcort_vim.atlases.fv{???????,side_number}.faces,...
        %     'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
        %     'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);   

end

hold on

if  also_plot_nonsgnf_elcs
    hscat_non_sgnf = scatter3(xyz_to_plot_nonsgn(:,1) + x_offset, xyz_to_plot_nonsgn(:,2), xyz_to_plot_nonsgn(:,3), 'filled',...
       'MarkerFaceAlpha',1,'MarkerFaceColor',plotcolor_nonsgn,'MarkerEdgeColor','k','LineWidth',0.01);
    hscat_non_sgnf.SizeData = marker_size_nonsgn;
end


hscat_sgnf = scatter3(xyz_to_plot(:,1) + x_offset, xyz_to_plot(:,2), xyz_to_plot(:,3), 'filled',...
   'MarkerFaceAlpha',1,'MarkerFaceColor',plotcolor,'MarkerEdgeColor','k','LineWidth',0.01);
hscat_sgnf.SizeData = marker_size;



set(gcf, 'Color', [1 1 1]); % white backgroud
view(view_angle(1),view_angle(2))
axis off; axis equal
camlight('headlight','infinite');
% % % % % % % % % % % scalebar(0,70,-50, 10, 'mm')

titlestr = param; 
title(titlestr,'interpreter', 'none')

% print(gcf,[PATH_ANALYSIS 'qqq.png'],'-dpng','-r300')

% At the end of your plotting code, after you've plotted both the significant
% and non-significant electrodes (if applicable), add the legend.
% Note: Adjust the legend entries based on your specific plot and preferences.

legendEntries = {}; % Initialize a cell array to hold legend entries.
legendHandles = []; % Initialize an array to hold legend handles (the plotted objects).

% Check if non-significant electrodes are plotted and add them to the legend.
if also_plot_nonsgnf_elcs
    legendHandles(end+1) = hscat_non_sgnf; % Add the handle for non-significant electrodes.
    legendEntries{end+1} = 'Non-significant electrodes'; % Add the legend entry.
end

% Add the significant electrodes to the legend.
legendHandles(end+1) = hscat_sgnf; % Add the handle for significant electrodes.
legendEntries{end+1} = ['Significant electrodes (alpha = ' num2str(p_thresh) ')']; % Add the legend entry with p_thresh.

% Create the legend.
legend(legendHandles, legendEntries, 'Location', 'bestoutside');

% Optionally, you can adjust the font size of the legend.
% set(legend,'FontSize',8);

