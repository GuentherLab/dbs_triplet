%%%% brainplot the top electrodes for a particular parameter 

% close all
% clear

% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

%% load electrode responses and mni coords

% load([PATH_RESULTS filesep 'resp_all_subjects'])
% 
n_elc = height(resp);

struct_to_plot = 'ctx';
% struct_to_plot = 'stn';
% struct_to_plot = 'thal';

%% Configuration Variables and Paths
switch struct_to_plot
    case 'ctx'
        %loading cortical reconstructions
        average_mni = load(PATH_AVERAGE_MNI);
    case 'stn'

    case 'thal'
        % subcort = load(PATH_SUBCORT_ATLAS);
        % subcort_vim = load(PATH_SUBCORT_ATLAS_VIM);
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








%% set params, make brainplot

inclusion_mode = 'thresh';
% inclusion_mode = 'proportion';

% p_thresh = 0.01; 
p_thresh = 0.00001; 
% p_thresh = 0.05 / 3; % bonf corrected 0.05

p_proportion = 0.01; 

% param = 'p_prod_cons_best_anypos';
% param = 'p_prod_vow_best_anypos';
% param = 'p_prod_syl_best_anypos';
% param = 'p_rank';
% param = 'p_prep';
% param = 'p_prep_syl_best_anypos';
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};
% param = {'p_prep_cons',1};
% param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
% param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; 
param = 'p_prep_syl_constit';

exclude_if_p_zero = 1; % exclude channels if they have p=0 for the key parameter



if exclude_if_p_zero
    excluded_rows = triplet_tablevar(resp, param) == 0; 
elseif ~exclude_if_p_zero
    excluded_rows = false(n_elc,1);
end

switch inclusion_mode
    case 'thresh'
        rows_to_plot = triplet_tablevar(resp, param) < p_thresh & ~excluded_rows;
    case 'proportion'
        varvals = triplet_tablevar(resp, param);
        varvals(excluded_rows) = nan; 
        [~, rows_ranked] = sort(varvals);
        rows_to_plot = rows_ranked( 1:round(p_proportion * n_elc) ); 
end


plotcolor = 'r';


snap_to_surf = 1; % if true, project eletrodes to nearest point on ctx surface

% shift electrodes so that they aren't covered by the brain surface
%%% gets applied after snapping to surface
%%% .... if snapping, offset of -1 should be enough to have points entirely above ctx surface (in L hem)
x_offset = -1;

elc_to_plot = resp(rows_to_plot,{'mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'}); 


figure;
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

xyz_to_plot_nonsnapped = [elc_to_plot.mni_nonlinear_x, elc_to_plot.mni_nonlinear_y, elc_to_plot.mni_nonlinear_z];
if snap_to_surf
    [~, surfpoint_idx] = min(pdist2(xyz_to_plot_nonsnapped,average_mni.Vertices), [], 2); % find nearest surf points
    xyz_to_plot = average_mni.Vertices(surfpoint_idx,:); 
elseif ~snap_to_surf
    xyz_to_plot = xyz_to_plot_nonsnapped;
end

hscat = scatter3(xyz_to_plot(:,1) + x_offset, xyz_to_plot(:,2), xyz_to_plot(:,3), 'filled',...
  'MarkerFaceAlpha',1,'MarkerFaceColor',plotcolor,'MarkerEdgeColor','k','LineWidth',0.01);
hscat.SizeData = 40;
set(gcf, 'Color', [1 1 1]); % white backgroud
view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
% % % % % % % % % % % scalebar(0,70,-50, 10, 'mm')

titlestr = param; 
title(titlestr,'interpreter', 'none')

% print(gcf,[PATH_ANALYSIS 'qqq.png'],'-dpng','-r300')

