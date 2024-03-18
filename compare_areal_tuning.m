 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas
  % load resp_all_subjects first

% close all

 %% params
vardefault('show_barplot',1);

newfig = 1; 

%%% define anatomical regions composed of smaller areas
regiondef = {   'mfg',  {'rostralmiddlefrontal' , 'caudalmiddlefrontal'};... middle frontal gyrus... maybe also inf front sulcus
                'ifg',  {'parsopercularis', 'parsorbitalis',  'parstriangularis'};... % inferior frontal gyrus
                'smc',  {'postcentral', 'precentral'};...                   % sensorimotor cortex
                'suptemp', {'superiortemporal', 'bankssts' }; ... % superior temporal
                % 'thal', {'Left-Thalamus-Proper' , 'Right-Thalamus-Proper' };... % thalamus.... use MOREL labels instead
                % 'wm',   {'Left-Cerebral-White-Matter', 'Right-Cerebral-White-Matter'}; ... % white matter.... use MOREL labels instead
                % 'ventdc', {'Left-VentralDC', 'Right-VentralDC'}... %% ? not sure what VentralDC is.... use MOREL labels instead
                 % 'sfg'   {'superiorfrontal'};... % superior frontal gyrus..... should be excluded - only has 2 electrodes
                 'stn', {'STh_L', 'STh_L'};...
                 'thal', {'VApc_L','VLa_L' , 'VLpv_L', 'VM_L' 'VM_R', 'VPM_L'};...
%                  'gp', {'GPe_L','GPi_sensorimotor_L'}; % <20 electrodes, so maybe not worth including
                };

param = 'p_rank'; 
% param = 'p_prep';

% param = {'p_stim_cons',1};
% param = {'p_stim_cons',2};
% param = {'p_stim_cons',3};
% param = {'p_stim_syl',1};
% param = {'p_stim_syl',2};
% param = {'p_stim_syl',3};
% param = {'p_stim_vow',1};
% param = {'p_stim_vow',2};
% param = {'p_stim_vow',3};

% param = {'p_prep_cons',1};
% param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};

% param = 'p_stim_cons_allpos';
% param = 'p_stim_vow_allpos';
% param = 'p_stim_syl_allpos';

% param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; 
% param = 'p_prep_syl_constit';

% param = 'p_prod_cons_allpos';
% param = 'p_prod_vow_allpos';
% param = 'p_prod_syl_allpos';


pthresh = 0.05; 

bar_face_color = [0.5 0.5 0.5]; 

%% analysis
[paramvals, param_name, full_param_string] = triplet_tablevar(resp, param); 
paramvalid = ~isnan(paramvals) & paramvals ~= 0; % electrodes with usable p values
paramsgn = paramvals < pthresh & paramvalid; % analyzable electrodes significantly tuned for param of interest

nelc = height(resp);
resp.region = cell(nelc,1); 

nregions = size(regiondef,1);
areastats = table(regiondef(:,1), regiondef(:,2), nan(nregions,2), 'VariableNames', {'region','subareas','ebar_lims'});

for iregion = 1:nregions
    thisregion = areastats.region{iregion};
    regionmatch1 = rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp,'InputVariables','fs_anatomy');
        regionmatch1 = any(regionmatch1{:,1},2);
    regionmatch2 = rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp,'InputVariables','MOREL_label_1');
        regionmatch2 = any(regionmatch2{:,1},2);
    regionmatch3 = rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp,'InputVariables','DISTAL_label_1');
        regionmatch3 = any(regionmatch3{:,1},2);
    regionmatch = regionmatch1 | regionmatch2 | regionmatch3; 
    resp.region(regionmatch) = repmat({thisregion},nnz(regionmatch),1);    areastats.nelc(iregion) = nnz(regionmatch); % total electrodes in this region
    areastats.nelc_valid(iregion) = nnz(regionmatch & paramvalid); % number of analyzable electrodes in this region for the param of interest 
    areastats.nelc_sgn(iregion) = nnz(regionmatch & paramsgn); % number of analyzable electrodes significantly tuned for param of interest in this region
    areastats.prop_sgn(iregion) = areastats.nelc_sgn(iregion) / areastats.nelc_valid(iregion); % proportion of tuned electrodes in this region

    %%%% compute error bar values - 95% confidence intervals using binomial test on each area independently
    % move fieldtrip version of binocdf to bottom of path so that we use matlab inbuilt version
    oldpath = path; path(oldpath, [PATH_FIELDTRIP_CODE filesep '\external\stats']); clear oldpath 
    alpha=.0001:.0001:.9999; 
    p = binocdf(areastats.nelc_sgn(iregion), areastats.nelc_valid(iregion), alpha);
    areastats.ebar_lims(iregion,1:2) = alpha([find(p>.975,1,'last'),find(p<.025,1,'first')]);
end


% get the overall proportion of tuned electrodes out of this with a region assigned and analyzable p values
%%%%% alternative computation: region_assigned = ~cellfun(@isempty,resp.region); proportion_signficant_overall = mean( paramsgn(paramvalid & region_assigned) ); 
proportion_signficant_overall = sum(areastats.nelc_sgn) / sum(areastats.nelc_valid); 

% number of tuned electrodes per region if they were randomly distributed across regions
expected_sgn_per_region_random = proportion_signficant_overall * areastats.nelc_valid; 

[chi_significant, chi_p, chi_stats] = chi2gof([1:nregions]', 'Frequency',areastats.nelc_sgn, 'Expected',expected_sgn_per_region_random, 'Emin',0);
% chi_p



%% plotting
if show_barplot

    if newfig 
        hfig = figure('color','w');
    end
    hbar = bar(areastats.prop_sgn);

    hold on

    ebar_neg =  areastats.prop_sgn - areastats.ebar_lims(:,1); 
    ebar_pos =  -areastats.prop_sgn + areastats.ebar_lims(:,2); 
    h_ebar = errorbar([1:nregions]', areastats.prop_sgn, ebar_neg, ebar_pos,'--');
    h_ebar.LineWidth = 0.8;
    h_ebar.LineStyle = 'none';
    h_ebar.Color = [0 0 0];

    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(pthresh);
    set(0, 'DefaultTextInterpreter', 'none')
    titlestr = [full_param_string, '..... p = ' num2str(chi_p)] ;
    htitle = title(titlestr); 

    hbar.FaceColor = bar_face_color; 

end

hold off


