 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas
  % load resp_all_subjects first

 %% params
vardefault('show_barplot',1);

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
                };

% param = 'p_prod_cons_mean';
% param = 'p_prod_vow_mean'; ;
% param = 'p_prod_syl_mean';;
% param = 'p_rank'; ;
% param = 'p_prep';
% param = 'p_prep_syl_mean';
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};
% param = {'p_prep_cons',1};
% param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; ;
% param = 'p_prep_syl_constit';;

pthresh = 0.05; 

%% analysis
paramvals = triplet_tablevar(resp, param); 
paramvalid = ~isnan(paramvals) & paramvals ~= 0; % electrodes with usable p values
paramsgn = paramvals < pthresh & paramvalid; % analyzable electrodes significantly tuned for param of interest

nelc = height(resp);
resp.region = cell(nelc,1); 

areastats = table(regiondef(:,1), regiondef(:,2), 'VariableNames', {'region','subareas'});

nregions = height(areastats);
for iregion = 1:nregions
    thisregion = areastats.region{iregion};
    regionmatch1 = rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp,'InputVariables','fs_anatomy');
        regionmatch1 = any(regionmatch1{:,1},2);
    regionmatch2 = rowfun(@(x)strcmp(x,areastats.subareas{iregion}),resp,'InputVariables','MOREL_label_1');
        regionmatch2 = any(regionmatch2{:,1},2);
    regionmatch = regionmatch1 | regionmatch2; 
    resp.region(regionmatch) = repmat({thisregion},nnz(regionmatch),1);    areastats.nelc(iregion) = nnz(regionmatch); % total electrodes in this region
    areastats.nelc_valid(iregion) = nnz(regionmatch & paramvalid); % number of analyzable electrodes in this region for the param of interest 
    areastats.nelc_sgn(iregion) = nnz(regionmatch & paramsgn); % number of analyzable electrodes significantly tuned for param of interest in this region
    areastats.prop_sgn(iregion) = areastats.nelc_sgn(iregion) / areastats.nelc_valid(iregion); % proportion of tuned electrodes in this region
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

    hfig = figure('color','w');
    hbar = bar(areastats.prop_sgn);
    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(pthresh);
    titlestr = {[param, '..... p = ' num2str(chi_p)] };
    title(titlestr); 

end




