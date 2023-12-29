 %%% compare tuning strengths across areas for multiple parameters

clearvars -except resp subs srt

% paramlist = { {'p_rank',1},                 'Rank Produced';... 
%                 {'p_prep_syl',1},              'Preparatory (syllable 1)';... 
%                 {'p_prod_syl',1},              'Production (syllable 1)';... 
%                 {'p_prod_syl',2},              'Production (syllable 2)';... 
%                 {'p_prod_syl',3},              'Production (syllable 3)';... 
% };
% subinds_to_plot = [1 2 4 5 6]; 

paramlist = { {'p_rank',1},                 'Rank Produced';... 
                {'p_prod_syl',1},              'Production (syllable 1)';... 
                {'p_prep_syl',1},              'Preparatory (syllable 1)';... 
                {'p_prod_syl',2},              'Production (syllable 2)';... 
                {'p_prod_syl',3},              'Production (syllable 3)';... 
};
subinds_to_plot = [1 2 3 4 6]; 

% ylimits = [0 0.25]; 

bar_face_color = 0.7 * ones(1,3); 

hfig = figure('Color',[1 1 1]);



nparams = size(paramlist,1);
for iparam = 1:nparams
    param = paramlist{iparam,1};
    compare_areal_tuning_20231218_ISSP_poster()

    subind = subinds_to_plot(iparam); 
    subplot(3,2,subind)
    hbar = bar(areastats.prop_sgn, 'FaceColor', bar_face_color);
    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(pthresh);
    titlestr = [paramlist{iparam,2}, '..... p = ' num2str(chi_p)];
%     titlestr = [paramlist{iparam,2}];
    title(titlestr); 
    ylabel({'proportion', 'sgnf. electrodes'})
    if exist('ylimits','var')
        ylim(ylimits)
    end
    box off
%     annostr = ['p = ' num2str(chi_p)];
%     hannot = annotation("textbox",'String',annostr,'Position',[.8 .1 .1 .1], 'Color', [0 0 0], 'FontWeight','bold');
end

%% brainplot
subplot(3,2,5)
hold off
plot_top_electrodes_mni_on_ctx_20231218_ISSP_poster()
title({'Electrodes selective for','production (syllable 1)'})


