 %%% compare tuning strengths across areas for multiple parameters


paramlist = { 'p_prod_cons_best_anypos', 'Produced consonant';...
                'p_prep_syl_best_anypos', 'Preparatory (any syllable)';...
                'p_prod_vow_best_anypos',  'Produced vowel';... 
                'p_prep_syl1',              'Preparatory (syllable 1)';... 
                'p_prod_syl_best_anypos',  'Produced syllable'; ...
                'p_prep_syl2',              'Preparatory (syllable 2)';... 
                'p_rank',                 'Rank Produced';... 
%                 'p_prep';...
                'p_prep_syl3',              'Preparatory (syllable 3)';... 
};

ylimits = [0 0.35]; 

hfig = figure;
% htile = tiledlayout(4,2)
% subplot(4,2)

nparams = size(paramlist,1);
for iparam = 1:nparams
    param = paramlist{iparam,1};
    compare_areal_tuning_20231211()

    subplot(4,2,iparam)
    hbar = bar(areastats.prop_sgn);
    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(pthresh);
    titlestr = [paramlist{iparam,2}, '..... p = ' num2str(chi_p)];
    title(titlestr); 
    ylabel({'proportion', 'sgnf. electrodes'})
    ylim(ylimits)
%     annostr = ['p = ' num2str(chi_p)];
%     hannot = annotation("textbox",'String',annostr,'Position',[.8 .1 .1 .1], 'Color', [0 0 0], 'FontWeight','bold');
%     nexttile
end