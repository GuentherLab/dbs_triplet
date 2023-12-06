%%%% create bar graphs showing the proportions of electrodes that are significantly tuned for different parameters

ylimit = [0 0.2]; 
ops.barcolor = [0.5 0.5 0.5];
ops.fontsize = 15;
ops.fontname = 'Arial';

resp.sgn_prod_syl_position = resp.p_prod_syl_position < 0.05;
resp.sgn_stim_syl_position = resp.p_stim_syl_position < 0.05;
resp.sgn_prod_cons_position = resp.p_prod_cons_position < 0.05;
resp.sgn_prod_vow_position = resp. p_prod_vow_position < 0.05;
resp.sgn_stim_cons_position = resp.p_stim_cons_position < 0.05;
resp.sgn_stim_vow_position = resp.p_stim_vow_position < 0.05;
resp.sgn_prep = resp.p_prep < 0.05;
resp.sgn_rank = resp.p_rank < 0.05;
resp.sgn_prep_syl1 = resp.p_prep_syl1 < 0.05;
resp.sgn_prep_syl2 = resp.p_prep_syl2 < 0.05;
resp.sgn_prep_syl3 = resp.p_prep_syl3 < 0.05;

close all

%%
hfig = figure; 
hbar = bar(mean(resp.sgn_prod_cons_position));
hax = gca; 
hax.XTickLabels = {'Cons 1', 'Cons 2', 'Cons 3'};
% title('production')
ylabel({'Proportion significantly', 'encoding electrodes'})
hbar.FaceColor = ops.barcolor;
hold on
yline(0.05)
% ylim(ylimit)
box off
set(gca,'FontName',ops.fontname)
set(gca,'FontSize',ops.fontsize)
set(gcf, 'Color', [1 1 1])

hfig = figure; 
hbar = bar(mean(resp.sgn_prod_vow_position));
hax = gca; 
hax.XTickLabels = {'vow 1', 'vow 2', 'vow 3'};
title('production')
hold on
yline(0.05)
ylim(ylimit)

hfig = figure; 
hbar = bar(mean(resp.sgn_prod_syl_position));

hax = gca; 
% hax.XTickLabels = {'Syllable 1', 'Syllable 2', 'Syllable 3'};
hax.XTickLabels = {'Syl 1', 'Syl 2', 'Syl 3'};
set(gca,'FontSize',ops.fontsize)
% title('production')
ylabel({'Proportion significantly', 'encoding electrodes'})
hbar.FaceColor = ops.barcolor;
hold on
yline(0.05)
% ylim(ylimit)
box off
set(gca,'FontName',ops.fontname)
set(gca,'FontSize',ops.fontsize)
set(gcf, 'Color', [1 1 1])


hfig = figure; 
hbar = bar(mean(resp.sgn_stim_cons_position));
hax = gca; 
hax.XTickLabels = {'cons 1', 'cons 2', 'cons 3'};
title('stim')
hold on
yline(0.05)
ylim(ylimit)

hfig = figure; 
hbar = bar(mean(resp.sgn_stim_vow_position));
hax = gca; 
hax.XTickLabels = {'vow 1', 'vow 2', 'vow 3'};
title('stim')
hold on
yline(0.05)
ylim(ylimit)

%%
close all

% ylimit = [0 0.7]; 
hfig = figure; 
hbar = bar([mean(resp.sgn_prep), mean(resp.sgn_rank)]);
hax = gca; 
hax.XTickLabels = {'prep', 'prod rank'};
% title('stim')
ylabel({'Proportion significantly', 'encoding electrodes'})
hold on
yline(0.05)
% ylim(ylimit)
box off
set(gca,'FontName',ops.fontname)
set(gcf, 'Color', [1 1 1])



%%
ylimit = [0 0.25]; 

inds = any(resp.sgn_prod_cons_position,2);

hfig = figure; 
hbar = bar(mean(resp.sgn_stim_cons_position(inds,:)));
hax = gca; 
hax.XTickLabels = {'stim cons 1', 'stim cons 2', 'stim cons 3'};
title('prod-cons encoding electrodes only')
hold on
yline(0.05)
ylim(ylimit)

%%
% ylimit = [0 0.25]; 

inds = resp.sgn_rank;

hfig = figure; 
hbar = bar(mean(resp.sgn_prod_cons_position(inds,:)));
hax = gca; 
hax.XTickLabels = {'prod cons 1', 'prod cons 2', 'prod cons 3'};
% title('rank-order encoding electrodes only')
hold on
yline(0.05)
% ylim(ylimit)

%%
ylimit = [0 0.25]; 
inds = resp.sgn_rank;

hfig = figure; 
hbar = bar(mean([resp.sgn_prep_syl1, resp.sgn_prep_syl2, resp.sgn_prep_syl3]));
hax = gca; 
hax.XTickLabels = {'prod syl 1', 'prod syl 2', 'prod syl 3'};
% title('preparatory activity')
ylabel({'Proportion significantly', 'encoding electrodes'})
hbar.FaceColor = ops.barcolor;
hold on
yline(0.05)
% ylim(ylimit)
box off
set(gca,'FontName',ops.fontname)
set(gca,'FontSize',ops.fontsize)
set(gcf, 'Color', [1 1 1])
