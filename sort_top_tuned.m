 %%%% find and plot the electrodes which are best tuned for a given parameter
 .... might need to run plot_top_electrodes_mni_on_ctx.m first


% param = 'p_prod_cons_best_anypos';
% param = 'p_prod_vow_best_anypos';
% param = 'p_prod_syl_best_anypos';
% param = 'p_rank';
% param = 'p_prep';
% param = 'p_prep_syl_best_anypos';
% param = {'p_prep_cons',1};
% param = {'p_prep_cons',2};
% param = {'p_prep_cons',3};
% param = {'p_prep_vow',1};
% param = {'p_prep_vow',2};
% param = {'p_prep_vow',3};
% param = {'p_prep_syl',1};
% param = {'p_prep_syl',2};
% param = {'p_prep_syl',3};

% param = {'p_prod_cons',1};
% param = {'p_prod_cons',2};
% param = {'p_prod_cons',3};
% param = {'p_prod_syl',1};
% param = {'p_prod_syl',2};
% param = {'p_prod_syl',3};

% param = 'p_prep_cons_constit';
% param = 'p_prep_vow_constit'; 
param = 'p_prep_syl_constit';

exclude_if_p_zero = 1; % delete channels if they have p=0 for the key parameter

%%
%  load('resp_all_subjects.mat'); 
resp = movevars(resp,{'sub','chan','fs_anatomy','MOREL_label_1'},'Before',1);

[srtvals, varname] = triplet_tablevar(resp, param);
[srtvals, idxorder] = sort(srtvals);

srt = resp(idxorder,:); 
srt = movevars(srt,{varname},'After','MOREL_label_1');

if exclude_if_p_zero
    pzero_rows = srtvals == 0; 
    srt = srt(~pzero_rows,:);    
        clear srtvals idxorder
end