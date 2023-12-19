 %%%% find and plot the electrodes which are best tuned for a given parameter
 .... might need to run plot_top_electrodes_mni_on_ctx.m first


% param = 'p_prod_cons_mean';
% param = 'p_prod_vow_mean';
% param = 'p_prod_syl_mean';
% param = 'p_rank';
% param = 'p_prep';
param = 'p_prep_syl_mean';
% param = 'p_prep_syl1';
% param = 'p_prep_syl2';
% param = 'p_prep_syl3';

% % % % % % % % % % % % % % delete channels if they have nan for the following parameter
% % % % % % % % % % % % % exclude_if_nan_param = 'p_prep';

exclude_if_p_zero = 1; % delete channels if they have p=0 for the key parameter

%%
%  load('resp_all_subjects.mat'); 
resp = movevars(resp,{'sub','chan'},'Before',1);

srt = sortrows(resp,param); 
srt = movevars(srt,{param},'After','chan');

if exclude_if_p_zero
    pzero_rows = srt{:,param} == 0; 
    srt = srt(~pzero_rows,:); 
end