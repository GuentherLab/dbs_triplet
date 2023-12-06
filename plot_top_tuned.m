 %%%% find and plot the electrodes which are best tuned for a given parameter

param = 'p_prep_syl3';

% delete channels if they have nan for the following parameter
exclude_if_nan_param = 'p_prep';

%  load('resp_all_subjects.mat'); 
resp = movevars(resp,{'sub','chan'},'Before',1);

srt = sortrows(resp,param); 


srt = movevars(resp,{param},'After','chan');