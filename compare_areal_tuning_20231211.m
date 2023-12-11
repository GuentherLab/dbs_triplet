 %%%% check whether there is a nonrandom distribution of significantly tuned electrodes across areas



regiondef = {   'sfg'   {'superiorfrontal'};... % superior frontal gyrus
                'mfg',  {'rostralmiddlefrontal' , 'caudalmiddlefrontal'};... middle frontal gyrus... maybe also inf front sulcus
                'ifg',  {'parsopercularis', 'parsorbitalis',  'parstriangularis'};... % inferior frontal gyrus
                'smc',  {'postcentral', 'precentral'};...                   % sensorimotor cortex
                'suptemp', {'superiortemporal', 'bankssts' }; ... % superior temporal
                'thal', {'Left-Thalamus-Proper' , 'Right-Thalamus-Proper' };... % thalamus
                'wm', {'Left-Cerebral-White-Matter'. 'Right-Cerebral-White-Matter'}... % white matter

% [chi_significant, chi_p, chi_stats] = chi2gof([1:nclusts]', 'Frequency',left_n_per_clust, 'Expected',proportion_left_overall*n_per_clust, 'Emin',0)