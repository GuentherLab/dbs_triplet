 %%% compare tuning strengths across areas for multiple parameters

%  clear
%  load('Z:\DBS\Analysis\triplet_results_am\resp_all_subjects.mat')

[paramlist] = expandarray('params_to_plot_20231217.csv');

ylimits = [0 0.25]; 

close all force

%%

hfig = figure;
set(0, 'DefaultTextInterpreter', 'none')
set(0,'DefaultFigureWindowStyle','docked')
nparams = numel(paramlist);
[figheight, figwidth] = size(paramlist); 
pt = paramlist';  % select from transposed params so that it match subplot ordering
for iparam = 1:nparams
    subplot(figheight, figwidth, iparam); 
    param = pt{iparam};
    compare_areal_tuning_20231211()

    hbar = bar(areastats.prop_sgn);
    hax = gca;
    hax.XTickLabels = areastats.region;
    hyline = yline(pthresh);
    titlestr = {[param{1,1} '_' num2str(param{1,2})], ['p = ' num2str(chi_p)] };
    title(titlestr); 
    ylabel({'proportion', 'sgnf. electrodes'})
    if exist('ylimits','var')
        ylim(ylimits)
    end

    ylabel('')

%     annostr = ['p = ' num2str(chi_p)];
%     hannot = annotation("textbox",'String',annostr,'Position',[.8 .1 .1 .1], 'Color', [0 0 0], 'FontWeight','bold');
end

alpha_bonferroni = 0.05 / nparams







%% function for importing and expanding params cell array
%%%%% expand each param for ranks 1,2, and 3 concatenating horizontally
%%%%% .... input is the filename of a table file containing param names
function [expandedCellArray] = expandarray(paramfile_name)

cellimport = table2cell(readtable(paramfile_name,'ReadVariableNames',0,'Delimiter','\t'));
[pheight, pwidth] = size(cellimport);

cel123 = cellfun(@(x){{x,1},{x,2},{x,3}},cellimport,'UniformOutput',false);

% Initialize the expanded cell array
expandedCellArray = cell(pheight, 3*pwidth);

% Loop through each cell in the original array
for i = 1:size(cel123, 1)
    for j = 1:size(cel123, 2)
        % Extract the 1x3 cell array from the current cell
        currentCell = cel123{i, j};
        
        % Expand and concatenate the elements horizontally
        for k = 1:3
            expandedCellArray{i, (j-1)*3 + k} = currentCell{k};
        end
    end
end

end