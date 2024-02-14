 %%% get specific subcolumn of a table variable
 % use this function to get e.g. the first/second/third consonant column from a trialtable or electrodes response table
%
%   [vals_out, varname, full_var_string]  = triplet_tablevar(tab, varname_and_subcol, rows)
%       vals_out = all values in the specified column of the specified variable
%       varname = string label of table variable [no column info included]
%       full_var_string = string containing table variable name plus column number... for labeling figures
%
%   tab = table input
%   varname_and_subcol = string or 1 or 2 value cell
%           first arg = table variable name (string), 2nd arg = subcolumn index(es) (number)
%           if 2nd value is not supplied, all subcolumns will be outputted
%   rows = table rows to output - defaults to all rows of the table

function [vals_out, varname, full_var_string] = triplet_tablevar(tab, varname_and_subcol, rows)

if iscell(varname_and_subcol)
    varname = varname_and_subcol{1};
    subcolidx = varname_and_subcol{2};
    full_var_string = [varname, ' ', num2str(subcolidx)];  
elseif isstring(varname_and_subcol) || ischar(varname_and_subcol)
    varname = varname_and_subcol;
    subcolidx = 1:size(tab{:,varname},2); % use all columns if not specified
    full_var_string = varname; 
else
    error('First argument must be 2-element cell or string')
end

% use all rows if not specified
vardefault('rows',1:height(tab));

vals_out = tab{:,varname}(rows,subcolidx); 