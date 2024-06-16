% Helper function to convert a table to a tab-delimited string
function clipboardString = convertTableToClipboardString(table)
    % Convert the table to a cell array
    cellArray = table2cell(table);
    
    % Add headers
    headers = table.Properties.VariableNames;
    cellArray = [headers; cellArray];
    
    % Convert the cell array to a tab-delimited string
    clipboardString = '';
    for i = 1:size(cellArray, 1)
        rowString = strjoin(cellArray(i, :), '\t');
        clipboardString = [clipboardString rowString '\n'];
    end
end
