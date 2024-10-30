% Function to create a box plot for each transition
function plot_boxplot_transitions(reactiontime_table, transitions, trans_field, dbsID, plot_title, phonotacticProbabilities)
    rt_data = [];  % Store all reaction times
    pp_labels = {};  % Store corresponding labels for the box plot

    figure;
    hold on;

    % Loop through each transition and gather data
    for i = 1:length(transitions)
        transition = transitions{i};
        
        % Get reaction times for this transition
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = reactiontime_table.reaction_time(valid_rows);
        
        % Skip transitions that don't have data for this DBS ID
        if isempty(rt_values)
            continue;
        end
        
        % Append reaction times and labels for the box plot
        rt_data = [rt_data; rt_values];
        
        % Get the phonotactic probability for this transition
        pp_value = phonotacticProbabilities.(transition);
        
        % Create the label with the phonotactic probability in parentheses
        pp_label = sprintf('%s (%.4f)', transition, pp_value);  % Add the probability with 4 decimal places
        
        % Add label for this transition
        pp_labels = [pp_labels; repmat({pp_label}, length(rt_values), 1)];
    end
    
    % Create the box plot
    boxplot(rt_data, pp_labels);
    
    % Add labels, title, and grid
    xlabel('Transitions', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Reaction Time (ms)', 'FontSize', 15, 'FontWeight', 'bold');
    title(plot_title, 'FontSize', 20, 'FontWeight', 'bold');
    grid on;

    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Boxplot_Transitions_%s_%s.png', dbsID, trans_field));
end
