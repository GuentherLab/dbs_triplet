% Function to create a scatter plot for a single transition field (trans1 or trans2)
function plot_single_transition_scatter(reactiontime_table, transitions, pp_values, markers, trans_field, dbsID, plot_title)
    n_transitions = length(transitions);  % Number of transitions

    % Use a larger colormap with more distinct colors
    colormap_name = parula(n_transitions);  % Replace with 'parula', 'jet', 'hsv', or any other colormap with more colors

    figure;
    hold on;
    
    % Loop through each transition and plot the data
    for i = 1:n_transitions
        transition = transitions{i};
        
        % Get reaction times for this transition
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = reactiontime_table.reaction_time(valid_rows);
        
        % Skip transitions that don't have data for this DBS ID
        if isempty(rt_values)
            continue;
        end
        
        % Get the phonotactic probability value for the current transition
        pp_value = pp_values(i);
        
        % Plot the scatter plot for this transition
        scatter(pp_value * ones(size(rt_values)), rt_values, 100, colormap_name(i, :), markers{mod(i - 1, length(markers)) + 1}, ...
            'filled', 'DisplayName', transition);  % Assign display name for legend
    end
    
    % Add labels, title, and legend
    xlabel('Phonotactic Probability', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Reaction Time (ms)', 'FontSize', 15, 'FontWeight', 'bold');
    title(plot_title, 'FontSize', 20, 'FontWeight', 'bold');
    legend('show', 'Location', 'bestoutside');  % Show legend with transition names
    grid on;
    
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.png', dbsID, trans_field));
end
