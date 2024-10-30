% Function to create a plot for a single transition field (trans1 or trans2)
function plot_single_transition(reactiontime_table, transitions, pp_values, markers, trans_field, dbsID, plot_title)
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
        
        % Calculate mean and std for this transition
        mean_rt = mean(rt_values);
        std_rt = std(rt_values);
        
        % Get the phonotactic probability value for the current transition
        pp_value = pp_values(i);
        
        % Cycle through markers if we have more transitions than markers
        marker_idx = mod(i - 1, length(markers)) + 1;
        
        % Plot with color coding from the parula colormap and different symbols
        % Increased marker size and set transparency for better visibility
        errorbar(pp_value, mean_rt, std_rt, markers{marker_idx}, 'Color', colormap_name(i, :), 'MarkerSize', 10, ...
            'MarkerFaceColor', colormap_name(i, :), 'MarkerEdgeColor', colormap_name(i, :), 'DisplayName', transition);  % Assign display name for legend
    end
    
    % Add labels, title, and legend
    xlabel('Phonotactic Probability', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Mean Reaction Time (ms)', 'FontSize', 15, 'FontWeight', 'bold');
    title(plot_title, 'FontSize', 20, 'FontWeight', 'bold');
    legend('show', 'Location', 'bestoutside');  % Show legend with transition names
    grid on;

    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.png', dbsID, trans_field));
end
