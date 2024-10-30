function plot_single_transition_modified(reactiontime_table, transitions, pp_values, colors, markers, trans_field, dbsID, plot_title)
    n_transitions = length(transitions);  % Number of transitions

    % Create a new figure with subplots
    figure;
    sgtitle(plot_title);  % Add a main title for the whole figure
    num_rows = ceil(sqrt(n_transitions));  % Rows of subplots
    num_cols = ceil(n_transitions / num_rows);  % Columns of subplots
    
    % Loop through each transition and plot the data in individual subplots
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
        
        % Get the phonotactic probability for this transition
        pp_value = pp_values(i);
        
        % Create a subplot for this transition
        subplot(num_rows, num_cols, i);
        hold on;
        
        % Plot with color coding and different symbols
        errorbar(pp_value, mean_rt, std_rt, markers{i}, 'Color', colors(i, :), 'MarkerSize', 8, ...
            'DisplayName', transition);  % Assign display name for legend
        
        % Add labels and title for the subplot
        xlabel('Phonotactic Probability');
        ylabel('Mean Reaction Time (ms)');
        title(transition);
        
        % Adjust x-axis limits: If pp_value is 0 or small, apply a small buffer
        xlim_buffer = 0.00005;  % Default buffer
        if pp_value <= 0
            xlim([-xlim_buffer, xlim_buffer * 10]);  % Adjust to give space for pp = 0
        else
            xlim([pp_value - xlim_buffer, pp_value + xlim_buffer]);
        end
        
        % Set font size for readability
        set(gca, 'FontSize', 12);
        grid on;
        hold off;
    end
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.png', dbsID, trans_field));
end
