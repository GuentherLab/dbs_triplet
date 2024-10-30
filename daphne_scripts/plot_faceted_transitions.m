%% % Function to create a faceted plot for each transition
function plot_faceted_transitions(reactiontime_table, transitions, pp_values, markers, trans_field, dbsID, plot_title)
    n_transitions = length(transitions);  % Number of transitions
    num_rows = ceil(sqrt(n_transitions));  % Calculate rows for subplots
    num_cols = ceil(n_transitions / num_rows);  % Calculate columns for subplots

    figure;
    sgtitle(plot_title);  % Main title for the whole figure
    
    % Loop through each transition and create a subplot for each
    for i = 1:n_transitions
        subplot(num_rows, num_cols, i);  % Create a subplot for this transition
        hold on;
        
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
        scatter(pp_value * ones(size(rt_values)), rt_values, 100, markers{mod(i - 1, length(markers)) + 1}, 'filled', ...
            'DisplayName', transition);
        
        % Add labels and title for each subplot
        xlabel('Phonotactic Probability');
        ylabel('Reaction Time');
        title(transition);
        grid on;
        hold off;
    end
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Faceted_Transitions_%s_%s.png', dbsID, trans_field));
end

