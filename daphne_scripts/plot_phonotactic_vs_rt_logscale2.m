% Function to plot transitions on x-axis (using phonotactic probabilities) and mean reaction time with error bars
% Phonotactic Probabilities will be shown in the x-axis (log scale) and transitions as legend labels
% Function to plot transitions on x-axis (using phonotactic probabilities) and mean reaction time with error bars
% Phonotactic Probabilities will be shown in the x-axis (log scale) and transitions as legend labels
function plot_phonotactic_vs_rt_logscale2(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = nan(size(transitions));
    std_rt = nan(size(transitions));
    
    % Colormap (Parula)
    colors = parula(length(transitions));  % Generate a Parula colormap with as many colors as there are transitions
    
    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
    % Filter valid transitions (those with sufficient data)
    valid_transitions = {};
    valid_pp_values = [];
    
    % Loop through each transition and calculate mean and std RT
    for i = 1:length(transitions)
        transition = transitions{i};
        
        % Get reaction times for the transition
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = reactiontime_table.reaction_time(valid_rows);
        
        % Calculate mean and std RT for this transition if data is available
        if numel(rt_values) > 1  % Ensure there are enough data points
            mean_rt(i) = mean(rt_values);  % Store the mean RT
            std_rt(i) = std(rt_values);    % Store the std deviation
            
            % Add to valid transitions and phonotactic probability values
            valid_transitions{end+1} = transition;
            valid_pp_values(end+1) = phonotacticProbabilities.(transition);
        end
    end
    
    % Skip if no valid transitions exist
    if isempty(valid_transitions)
        disp(['No valid transitions for DBS ID: ', dbsID, ', skipping plot.']);
        return;  % Do not plot if there are no valid transitions
    end
    
    % Sort the transitions and phonotactic probabilities in ascending order
    [sorted_pp_values, sort_idx] = sort(valid_pp_values);  % Sort phonotactic probabilities
    sorted_transitions = valid_transitions(sort_idx);  % Re-order transitions based on sorted indices
    
    % Create the plot
    figure;
    hold on;
    
    % Replace zero phonotactic probabilities with a small positive number
    sorted_pp_values(sorted_pp_values == 0) = eps;
    
    % Plot mean RT with error bars (std deviation) using scatter plot with different colors from the Parula colormap
    for i = 1:length(sorted_transitions)
        % Format the phonotactic probability for the legend (display as '0' if zero, otherwise four decimals)
        pp_value = sorted_pp_values(i);
        if pp_value == eps
            pp_display = '0';  % Show "0" for probabilities that were originally zero
        else
            pp_display = sprintf('%.4f', pp_value);
        end
        
        % Plot with error bars at the phonotactic probability (x-axis) and mean RT (y-axis)
        errorbar(pp_value, mean_rt(sort_idx(i)), std_rt(sort_idx(i)), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
            sprintf('%s', sorted_transitions{i}));  % Legend shows transition ID
    end
    
    % Add labels, title, and grid
    xlabel('Phonotactic Probability (Log Scale)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Mean Reaction Time (s)', 'FontSize', 15, 'FontWeight', 'bold');
    
    % Update the title based on whether it's trans1 or trans2
    if strcmp(trans_field, 'trans1')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for First Transitions in Subject: %s', dbsID);
    elseif strcmp(trans_field, 'trans2')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for Second Transitions in Subject: %s', dbsID);
    end
    
    % Set the plot title
    title(title_str, 'FontSize', 20, 'FontWeight', 'bold');
    
    % Set x-axis to log scale
    set(gca, 'XScale', 'log');
    
    % Adjust x-axis limits to avoid cluttering near very small values
    xlim([min(sorted_pp_values)*0.8, max(sorted_pp_values)*1.2]);  % Add buffer
    
    % Add a legend to identify the transition ID for each point
    legend('show', 'Location', 'bestoutside');
    
    grid on;
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_LogScale_%s_%s.png', dbsID, trans_field));
end



