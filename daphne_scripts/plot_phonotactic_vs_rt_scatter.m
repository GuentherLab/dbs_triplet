% Function to plot scatter plot with phonotactic probabilities on the x-axis and reaction time on the y-axis
% Phonotactic probabilities will be shown on the x-axis, and transitions will be represented by different markers
function plot_phonotactic_vs_rt_scatter(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = nan(size(transitions));
    std_rt = nan(size(transitions));
    pp_values = nan(size(transitions));  % Array to store phonotactic probabilities
    
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
            pp_values(i) = phonotacticProbabilities.(transition);  % Store the phonotactic probability
            
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
    
    % Create the plot
    figure;
    hold on;
    
    % Plot mean RT with error bars (std deviation) using scatter plot with different colors from the Parula colormap
    for i = 1:length(valid_transitions)
        % Format the phonotactic probability for the legend (display as '0' if zero, otherwise four decimals)
        pp_value = valid_pp_values(i);
        if pp_value == 0
            pp_display = '0';
        else
            pp_display = sprintf('%.4f', pp_value);
        end
        
        % Plot with error bars
        errorbar(pp_values(i), mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
            sprintf('%s (%s)', valid_transitions{i}, pp_display));  % Legend shows transition and formatted phonotactic probability
    end
    
    % Set the x-axis to be logarithmic (if needed)
    set(gca, 'XScale', 'log');
    
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
    
    % Add a legend to identify the transition ID and phonotactic probability for each point
    legend('show', 'Location', 'bestoutside');
    
    grid on;
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_Scatter_%s_%s.png', dbsID, trans_field));
end
