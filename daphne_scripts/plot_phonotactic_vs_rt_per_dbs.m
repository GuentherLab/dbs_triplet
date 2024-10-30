% Function to plot phonotactic probability vs mean reaction time for each DBS ID with shift and log scale x-axis
function plot_phonotactic_vs_rt_per_dbs(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = nan(size(transitions));
    std_rt = nan(size(transitions));
    pp_values = zeros(size(transitions));  % Store phonotactic probabilities
    
    % Small shift to avoid issues with log scale
    small_shift = 0.001;
    
    % Colormap (Parula)
    colors = parula(length(transitions));  % Generate a Parula colormap with as many colors as there are transitions
    
    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
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
            pp_values(i) = phonotacticProbabilities.(transition) + small_shift;  % Add the small shift to avoid log scale issues
        end
    end
    
    % Skip if no valid transitions exist
    if all(isnan(mean_rt))
        disp(['No valid transitions for DBS ID: ', dbsID, ', skipping plot.']);
        return;  % Do not plot if there are no valid transitions
    end
    
    % Create the plot for this DBS ID
    figure;
    hold on;
    
    % Plot mean RT with error bars (std deviation) using scatter plot with different colors from the Parula colormap
    for i = 1:length(transitions)
        if ~isnan(mean_rt(i))
            % Plot with error bars
            errorbar(pp_values(i), mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
                'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
                sprintf('%s (%.4f)', transitions{i}, pp_values(i) - small_shift));  % Display original phonotactic probability in legend
        end
    end
    
    % Set the x-axis to be logarithmic
    set(gca, 'XScale', 'log');
    
    % Manually set ticks for better log-scale representation
    xticks([0.001 0.01 0.1 1 10]);  % Adjust the tick marks as needed
    
    % Add labels, title, and grid
    xlabel('Phonotactic Probability (Log Scale)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Mean Reaction Time (s)', 'FontSize', 15, 'FontWeight', 'bold');
    
    % Update the title based on whether it's trans1 or trans2
    if strcmp(trans_field, 'trans1')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for First Transitions - DBS ID: %s', dbsID);
    elseif strcmp(trans_field, 'trans2')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for Second Transitions - DBS ID: %s', dbsID);
    end
    
    % Set the plot title
    title(title_str, 'FontSize', 20, 'FontWeight', 'bold');
    
    % Add a legend to identify the phonotactic probability for each transition
    legend('show', 'Location', 'bestoutside');
    
    grid on;
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.png', dbsID, trans_field));
end
