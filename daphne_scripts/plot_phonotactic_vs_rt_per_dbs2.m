% Function to plot transitions with phonotactic probabilities (shifted) and reaction times in log scale
function plot_phonotactic_vs_rt_per_dbs2(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
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
            % Add 0.0001 to the phonotactic probability for better log scale visualization
            valid_pp_values(end+1) = phonotacticProbabilities.(transition) + 0.0001;  
        else
            % Debugging message for transitions with insufficient data
            disp(['No valid data for transition: ', transition]);
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
        pp_display = sprintf('%.4f', pp_value);  % Ensure correct formatting of the phonotactic probability

        % Plot with error bars
        errorbar(pp_value, mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
            sprintf('%s (%s)', valid_transitions{i}, pp_display));  % Legend shows transition and formatted phonotactic probability
    end

    % Add labels, title, and grid
    xlabel('Phonotactic Probability (Log Scale)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Mean Reaction Time (s)', 'FontSize', 15, 'FontWeight', 'bold');
    
    % Set the plot title
    if strcmp(trans_field, 'trans1')
        title_str = sprintf('Mean Reaction Time for First Transitions - DBS ID: %s', dbsID);
    else
        title_str = sprintf('Mean Reaction Time for Second Transitions - DBS ID: %s', dbsID);
    end
    title(title_str, 'FontSize', 20, 'FontWeight', 'bold');

    % Add a legend to identify the phonotactic probability for each transition
    legend('show', 'Location', 'bestoutside');

    % Set the x-axis to log scale and define custom limits
    set(gca, 'XScale', 'log');
    xlim([1e-4, 1e-2]);  % Set the x-axis limits between 0.0001 and 0.01 (log scale)

    % Add grid and finish plot
    grid on;
    hold off;

    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.png', dbsID, trans_field));
end
