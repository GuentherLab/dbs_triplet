function plot_phonotactic_vs_rt_per_dbs5(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)

    % Adjust phonotactic probabilities by adding 0.0001 to all transitions to avoid log(0)
    adjusted_phonotacticProbabilities = phonotacticProbabilities;
    adjustment = 0.0001;  % Add this to all probabilities
    transition_fields = fieldnames(phonotacticProbabilities);
    
    for i = 1:numel(transition_fields)
        adjusted_phonotacticProbabilities.(transition_fields{i}) = phonotacticProbabilities.(transition_fields{i}) + adjustment;
    end
    
    % Initialize arrays to store valid transitions and phonotactic probabilities
    valid_transitions = {};
    valid_pp_values = [];
    mean_rt = [];
    std_rt = [];

    % Loop through each transition and calculate mean and std RT
    for i = 1:length(transitions)
        transition = transitions{i};
        
        % Get reaction times for the transition
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = reactiontime_table.reaction_time(valid_rows);
        
        % Calculate mean and std RT for this transition if data is available
        if numel(rt_values) > 1  % Ensure there are enough data points
            mean_rt(end+1) = abs(mean(rt_values));  % Store the positive mean RT
            std_rt(end+1) = std(rt_values);    % Store the std deviation
            
            % Add to valid transitions and adjusted phonotactic probability values
            valid_transitions{end+1} = transition;
            valid_pp_values(end+1) = phonotacticProbabilities.(transition) + 0.0001;  % Adding 0.0001 to phonotactic probabilities to avoid log issues
        end
    end

    % Skip if no valid transitions exist
    if isempty(valid_transitions)
        disp(['No valid transitions for DBS ID: ', dbsID, ', skipping plot.']);
        return;  % Do not plot if there are no valid transitions
    end

    % Manually define a set of distinct colors
    colors = [...
        0 0.4470 0.7410;  % Blue
        0.8500 0.3250 0.0980;  % Orange
        0.9290 0.6940 0.1250;  % Yellow (replaced below)
        0.4940 0.1840 0.5560;  % Purple
        0.4660 0.6740 0.1880;  % Green
        0.3010 0.7450 0.9330;  % Light Blue
        0.6350 0.0780 0.1840;  % Dark Red
        0.5 0.5 0;  % Olive
        0.75 0.75 0;  % Yellow-green
        0.75 0 0.75;  % Magenta
        0.25 0.25 0.75;  % Navy Blue
        0.83 0.65 0.15];  % Mustard Yellow

    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
    % Create the plot
    figure;
    hold on;

    % Scatter plot your data points with jitter
    for i = 1:length(valid_transitions)  % Use valid transitions for the plot
        % Apply jitter to phonotactic probabilities
        jittered_pp = valid_pp_values(i) + randn * 0.00005;  % Increase the jitter magnitude for better visualization
        
        % Plot with error bars
        errorbar(jittered_pp, mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
            valid_transitions{i});  % Only display the transition in the legend
    end

    % ---- Add the linear fit line ----
    % Remove NaN values from valid_pp_values and mean_rt
    valid_idx = ~isnan(valid_pp_values) & ~isnan(mean_rt);
    filtered_pp_values = valid_pp_values(valid_idx);
    filtered_mean_rt = mean_rt(valid_idx);
    
    % Log-transform the phonotactic probabilities for the linear fit
    log_pp_values = log(filtered_pp_values);

    % Add the linear fit line
    p = polyfit(log_pp_values, filtered_mean_rt, 1);  % Fit a linear model to log-transformed pp values
    yfit = polyval(p, log_pp_values);  % Get the fitted values

    % Plot the fitted line (re-convert log_pp_values to the original scale for plotting)
    fitLine = plot(exp(log_pp_values), yfit, 'k-', 'LineWidth', 3.0);

    % Remove the line from the legend explicitly
    set(get(get(fitLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

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

    % Add a legend to identify the transitions (only for scatter points)
    legend('show', 'Location', 'bestoutside');

    % Set x-axis to log scale
    set(gca, 'XScale', 'log');

    % Set x-axis limits with a buffer
    xlim([min(valid_pp_values) - 0.00005, max(valid_pp_values) + 0.0002]);

    grid on;
    hold off;

  function plot_phonotactic_vs_rt_per_dbs5(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)

    % Adjust phonotactic probabilities by adding 0.0001 to all transitions to avoid log(0)
    adjusted_phonotacticProbabilities = phonotacticProbabilities;
    adjustment = 0.0001;  % Add this to all probabilities
    transition_fields = fieldnames(phonotacticProbabilities);
    
    for i = 1:numel(transition_fields)
        adjusted_phonotacticProbabilities.(transition_fields{i}) = phonotacticProbabilities.(transition_fields{i}) + adjustment;
    end
    
    % Initialize arrays to store valid transitions and phonotactic probabilities
    valid_transitions = {};
    valid_pp_values = [];
    mean_rt = [];
    std_rt = [];

    % Loop through each transition and calculate mean and std RT
    for i = 1:length(transitions)
        transition = transitions{i};
        
        % Get reaction times for the transition
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = reactiontime_table.reaction_time(valid_rows);
        
        % Calculate mean and std RT for this transition if data is available
        if numel(rt_values) > 1  % Ensure there are enough data points
            mean_rt(end+1) = abs(mean(rt_values));  % Store the positive mean RT
            std_rt(end+1) = std(rt_values);    % Store the std deviation
            
            % Add to valid transitions and adjusted phonotactic probability values
            valid_transitions{end+1} = transition;
            valid_pp_values(end+1) = phonotacticProbabilities.(transition) + 0.0001;  % Adding 0.0001 to phonotactic probabilities to avoid log issues
        end
    end

    % Skip if no valid transitions exist
    if isempty(valid_transitions)
        disp(['No valid transitions for DBS ID: ', dbsID, ', skipping plot.']);
        return;  % Do not plot if there are no valid transitions
    end

    % Manually define a set of distinct colors
    colors = [...
        0 0.4470 0.7410;  % Blue
        0.8500 0.3250 0.0980;  % Orange
        0.9290 0.6940 0.1250;  % Yellow (replaced below)
        0.4940 0.1840 0.5560;  % Purple
        0.4660 0.6740 0.1880;  % Green
        0.3010 0.7450 0.9330;  % Light Blue
        0.6350 0.0780 0.1840;  % Dark Red
        0.5 0.5 0;  % Olive
        0.75 0.75 0;  % Yellow-green
        0.75 0 0.75;  % Magenta
        0.25 0.25 0.75;  % Navy Blue
        0.83 0.65 0.15];  % Mustard Yellow

    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
    % Create the plot
    figure;
    hold on;

    % Scatter plot your data points with jitter
    for i = 1:length(valid_transitions)  % Use valid transitions for the plot
        % Apply jitter to phonotactic probabilities
        jittered_pp = valid_pp_values(i) + randn * 0.00005;  % Increase the jitter magnitude for better visualization
        
        % Plot with error bars
        errorbar(jittered_pp, mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), 'DisplayName', ...
            valid_transitions{i});  % Only display the transition in the legend
    end

    % ---- Add the linear fit line ----
    % Remove NaN values from valid_pp_values and mean_rt
    valid_idx = ~isnan(valid_pp_values) & ~isnan(mean_rt);
    filtered_pp_values = valid_pp_values(valid_idx);
    filtered_mean_rt = mean_rt(valid_idx);
    
    % Log-transform the phonotactic probabilities for the linear fit
    log_pp_values = log(filtered_pp_values);

    % Add the linear fit line
    p = polyfit(log_pp_values, filtered_mean_rt, 1);  % Fit a linear model to log-transformed pp values
    yfit = polyval(p, log_pp_values);  % Get the fitted values

    % Plot the fitted line (re-convert log_pp_values to the original scale for plotting)
    fitLine = plot(exp(log_pp_values), yfit, 'k-', 'LineWidth', 3.0);

    % Remove the line from the legend explicitly
    set(get(get(fitLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

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

    % Add a legend to identify the transitions (only for scatter points)
    legend('show', 'Location', 'bestoutside');

    % Set x-axis to log scale
    set(gca, 'XScale', 'log');

    % Set x-axis limits with a buffer
    xlim([min(valid_pp_values) - 0.00005, max(valid_pp_values) + 0.0002]);

    grid on;
    hold off;

    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.jpeg', dbsID, trans_field));

end

end
