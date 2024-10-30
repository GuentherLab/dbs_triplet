% Function to plot phonotactic probabilities vs reaction time with log scale, jitter, and a linear fit line
function plot_phonotactic_vs_rt_per_dbs4(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)

    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = [];
    std_rt = [];
    valid_pp_values = [];
    valid_transitions = {};

    % Adjust phonotactic probabilities by adding 0.0001 to all transitions to avoid log(0)
    adjustment = 0.0001;

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
            valid_transitions{end+1} = transition;
            valid_pp_values(end+1) = phonotacticProbabilities.(transition) + adjustment;  % Avoid log(0) issue
        end
    end

    % Diagnostic message to check if data exists
    if isempty(valid_transitions)
        disp(['No valid transitions for DBS ID: ', dbsID, ', skipping plot.']);
        return;
    else
        disp('Valid transitions found. Proceeding to plot.');
    end

    % Check if we have valid data to fit a line
    if isempty(valid_pp_values) || isempty(mean_rt)
        disp('No valid data for linear fit.');
        return;
    end

    % Manually define a set of distinct colors
    colors = lines(length(valid_transitions));

    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};

    % Create the plot
    figure;
    hold on;

    % Scatter plot the data points with jitter
    for i = 1:length(valid_transitions)
        jittered_pp = valid_pp_values(i) + randn * 0.00005;  % Add jitter
        scatter(jittered_pp, mean_rt(i), 80, colors(i, :), 'filled');  % Scatter plot
    end

    % Print values to verify
    disp('Valid pp values:');
    disp(valid_pp_values);
    disp('Mean RT values:');
    disp(mean_rt);

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
plot(exp(log_pp_values), yfit, 'k-', 'LineWidth', 3.0);  % Red line, adjusted for log scale


    % Add labels, title, and grid
    xlabel('Phonotactic Probability (Log Scale)', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Mean Reaction Time (s)', 'FontSize', 15, 'FontWeight', 'bold');

    % Update the title based on whether it's trans1 or trans2
    if strcmp(trans_field, 'trans1')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for First Transitions - DBS ID: %s', dbsID);
    elseif strcmp(trans_field, 'trans2')
        title_str = sprintf('Phonotactic Probability vs Mean Reaction Time for Second Transitions - DBS ID: %s', dbsID);
    end
    title(title_str, 'FontSize', 20, 'FontWeight', 'bold');

    % Set x-axis to log scale
    set(gca, 'XScale', 'log');

    % Set x-axis limits with a buffer
    xlim([min(valid_pp_values) - 0.00005, max(valid_pp_values) + 0.0002]);

    grid on;
    hold off;

    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_%s_%s.jpeg', dbsID, trans_field));

end
