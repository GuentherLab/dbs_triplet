% Function to create histogram-like plots based on individual transitions
function plot_phonotactic_vs_rt_histogram_individual(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = nan(size(transitions));
    std_rt = nan(size(transitions));
    
    % Colormap (Parula)
    colors = parula(length(transitions));  % Generate a Parula colormap with as many colors as there are transitions
    
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
    
    % Create the plot
    figure;
    hold on;
    
    % Plot histogram with error bars for each transition
    for i = 1:length(valid_transitions)
        % Get mean reaction time and std for this transition
        pp_value = valid_pp_values(i);
        mean_value = mean_rt(i);
        std_value = std_rt(i);
        
        % Bar plot for mean reaction time
        bar(pp_value, mean_value, 'FaceColor', colors(i, :), 'EdgeColor', colors(i, :));
        
        % Error bars representing the standard deviation
        errorbar(pp_value, mean_value, std_value, 'k', 'LineWidth', 1.5);
    end
    
    % Set the x-axis to represent phonotactic probability (log scale for better representation)
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
    
    grid on;
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Phonotactic_Probability_vs_RT_Histogram_%s_%s.png', dbsID, trans_field));
end
