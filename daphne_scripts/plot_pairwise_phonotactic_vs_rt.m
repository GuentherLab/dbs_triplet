function plot_pairwise_phonotactic_vs_rt(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Prepare arrays to store mean and std RTs for the transitions
    mean_rt = nan(size(transitions));
    std_rt = nan(size(transitions));
    
    % Colormap (Parula)
    colors = parula(length(transitions));  % Generate a Parula colormap with as many colors as there are transitions
    
    % Marker styles for transitions
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
    % Filter valid transitions and assign to different probability ranges
    transitions_group_1 = {}; % 0.0000 - 0.0005
    transitions_group_2 = {}; % 0.0005 - 0.0010
    
    for i = 1:length(transitions)
        transition = transitions{i};
        pp_value = phonotacticProbabilities.(transition);
        
        % Get reaction times for the transition and take the absolute values
        valid_rows = strcmp(reactiontime_table.(trans_field), transition);
        rt_values = abs(reactiontime_table.reaction_time(valid_rows));  % Take absolute reaction times
        
        % Calculate mean and std RT for this transition if data is available
        if numel(rt_values) > 1  % Ensure there are enough data points
            mean_rt(i) = mean(rt_values);  % Store the mean RT (absolute)
            std_rt(i) = std(rt_values);    % Store the std deviation (absolute)
            
            % Assign transitions to their respective ranges
            if pp_value >= 0 && pp_value <= 0.0005
                transitions_group_1{end+1} = transition;
            elseif pp_value > 0.0005 && pp_value <= 0.0010
                transitions_group_2{end+1} = transition;
            end
        end
    end
    
    % Create figure for pairwise plots
    figure;
    
    % Subplot 1: Phonotactic probabilities between 0.0000 - 0.0005
    subplot(1, 2, 1);
    hold on;
    for i = 1:length(transitions_group_1)
        transition = transitions_group_1{i};
        pp_value = phonotacticProbabilities.(transition);
        
        % Plot with error bars
        errorbar(pp_value + randn * 0.00001, mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), ...
            'DisplayName', transition);  % Legend shows transition and formatted phonotactic probability
    end
    title('Phonotactic Probabilities (0 - 0.0005)');
    xlabel('Phonotactic Probability (Log Scale)');
    ylabel('Mean Reaction Time (s)');
    legend('show', 'Location', 'bestoutside');
    set(gca, 'XScale', 'log');
    grid on;
    hold off;
    
    % Subplot 2: Phonotactic probabilities between 0.0005 - 0.0010
    subplot(1, 2, 2);
    hold on;
    for i = 1:length(transitions_group_2)
        transition = transitions_group_2{i};
        pp_value = phonotacticProbabilities.(transition);
        
        % Plot with error bars
        errorbar(pp_value + randn * 0.00001, mean_rt(i), std_rt(i), markers{mod(i-1, length(markers)) + 1}, ...
            'MarkerSize', 12, 'MarkerFaceColor', colors(i, :), 'Color', colors(i, :), ...
            'DisplayName', transition);  % Legend shows transition and formatted phonotactic probability
    end
    title('Phonotactic Probabilities (0.0005 - 0.0010)');
    xlabel('Phonotactic Probability (Log Scale)');
    ylabel('Mean Reaction Time (s)');
    legend('show', 'Location', 'bestoutside');
    set(gca, 'XScale', 'log');
    grid on;
    hold off;
    
    % Optionally, save the figure
    saveas(gcf, sprintf('Pairwise_Phonotactic_Probability_vs_RT_Abs_%s_%s.png', dbsID, trans_field));
end
