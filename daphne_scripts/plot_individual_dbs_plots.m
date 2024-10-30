function plot_individual_dbs_plots(vowel, significant_dbsIDs, trans_field, vowel_groups, plot_title, subtable)
    transitions = vowel_groups.(vowel);  % Extract transitions for the given vowel
    n_transitions = length(transitions);
    
    % Loop over each DBS ID to create individual plots
    for i_sub = 1:numel(significant_dbsIDs)
        dbsID = significant_dbsIDs{i_sub};
        reactiontime_table = subtable.reactiontime_tables{i_sub};
        
        % Check if reactiontime_table is empty or not a table, and skip if so
        if isempty(reactiontime_table) || ~istable(reactiontime_table)
            continue;
        end
        
        % Check if there's any valid data for this DBS ID across all transitions
        valid_transitions_count = 0;  % Counter to track how many transitions have data
        for i = 1:n_transitions
            transition = transitions{i};
            valid_rows = strcmp(reactiontime_table.(trans_field), transition);
            rt_values = reactiontime_table.reaction_time(valid_rows);
            if ~isempty(rt_values)
                valid_transitions_count = valid_transitions_count + 1;
            end
        end
        
        % Skip this DBS ID if no valid transitions found
        if valid_transitions_count == 0
            disp(['No significant transitions found for DBS ID: ', dbsID, ', skipping plot.']);
            continue;
        end
        
        % If there's valid data, proceed to create the plot
        figure;
        sgtitle([plot_title, ' (DBS ID: ', dbsID, ')']);  % Main title with DBS ID
        
        % Define colors and markers
        colors = lines(n_transitions);  % One color for each transition
        markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};  % Marker styles
        
        % Create subplots for each transition
        for i = 1:n_transitions
            subplot(ceil(sqrt(n_transitions)), ceil(sqrt(n_transitions)), i);  % Arrange subplots in a grid
            hold on;
            
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
            
            % Define the phonotactic probability values
            phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                              'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                              'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
            pp_value = phonotacticProbabilities.(transition);  % Get phonotactic probability for this transition
            
            % Add a small jitter to the phonotactic probability values to avoid overlap
            jittered_pp = pp_value + randn * 0.0001;  % Adjust the jitter magnitude as needed
            
            % Plot with color coding and different symbols for each transition
            errorbar(jittered_pp, mean_rt, std_rt, markers{i}, 'Color', colors(i, :), 'MarkerSize', 8, ...
                'DisplayName', transition);  % Label for each transition
            
            % Add labels and title for each subplot
            xlabel('Phonotactic Probability');
            ylabel('Mean Reaction Time (ms)');
            title(transition);
            
            % Adjust x-axis limits: If the range is zero, apply a small buffer
            xlim_buffer = 0.00005;  % Set a small default buffer
            
            if pp_value > 0
                xlim([pp_value - xlim_buffer, pp_value + xlim_buffer]);
            else
                xlim([-xlim_buffer, xlim_buffer]);  % Handle zero or near-zero values
            end
            
            % Set font size for readability
            set(gca, 'FontSize', 12);
            grid on;
            hold off;
        end
    end
end
