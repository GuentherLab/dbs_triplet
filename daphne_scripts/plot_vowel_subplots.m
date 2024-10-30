function plot_vowel_subplots(vowel, significant_dbsIDs, trans_field, vowel_groups, plot_title, subtable)
    transitions = vowel_groups.(vowel);  % Extract transitions for the given vowel
    n_subplots = length(transitions);
    
    % Create a figure for subplots
    figure;
    sgtitle(plot_title);  % Add a main title for the subplot figure
    
    % Define the transition strings and their phonotactic probabilities
    phonotacticProbabilities = struct('aht', 0.0001, 'ahv', 0.0001, 'ahgh', 0, 'ahs', 0.0007, ...
                                      'oot', 0, 'oov', 0, 'oogh', 0, 'oos', 0, ...
                                      'eet', 0.0002, 'eev', 0.0007, 'eegh', 0.0005, 'ees', 0.0003);
    
    % Preallocate pp_values
    pp_values = zeros(1, n_subplots);
    
    % Retrieve phonotactic probabilities for each transition
    for i = 1:n_subplots
        transition = transitions{i};
        pp_values(i) = phonotacticProbabilities.(transition);  % Get the phonotactic probability for each transition
    end
    
    % Define a colormap and marker options for each DBS ID
    colors = lines(numel(significant_dbsIDs));  % One color for each DBS ID
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};  % Marker styles
    
    % Create subplots for each transition
    for i = 1:n_subplots
        subplot(ceil(sqrt(n_subplots)), ceil(sqrt(n_subplots)), i);  % Arrange subplots in a grid
        hold on;
        
        transition = transitions{i};
        
        % Loop over DBS IDs to plot each one with a distinct color and marker
        for i_sub = 1:numel(significant_dbsIDs)
            dbsID = significant_dbsIDs{i_sub};
            reactiontime_table = subtable.reactiontime_tables{i_sub};

            % Check if reactiontime_table is empty or not a table, and skip if so
            if isempty(reactiontime_table) || ~istable(reactiontime_table)
                continue;
            end
            
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
            
            % Add a small jitter to the phonotactic probability values to avoid overlap
            jittered_pp = pp_values(i) + randn * 0.0001;  % Adjust the jitter magnitude as needed
            
            % Plot with distinct color and marker for each DBS ID
            errorbar(jittered_pp, mean_rt, std_rt, markers{i_sub}, 'Color', colors(i_sub, :), 'MarkerSize', 8, ...
                'DisplayName', dbsID);  % Assign DBS ID to the legend
        end
        
        % Add labels, title, and adjust axis limits
        xlabel('Phonotactic Probability');
        ylabel('Mean Reaction Time (ms)');
        title(transition);
        
        % Adjust x-axis limits with a buffer, but check that the range is valid
        if max(pp_values) ~= min(pp_values)
            xlim_buffer = (max(pp_values) - min(pp_values)) * 0.1;  % 10% of the range
            xlim([min(pp_values) - xlim_buffer, max(pp_values) + xlim_buffer]);
        else
            % If min and max are the same, set a small default range
            xlim([pp_values(i) - 0.0001, pp_values(i) + 0.0001]);
        end
        
        % Set font size for readability
        set(gca, 'FontSize', 12);  % Increase font size for axes
        
        % Show legend with DBS IDs
        legend('show', 'Location', 'bestoutside');  % Show legend with DBS IDs
        grid on;
        hold off;
    end
end
