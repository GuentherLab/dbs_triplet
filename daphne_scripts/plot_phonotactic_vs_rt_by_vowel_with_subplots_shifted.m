% Function to create subplots for each vowel-based transition group
function plot_phonotactic_vs_rt_by_vowel_with_subplots_shifted(reactiontime_table, phonotacticProbabilities, transitions, dbsID, trans_field)
    % Correctly initialize vowel_groups with empty cell arrays for 'a', 'o', and 'e'
    vowel_groups.a = {};
    vowel_groups.o = {};
    vowel_groups.e = {};
    
    % Small shift to apply to the phonotactic probabilities
    small_shift = 0.0001;
    
    % Group transitions based on the first vowel
    for i = 1:length(transitions)
        transition = transitions{i};
        first_vowel = transition(1);  % Extract the first letter of the transition
        
        switch first_vowel
            case 'a'
                vowel_groups.a{end+1} = transition;
            case 'o'
                vowel_groups.o{end+1} = transition;
            case 'e'
                vowel_groups.e{end+1} = transition;
        end
    end

    % Create subplots for each vowel group
    figure;
    
    % Define marker styles for the plot
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*', 'x'};
    
    vowels = fieldnames(vowel_groups);  % Get the field names ('a', 'o', 'e')
    
    for v = 1:length(vowels)
        vowel = vowels{v};
        transitions_for_vowel = vowel_groups.(vowel);
        
        % Prepare arrays to store mean and std RTs, phonotactic probabilities
        mean_rt = nan(size(transitions_for_vowel));
        std_rt = nan(size(transitions_for_vowel));
        pp_values = zeros(size(transitions_for_vowel));  % Store phonotactic probabilities
        
        % Calculate mean, std RTs, and sort by phonotactic probabilities
        for i = 1:length(transitions_for_vowel)
            transition = transitions_for_vowel{i};
            
            % Get reaction times for the transition
            valid_rows = strcmp(reactiontime_table.(trans_field), transition);
            rt_values = reactiontime_table.reaction_time(valid_rows);
            
            % Store mean and std RTs for valid transitions
            if numel(rt_values) > 1
                mean_rt(i) = mean(rt_values);
                std_rt(i) = std(rt_values);
                % Apply the small shift to phonotactic probabilities
                pp_values(i) = phonotacticProbabilities.(transition) + small_shift;
            end
        end
        
        % Sort by phonotactic probability
        [pp_values, sort_idx] = sort(pp_values);
        mean_rt = mean_rt(sort_idx);
        std_rt = std_rt(sort_idx);
        transitions_for_vowel = transitions_for_vowel(sort_idx);
        
        % Create a subplot for this vowel group
        subplot(1, length(vowels), v);
        hold on;
        
        % Plot mean RT with error bars (std deviation)
        for i = 1:length(transitions_for_vowel)
            pp_value = pp_values(i);
            marker_style = markers{mod(i-1, length(markers)) + 1};
            
            % Plot each transition
            errorbar(pp_value, mean_rt(i), std_rt(i), marker_style, ...
                'MarkerSize', 12, 'MarkerFaceColor', 'auto', 'DisplayName', ...
                sprintf('%s (%.4f)', transitions_for_vowel{i}, pp_value));
        end
        
        % Set x-axis to log scale and add labels
        set(gca, 'XScale', 'log');
        xlabel('Phonotactic Probability (Log Scale)', 'FontSize', 12);
        ylabel('Mean Reaction Time (s)', 'FontSize', 12);
        title(sprintf('Transitions starting with "%s"', vowel), 'FontSize', 14);
        
        % Add grid and legend
        grid on;
        legend('show', 'Location', 'bestoutside');
        
        hold off;
    end
    
    % Add a main title for the entire figure
    sgtitle(sprintf('Phonotactic Probability vs Reaction Time by Vowel Group - DBS ID: %s', dbsID), 'FontSize', 16);
end
