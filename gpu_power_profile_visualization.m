%% GPU Power Monitoring Visualization with State Annotations
% This script creates publication-quality figures from GPU monitoring CSV files
% with detailed state annotations for training and inference phases

clear; close all; clc;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% File paths (update these to your CSV file locations)
training_csv = 'training_metrics_with_states.csv';
inference_csv = 'inference_metrics_with_states.csv';

% Color palettes (ColorBrewer)
colors_qualitative = {'#66c2a5', '#fc8d62', '#8da0cb'};
colors_sequential = {'#e5f5f9', '#99d8c9', '#2ca25f'};

% State colors for training
state_colors_train = struct();
state_colors_train.initiation = '#8da0cb';      % Blue
state_colors_train.forward = '#66c2a5';         % Green
state_colors_train.backward = '#fc8d62';        % Orange
state_colors_train.communication = '#e78ac3';   % Pink

% State colors for inference
state_colors_infer = struct();
state_colors_infer.processing = '#fc8d62';      % Orange
state_colors_infer.waiting = '#8da0cb';         % Blue

%% ========================================================================
%  LOAD DATA
%  ========================================================================

fprintf('Loading data...\n');

% Load training data
try
    train_data = readtable(training_csv);
    fprintf('  Training data: %d samples loaded\n', height(train_data));
catch
    error('Could not load training CSV file: %s', training_csv);
end

% Load inference data
try
    infer_data = readtable(inference_csv);
    fprintf('  Inference data: %d samples loaded\n', height(infer_data));
catch
    error('Could not load inference CSV file: %s', inference_csv);
end

%% ========================================================================
%  TRAINING PHASE VISUALIZATION
%  ========================================================================

fprintf('\nCreating training phase visualization...\n');

fig1 = figure('Position', [100, 100, 1600, 1000], 'Color', 'w');

% Convert time to seconds
time_train = train_data.time_ms / 1000;

% Create 2x2 subplot layout
subplot(2, 2, 1); hold on; grid on; box on;
subplot(2, 2, 2); hold on; grid on; box on;
subplot(2, 2, 3); hold on; grid on; box on;
subplot(2, 2, 4); hold on; grid on; box on;

%% --- POWER DRAW SUBPLOT ---
subplot(2, 2, 1); hold on; grid on; box on;

% Plot all data
plot(time_train, train_data.power_draw_w, 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 2);

% Find and mark representative states
states_to_mark = {
    'training_initiation', 'Initiation\n(Warmup)', state_colors_train.initiation;
    'epoch_1_batch_1_forward_pass', 'Forward Pass\n(Compute predictions)', state_colors_train.forward;
    'epoch_1_batch_1_backward_pass', 'Backward Pass\n(Compute gradients)', state_colors_train.backward;
    'epoch_1_batch_1_communication', 'Communication\n(GPU sync)', state_colors_train.communication;
};

for i = 1:size(states_to_mark, 1)
    state_name = states_to_mark{i, 1};
    state_label = states_to_mark{i, 2};
    state_color = hex2rgb(states_to_mark{i, 3});
    
    % Find samples with this state
    idx = strcmp(train_data.gpu_state, state_name);
    
    if any(idx)
        % Step 1: Find start of first occurrence
        first_idx = find(idx, 1, 'first');
        t_start = time_train(first_idx);
        p_start = train_data.power_draw_w(first_idx);
        
        % Step 2: Find end by detecting state transition
        end_idx = first_idx;
        while end_idx < length(train_data.gpu_state) && ...
              strcmp(train_data.gpu_state(end_idx + 1), state_name)
            end_idx = end_idx + 1;
        end
        
        t_end = time_train(end_idx);
        p_end = train_data.power_draw_w(end_idx);
        
        % Now place markers at start and end of SAME occurrence
        plot(t_start, p_start, 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', state_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold on;
        plot(t_end, p_end, '^', 'MarkerSize', 8, ...
             'MarkerFaceColor', state_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

        % Add annotation (text and arrow)
        % if i == 1  % Initiation - annotate from left
        %     text(t_marker + 0.3, p_marker + 2, sprintf(state_label), ...
        %          'FontSize', 18, 'FontName', 'Times New Roman', ...
        %          'Interpreter', 'tex', 'Color', state_color, 'FontWeight', 'bold');
        % elseif i == 2  % Forward pass - annotate from top
        %     annotation('textarrow', [0.25, 0.25], [0.82, 0.87], ...
        %               'String', sprintf(state_label), 'FontSize', 16, ...
        %               'FontName', 'Times New Roman', 'Color', state_color, ...
        %               'LineWidth', 1.5, 'HeadStyle', 'cback1');
        % elseif i == 3  % Backward pass - annotate from top
        %     annotation('textarrow', [0.30, 0.30], [0.75, 0.80], ...
        %               'String', sprintf(state_label), 'FontSize', 16, ...
        %               'FontName', 'Times New Roman', 'Color', state_color, ...
        %               'LineWidth', 1.5, 'HeadStyle', 'cback1');
        % else  % Communication - annotate from bottom
        %     text(t_marker + 0.2, p_marker - 2, sprintf(state_label), ...
        %          'FontSize', 18, 'FontName', 'Times New Roman', ...
        %          'Interpreter', 'tex', 'Color', state_color, 'FontWeight', 'bold');
        % end
    end
end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Power Draw (W)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Training: Power Consumption', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- TEMPERATURE SUBPLOT ---
subplot(2, 2, 2); hold on; grid on; box on;

plot(time_train, train_data.temperature_c, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark, 1)
%     state_name = states_to_mark{i, 1};
%     state_color = hex2rgb(states_to_mark{i, 3});
% 
%     idx = strcmp(train_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_train(first_idx), train_data.temperature_c(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Temperature ($^\circ$C)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Training: GPU Temperature', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- GPU UTILIZATION SUBPLOT ---
subplot(2, 2, 3); hold on; grid on; box on;

plot(time_train, train_data.utilization_pct, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark, 1)
%     state_name = states_to_mark{i, 1};
%     state_color = hex2rgb(states_to_mark{i, 3});
% 
%     idx = strcmp(train_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_train(first_idx), train_data.utilization_pct(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('GPU Utilization (\%)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Training: GPU Utilization', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- MEMORY USAGE SUBPLOT ---
subplot(2, 2, 4); hold on; grid on; box on;

memory_gb_train = train_data.memory_used_mb / 1024;
plot(time_train, memory_gb_train, 'Color', [0.4 0.2 0.8], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark, 1)
%     state_name = states_to_mark{i, 1};
%     state_color = hex2rgb(states_to_mark{i, 3});
% 
%     idx = strcmp(train_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_train(first_idx), memory_gb_train(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Memory Used (GB)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Training: GPU Memory', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

% Save figure
saveas(fig1, 'training_metrics_annotated.png');
saveas(fig1, 'training_metrics_annotated.fig');
fprintf('  Saved: training_metrics_annotated.png/fig\n');


%% --- MEMORY USAGE SUBPLOT (your last subplot) ---
subplot(2, 2, 4); hold on; grid on; box on;
% ... (your existing subplot 4 code) ...
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- CREATE SHARED LEGEND FOR ALL SUBPLOTS ---
% Create dummy handles for legend (not visible in any subplot)
h_legend = [];
legend_labels = {};

% Colors for each state
col_init = hex2rgb(state_colors_train.initiation);
col_forward = hex2rgb(state_colors_train.forward);
col_backward = hex2rgb(state_colors_train.backward);
col_comm = hex2rgb(state_colors_train.communication);

% Create invisible dummy plots for legend
h_legend(1) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_init, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{1} = 'Initiation';

h_legend(2) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_forward, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{2} = 'Forward Pass';

h_legend(3) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_backward, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{3} = 'Backward Pass';

h_legend(4) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_comm, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{4} = 'Communication';

% Add marker type indicators
h_legend(5) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{5} = 'State Start';

h_legend(6) = plot(NaN, NaN, '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels{6} = 'State End';

% Create legend in a good position
leg = legend(h_legend, legend_labels, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'Interpreter', 'latex', ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'NumColumns', 3);

% Optional: Adjust legend position manually if needed
% leg.Position = [0.15, 0.02, 0.7, 0.04];  % [left, bottom, width, height]

%% ========================================================================
%  INFERENCE PHASE VISUALIZATION
%  ========================================================================

fprintf('\nCreating inference phase visualization...\n');

fig2 = figure('Position', [150, 150, 1600, 1000], 'Color', 'w');

% Convert time to seconds
time_infer = infer_data.time_ms / 1000;

%% --- POWER DRAW SUBPLOT ---
subplot(2, 2, 1); hold on; grid on; box on;

% Plot all data
plot(time_infer, infer_data.power_draw_w, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

% Find and mark representative states
states_to_mark_infer = {
    'processing_batch_1', 'Processing Query 1\n(Inference compute)', state_colors_infer.processing;
    'waiting_for_queries', 'Waiting for Queries\n(GPU idle)', state_colors_infer.waiting;
};

for i = 1:size(states_to_mark_infer, 1)
    state_name = states_to_mark_infer{i, 1};
    state_label = states_to_mark_infer{i, 2};
    state_color = hex2rgb(states_to_mark_infer{i, 3});
    
    % Find samples with this state
    idx = strcmp(infer_data.gpu_state, state_name);
    
    if any(idx)
        % Step 1: Find start of first occurrence
        first_idx = find(idx, 1, 'first');
        t_start = time_infer(first_idx);
        p_start = infer_data.power_draw_w(first_idx);
        
        % Step 2: Find end by detecting state transition
        end_idx = first_idx;
        while end_idx < length(infer_data.gpu_state) && ...
              strcmp(infer_data.gpu_state(end_idx + 1), state_name)
            end_idx = end_idx + 1;
        end
        
        t_end = time_infer(end_idx);
        p_end = infer_data.power_draw_w(end_idx);
        
        % Now place markers at start and end of SAME occurrence
        plot(t_start, p_start, 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', state_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold on;
        plot(t_end, p_end, '^', 'MarkerSize', 8, ...
             'MarkerFaceColor', state_color, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

        % Add annotation with arrow
        % if i == 1  % Processing - annotate from top
        %     annotation('textarrow', [0.18, 0.18], [0.78, 0.85], ...
        %               'String', sprintf(state_label), 'FontSize', 18, ...
        %               'FontName', 'Times New Roman', 'Color', state_color, ...
        %               'LineWidth', 1.5, 'HeadStyle', 'cback1', 'FontWeight', 'bold');
        % else  % Waiting - annotate from bottom
        %     annotation('textarrow', [0.23, 0.23], [0.70, 0.63], ...
        %               'String', sprintf(state_label), 'FontSize', 18, ...
        %               'FontName', 'Times New Roman', 'Color', state_color, ...
        %               'LineWidth', 1.5, 'HeadStyle', 'cback1', 'FontWeight', 'bold');
        % end
    end
end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Power Draw (W)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Inference: Power Consumption', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- TEMPERATURE SUBPLOT ---
subplot(2, 2, 2); hold on; grid on; box on;

plot(time_infer, infer_data.temperature_c, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark_infer, 1)
%     state_name = states_to_mark_infer{i, 1};
%     state_color = hex2rgb(states_to_mark_infer{i, 3});
% 
%     idx = strcmp(infer_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_infer(first_idx), infer_data.temperature_c(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Temperature ($^\circ$C)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Inference: GPU Temperature', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- GPU UTILIZATION SUBPLOT ---
subplot(2, 2, 3); hold on; grid on; box on;

plot(time_infer, infer_data.utilization_pct, 'Color', [0.2 0.6 0.2], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark_infer, 1)
%     state_name = states_to_mark_infer{i, 1};
%     state_color = hex2rgb(states_to_mark_infer{i, 3});
% 
%     idx = strcmp(infer_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_infer(first_idx), infer_data.utilization_pct(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('GPU Utilization (\%)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Inference: GPU Utilization', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- MEMORY USAGE SUBPLOT ---
subplot(2, 2, 4); hold on; grid on; box on;

memory_gb_infer = infer_data.memory_used_mb / 1024;
plot(time_infer, memory_gb_infer, 'Color', [0.4 0.2 0.8], 'LineWidth', 2);

% Mark the same representative states
% for i = 1:size(states_to_mark_infer, 1)
%     state_name = states_to_mark_infer{i, 1};
%     state_color = hex2rgb(states_to_mark_infer{i, 3});
% 
%     idx = strcmp(infer_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_infer(first_idx), memory_gb_infer(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Memory Used (GB)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Inference: GPU Memory', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- CREATE SHARED LEGEND FOR ALL INFERENCE SUBPLOTS ---
% Create dummy handles for legend (not visible in any subplot)
h_legend_infer = [];
legend_labels_infer = {};

% Colors for each state
col_processing = hex2rgb(state_colors_infer.processing);
col_waiting = hex2rgb(state_colors_infer.waiting);

% Create invisible dummy plots for legend
h_legend_infer(1) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_processing, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels_infer{1} = 'Processing Query';

h_legend_infer(2) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', col_waiting, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels_infer{2} = 'Waiting for Queries';

% Add marker type indicators
h_legend_infer(3) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels_infer{3} = 'State Start';

h_legend_infer(4) = plot(NaN, NaN, '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend_labels_infer{4} = 'State End';

% Create legend in a good position
leg_infer = legend(h_legend_infer, legend_labels_infer, ...
    'FontSize', 24, ...
    'FontName', 'Times New Roman', ...
    'Interpreter', 'latex', ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', ...
    'NumColumns', 2);

% Optional: Adjust legend position manually if needed
% leg_infer.Position = [0.15, 0.02, 0.7, 0.04];  % [left, bottom, width, height]



% Save figure
saveas(fig2, 'inference_metrics_annotated.png');
saveas(fig2, 'inference_metrics_annotated.fig');
fprintf('  Saved: inference_metrics_annotated.png/fig\n');

%% ========================================================================
%  COMBINED POWER COMPARISON FIGURE
%  ========================================================================

fprintf('\nCreating combined power comparison...\n');

fig3 = figure('Position', [200, 200, 1800, 600], 'Color', 'w');

%% --- TRAINING POWER (LEFT) ---
subplot(1, 2, 1); hold on; grid on; box on;

% Plot with state-based coloring
unique_states_train = unique(train_data.gpu_state, 'stable');
for i = 1:length(unique_states_train)
    state = unique_states_train{i};
    idx = strcmp(train_data.gpu_state, state);
    
    % Determine color based on state type
    if contains(state, 'forward')
        color = hex2rgb(state_colors_train.forward);
    elseif contains(state, 'backward')
        color = hex2rgb(state_colors_train.backward);
    elseif contains(state, 'communication')
        color = hex2rgb(state_colors_train.communication);
    else
        color = [0.7 0.7 0.7];
    end
    
    plot(time_train(idx), train_data.power_draw_w(idx), '.', ...
         'Color', color, 'MarkerSize', 15);
end

% Add representative markers
% for i = 1:size(states_to_mark, 1)
%     state_name = states_to_mark{i, 1};
%     state_color = hex2rgb(states_to_mark{i, 3});
% 
%     idx = strcmp(train_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_train(first_idx), train_data.power_draw_w(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 2);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Power Draw (W)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Training Phase', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');

% Add legend
legend_entries = {'Forward Pass', 'Backward Pass', 'Communication'};
legend_colors = {state_colors_train.forward, state_colors_train.backward, ...
                 state_colors_train.communication};
h_legend = [];
for i = 1:length(legend_entries)
    h_legend(i) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
                       'MarkerFaceColor', hex2rgb(legend_colors{i}), ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end
legend(h_legend, legend_entries, 'FontSize', 24, 'FontName', 'Times New Roman', ...
       'Location', 'best', 'Interpreter', 'latex');

set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

%% --- INFERENCE POWER (RIGHT) ---
subplot(1, 2, 2); hold on; grid on; box on;

% Plot with state-based coloring
unique_states_infer = unique(infer_data.gpu_state, 'stable');
for i = 1:length(unique_states_infer)
    state = unique_states_infer{i};
    idx = strcmp(infer_data.gpu_state, state);
    
    % Determine color based on state type
    if contains(state, 'processing')
        color = hex2rgb(state_colors_infer.processing);
    elseif contains(state, 'waiting')
        color = hex2rgb(state_colors_infer.waiting);
    else
        color = [0.7 0.7 0.7];
    end
    
    plot(time_infer(idx), infer_data.power_draw_w(idx), '.', ...
         'Color', color, 'MarkerSize', 15);
end

% Add representative markers
% for i = 1:size(states_to_mark_infer, 1)
%     state_name = states_to_mark_infer{i, 1};
%     state_color = hex2rgb(states_to_mark_infer{i, 3});
% 
%     idx = strcmp(infer_data.gpu_state, state_name);
%     if any(idx)
%         first_idx = find(idx, 1, 'first');
%         plot(time_infer(first_idx), infer_data.power_draw_w(first_idx), 'o', ...
%              'MarkerSize', 8, 'MarkerFaceColor', state_color, ...
%              'MarkerEdgeColor', 'k', 'LineWidth', 2);
%     end
% end

xlabel('Time (s)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Power Draw (W)', 'FontSize', 30, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title('Inference Phase', 'FontSize', 32, 'FontName', 'Times New Roman', ...
      'Interpreter', 'latex', 'FontWeight', 'bold');

% Add legend
legend_entries_infer = {'Processing Query', 'Waiting for Query'};
legend_colors_infer = {state_colors_infer.processing, state_colors_infer.waiting};
h_legend_infer = [];
for i = 1:length(legend_entries_infer)
    h_legend_infer(i) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
                             'MarkerFaceColor', hex2rgb(legend_colors_infer{i}), ...
                             'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end
legend(h_legend_infer, legend_entries_infer, 'FontSize', 24, ...
       'FontName', 'Times New Roman', 'Location', 'best', 'Interpreter', 'latex');

set(gca, 'FontSize', 30, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

% Save figure
saveas(fig3, 'power_comparison_annotated.png');
saveas(fig3, 'power_comparison_annotated.fig');
fprintf('  Saved: power_comparison_annotated.png/fig\n');

%% ========================================================================
%  HELPER FUNCTION: HEX TO RGB
%  ========================================================================

function rgb = hex2rgb(hex)
    % Convert hex color code to RGB triplet
    % Input: '#RRGGBB' or 'RRGGBB'
    % Output: [R G B] with values 0-1
    
    % Remove '#' if present
    hex = strrep(hex, '#', '');
    
    % Convert to RGB
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;
    
    rgb = [r g b];
end

fprintf('\n=== All figures created successfully! ===\n');
fprintf('Output files:\n');
fprintf('  - training_metrics_annotated.png/fig\n');
fprintf('  - inference_metrics_annotated.png/fig\n');
fprintf('  - power_comparison_annotated.png/fig\n');