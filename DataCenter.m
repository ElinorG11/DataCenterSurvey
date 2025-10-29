% Enhanced simulation code for AI Data Center Power Systems
% This version includes:
% 1) Phase-plane portraits
% 2) Parametric trade-off analysis
% 3) Quantitative stability metrics
% 4) Multi-parameter sensitivity analysis
% 5) Time-varying load from CSV file

clear; clc; close all;

%% ========================================================================
%  LOAD CSV DATA FOR TIME-VARYING LOAD
%  ========================================================================

% Load the CSV file
csv_data = readtable('training_metrics_with_states.csv');

% Extract and process time and power data
time_csv = csv_data.time_ms / 1000;              % Convert ms to seconds
time_csv = time_csv - time_csv(1);               % Start from t=0
PLstep_raw = csv_data.power_draw_w * 100 * 1e3;  % Scale power values

% Create input for Simulink (Nx2 matrix: [time, power])
PLstep = [time_csv, PLstep_raw];

% Verify data
fprintf('========================================\n');
fprintf('CSV DATA LOADED\n');
fprintf('========================================\n');
fprintf('Loaded %d data points\n', length(time_csv));
fprintf('Time range: %.3f to %.3f seconds\n', min(time_csv), max(time_csv));
fprintf('Power range: %.2f to %.2f MW\n', min(PLstep_raw)/1e6, max(PLstep_raw)/1e6);
fprintf('\n');

%% ========================================================================
%  PART 1: PARAMETER SWEEP SIMULATIONS
%  ========================================================================

% Base simulation parameters
Px = 2*max(PLstep_raw);      % Power transfer capability [W]
Prt = Px;                     % Rated generator power [W]
fs = 60;                      % Nominal frequency [Hz]
ws = 2 * pi * fs;             % Nominal frequency [rad/s]
K = 2.2e-04 * ws^2/Prt;       % Generator's inertia constant [1/(W*s^2)]
Pref = 0.5*mean(PLstep_raw);  % Generator's reference power [W]

% Simulink parameters
SimTime = max(time_csv) + 1;  % Simulation time [s]
RelTol = 1e-4;                % Simulation accuracy
MaxStep = 1e-3;               % Max step size [s]

% Define ranges for parameter sweeps
alpha_values = [0.01, 0.1, 1, 10, 100];  % Damping parameter sweep [1/s]
num_alpha = length(alpha_values);

% Storage for results
results = struct();

% Run simulations for different alpha values
disp('========================================');
disp('Running Parameter Sweep for Alpha');
disp('========================================');

for i = 1:num_alpha
    alpha = alpha_values(i);
    
    fprintf('Simulation %d/%d: alpha = %.2f [1/s]\n', i, num_alpha, alpha);
    
    % Run Simulink
    DataCenterSim;
    sim(bdroot);
    
    % Interpolate simulation outputs to match CSV time points
    % This ensures all vectors have the same length
    P_interp = interp1(ts, P, time_csv, 'linear', 'extrap');
    DeltaOmega_interp = interp1(ts, DeltaOmega, time_csv, 'linear', 'extrap');
    delta_interp = interp1(ts, delta, time_csv, 'linear', 'extrap');
    f_interp = (DeltaOmega_interp + ws)/(2*pi);
    
    % Store results - all vectors now have length = length(time_csv)
    results(i).alpha = alpha;
    results(i).ts = time_csv;              % CSV time vector
    results(i).PL = PLstep_raw;            % Load from CSV
    results(i).P = P_interp;               % Interpolated generator power
    results(i).Pg = PLstep_raw - P_interp; % Grid power (correct size!)
    results(i).DeltaOmega = DeltaOmega_interp;
    results(i).delta = delta_interp;
    results(i).f = f_interp;
    
    % Calculate quantitative metrics
    results(i).metrics = calculate_metrics(time_csv, P_interp, f_interp, ...
                                          delta_interp, Pref, fs, max(PLstep_raw));
end

disp('Parameter sweep completed.');
disp(' ');

%% ========================================================================
%  PART 2: ADDITIONAL PARAMETER SWEEPS (Grid Strength, Inertia, Pref)
%  ========================================================================

% 2.1) Grid Strength Variation (Px)
disp('========================================');
disp('Running Parameter Sweep for Grid Strength (Px)');
disp('========================================');

Px_values = [1.5*max(PLstep_raw), 2*max(PLstep_raw), 3*max(PLstep_raw), 4*max(PLstep_raw)];
num_Px = length(Px_values);
results_Px = struct();

alpha = 1;  % Fixed alpha for this sweep

for i = 1:num_Px
    Px = Px_values(i);
    
    fprintf('Simulation %d/%d: Px = %.0f MW\n', i, num_Px, Px/1e6);
    
    DataCenterSim;
    sim(bdroot);
    
    % Interpolate simulation outputs to match CSV time points
    P_interp = interp1(ts, P, time_csv, 'linear', 'extrap');
    DeltaOmega_interp = interp1(ts, DeltaOmega, time_csv, 'linear', 'extrap');
    delta_interp = interp1(ts, delta, time_csv, 'linear', 'extrap');
    f_interp = (DeltaOmega_interp + ws)/(2*pi);
    
    % Store results
    results_Px(i).Px = Px;
    results_Px(i).ts = time_csv;
    results_Px(i).PL = PLstep_raw;
    results_Px(i).P = P_interp;
    results_Px(i).Pg = PLstep_raw - P_interp;
    results_Px(i).f = f_interp;
    results_Px(i).delta = delta_interp;
    results_Px(i).DeltaOmega = DeltaOmega_interp;
    results_Px(i).metrics = calculate_metrics(time_csv, P_interp, f_interp, ...
                                             delta_interp, Pref, fs, max(PLstep_raw));
end

% Reset Px to base value
Px = 2*max(PLstep_raw);

% 2.2) Inertia Variation (K)
disp('========================================');
disp('Running Parameter Sweep for Inertia (K)');
disp('========================================');

K_multipliers = [0.5, 1, 2, 4];  % Multipliers for base K value
num_K = length(K_multipliers);
results_K = struct();

K_base = 2.2e-04 * ws^2/Prt;
alpha = 1;  % Fixed alpha for this sweep

for i = 1:num_K
    K = K_multipliers(i) * K_base;
    
    fprintf('Simulation %d/%d: K = %.2fx base\n', i, num_K, K_multipliers(i));
    
    DataCenterSim;
    sim(bdroot);
    
    % Interpolate simulation outputs to match CSV time points
    P_interp = interp1(ts, P, time_csv, 'linear', 'extrap');
    DeltaOmega_interp = interp1(ts, DeltaOmega, time_csv, 'linear', 'extrap');
    delta_interp = interp1(ts, delta, time_csv, 'linear', 'extrap');
    f_interp = (DeltaOmega_interp + ws)/(2*pi);
    
    % Store results
    results_K(i).K_mult = K_multipliers(i);
    results_K(i).ts = time_csv;
    results_K(i).PL = PLstep_raw;
    results_K(i).P = P_interp;
    results_K(i).Pg = PLstep_raw - P_interp;
    results_K(i).f = f_interp;
    results_K(i).delta = delta_interp;
    results_K(i).DeltaOmega = DeltaOmega_interp;
    results_K(i).metrics = calculate_metrics(time_csv, P_interp, f_interp, ...
                                            delta_interp, Pref, fs, max(PLstep_raw));
end

% Reset K to base value
K = K_base;

% 2.3) Reference Power Variation (Pref)
disp('========================================');
disp('Running Parameter Sweep for Pref');
disp('========================================');

Pref_fractions = [0.25, 0.5, 0.75, 1.0];  % Fractions of mean load
num_Pref = length(Pref_fractions);
results_Pref = struct();

alpha = 1;  % Fixed alpha for this sweep

for i = 1:num_Pref
    Pref = Pref_fractions(i) * mean(PLstep_raw);
    
    fprintf('Simulation %d/%d: Pref = %.2f * mean(PLstep)\n', i, num_Pref, ...
            Pref_fractions(i));
    
    DataCenterSim;
    sim(bdroot);
    
    % Interpolate simulation outputs to match CSV time points
    P_interp = interp1(ts, P, time_csv, 'linear', 'extrap');
    DeltaOmega_interp = interp1(ts, DeltaOmega, time_csv, 'linear', 'extrap');
    delta_interp = interp1(ts, delta, time_csv, 'linear', 'extrap');
    f_interp = (DeltaOmega_interp + ws)/(2*pi);
    
    % Store results
    results_Pref(i).Pref_frac = Pref_fractions(i);
    results_Pref(i).ts = time_csv;
    results_Pref(i).PL = PLstep_raw;
    results_Pref(i).P = P_interp;
    results_Pref(i).Pg = PLstep_raw - P_interp;
    results_Pref(i).f = f_interp;
    results_Pref(i).delta = delta_interp;
    results_Pref(i).DeltaOmega = DeltaOmega_interp;
    results_Pref(i).metrics = calculate_metrics(time_csv, P_interp, f_interp, ...
                                               delta_interp, Pref, fs, max(PLstep_raw));
end

% Reset Pref to base value
Pref = 0.5*mean(PLstep_raw);

disp('All parameter sweeps completed.');
disp(' ');

%% ========================================================================
%  PART 3: FIGURE GENERATION
%  ========================================================================

% Set default plotting parameters
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

colors = lines(max([num_alpha, num_Px, num_K, num_Pref]));

%% FIGURE 1: Original Time-Domain Response (Enhanced)
figure('Position', [100, 100, 1200, 900], 'Color', 'white');
label_fontsize = 30;

subplot(4,1,1);
for i = 1:num_alpha
    plot(results(i).ts, results(i).PL/1e6, 'k-', 'LineWidth', 2);
    hold on;
end
ylabel('$P_L$ [MW]', 'FontSize', label_fontsize);
% title('Load, Generator Power, Frequency, and Rotor Angle vs Time', ...
%       'FontSize', 18);
grid on;

subplot(4,1,2);
for i = 1:num_alpha
    plot(results(i).ts, results(i).P/1e6, 'LineWidth', 2, ...
         'Color', colors(i,:), 'DisplayName', sprintf('$\\alpha=%.2f$', ...
         results(i).alpha));
    hold on;
end
ylabel('$P$ [MW]', 'FontSize', label_fontsize);
legend('Location', 'best', 'FontSize', label_fontsize);
grid on;

subplot(4,1,3);
for i = 1:num_alpha
    plot(results(i).ts, results(i).f, 'LineWidth', 2, 'Color', colors(i,:));
    hold on;
end
yline(fs, 'k--', 'LineWidth', 1.5);
ylabel('$f$ [Hz]', 'FontSize', label_fontsize);
grid on;

subplot(4,1,4);
for i = 1:num_alpha
    plot(results(i).ts, results(i).delta * (180/pi), 'LineWidth', 2, ...
         'Color', colors(i,:));
    hold on;
end
ylabel('$\delta$ [deg]', 'FontSize', label_fontsize);
xlabel('Time [s]', 'FontSize', label_fontsize);
grid on;

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', label_fontsize);

%% FIGURE 2: Phase-Plane Portrait
figure('Position', [150, 150, 800, 700], 'Color', 'white');

fig2_label_fontsize = 30;

for i = 1:num_alpha
    % Convert to deviations from equilibrium
    % Find equilibrium values (after settling)
    delta_eq = results(i).delta(end);
    omega_eq = results(i).DeltaOmega(end);
    
    % Calculate deviations
    delta_dev = (results(i).delta - delta_eq) * (180/pi);  % Convert to degrees
    omega_dev = results(i).DeltaOmega - omega_eq;
    
    % Plot trajectory
    plot(delta_dev, omega_dev, 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', sprintf('$\\alpha=%.2f$ [1/s]', results(i).alpha));
    hold on;
    
    % Mark initial condition
    plot(delta_dev(1), omega_dev(1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
         'HandleVisibility', 'off');
end

% Mark equilibrium point
plot(0, 0, 'kx', 'MarkerSize', 15, 'LineWidth', 3, ...
     'DisplayName', 'Equilibrium');

xlabel('$\delta - \bar{\delta}$ [deg]', 'FontSize', fig2_label_fontsize);
ylabel('$\Delta\omega - \bar{\Delta\omega}$ [rad/s]', 'FontSize', fig2_label_fontsize);
%title('Phase-Plane Portrait', 'FontSize', 18);
legend('Location', 'best', 'FontSize', fig2_label_fontsize);
grid on;
axis equal;

% Set font sizes
set(gca, 'FontSize', fig2_label_fontsize);

%% FIGURE 3: Parametric Trade-off Analysis (Alpha)
figure('Position', [200, 200, 1200, 800], 'Color', 'white');

fig3_label_fontsize = 30;

% Extract metrics for plotting
alpha_plot = [results.alpha];
settling_times = arrayfun(@(x) x.metrics.settling_time, results);
peak_powers = arrayfun(@(x) x.metrics.peak_power/1e6, results);
peak_overshoots = arrayfun(@(x) x.metrics.peak_power_overshoot_pct, results);
freq_nadirs = arrayfun(@(x) x.metrics.frequency_nadir, results);
max_angles = arrayfun(@(x) x.metrics.max_rotor_angle, results);

subplot(2,3,1);
semilogx(alpha_plot, settling_times, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\alpha$ [1/s]', 'FontSize', fig3_label_fontsize);
ylabel('Settling Time [s]', 'FontSize', fig3_label_fontsize);
grid on;

subplot(2,3,2);
semilogx(alpha_plot, peak_powers, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\alpha$ [1/s]', 'FontSize', fig3_label_fontsize);
ylabel('Peak Power [MW]', 'FontSize', fig3_label_fontsize);
grid on;

subplot(2,3,3);
semilogx(alpha_plot, peak_overshoots, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\alpha$ [1/s]', 'FontSize', fig3_label_fontsize);
ylabel('Peak Overshoot [\%]', 'FontSize', fig3_label_fontsize);
grid on;

subplot(2,3,4);
semilogx(alpha_plot, freq_nadirs, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\alpha$ [1/s]', 'FontSize', fig3_label_fontsize);
ylabel('Frequency Nadir [Hz]', 'FontSize', fig3_label_fontsize);
yline(fs, 'k--', 'LineWidth', 1.5);
grid on;

subplot(2,3,5);
semilogx(alpha_plot, max_angles, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('$\alpha$ [1/s]', 'FontSize', fig3_label_fontsize);
ylabel('Max Rotor Angle [deg]', 'FontSize', fig3_label_fontsize);
grid on;

subplot(2,3,6);
% Trade-off plot: settling time vs peak overshoot
plot(peak_overshoots, settling_times, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Peak Overshoot [\%]', 'FontSize', fig3_label_fontsize);
ylabel('Settling Time [s]', 'FontSize', fig3_label_fontsize);
%title('Trade-off', 'FontSize', 16);
grid on;
% Add labels for each point
for i = 1:length(alpha_plot)
    text(peak_overshoots(i), settling_times(i), ...
         sprintf('  $\\alpha$=%.2f', alpha_plot(i)), ...
         'FontSize', 12);
end

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', fig3_label_fontsize);

%% FIGURE 4: Damping Comparison (Selected Alpha Values)
figure('Position', [250, 250, 1200, 400], 'Color', 'white');

fig4_label_fontsize = 30;

% Select representative alpha values for comparison
selected_idx = [1, 3, 5];  % Low, medium, high damping

subplot(1,3,1);
for idx = selected_idx
    plot(results(idx).ts, results(idx).P/1e6, 'LineWidth', 2, ...
         'Color', colors(idx,:), 'DisplayName', ...
         sprintf('$\\alpha=%.2f$', results(idx).alpha));
    hold on;
end
ylabel('$P$ [MW]', 'FontSize', fig4_label_fontsize);
xlabel('Time [s]', 'FontSize', fig4_label_fontsize);
%title('Generator Power', 'FontSize', 16);
legend('Location', 'best', 'FontSize', fig4_label_fontsize);
grid on;

subplot(1,3,2);
for idx = selected_idx
    plot(results(idx).ts, results(idx).f, 'LineWidth', 2, ...
         'Color', colors(idx,:));
    hold on;
end
yline(fs, 'k--', 'LineWidth', 1.5);
ylabel('$f$ [Hz]', 'FontSize', fig4_label_fontsize);
xlabel('Time [s]', 'FontSize', fig4_label_fontsize);
%title('Frequency', 'FontSize', 16);
grid on;

subplot(1,3,3);
for idx = selected_idx
    plot(results(idx).ts, results(idx).delta * (180/pi), 'LineWidth', 2, ...
         'Color', colors(idx,:));
    hold on;
end
ylabel('$\delta$ [deg]', 'FontSize', fig4_label_fontsize);
xlabel('Time [s]', 'FontSize', fig4_label_fontsize);
%title('Rotor Angle', 'FontSize', 16);
grid on;

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', fig4_label_fontsize);

%% FIGURE 5: Grid Strength (Px) Parameter Sweep
figure('Position', [300, 300, 1200, 800], 'Color', 'white');

fig5_label_fontsize = 30;

subplot(2,2,1);
for i = 1:num_Px
    plot(results_Px(i).ts, results_Px(i).P/1e6, 'LineWidth', 2, ...
         'Color', colors(i,:), 'DisplayName', ...
         sprintf('$P_x=%.0f$ MW', results_Px(i).Px/1e6));
    hold on;
end
ylabel('$P$ [MW]', 'FontSize', fig5_label_fontsize);
xlabel('Time [s]', 'FontSize', fig5_label_fontsize);
%title('Effect of Grid Strength', 'FontSize', 16);
legend('Location', 'best', 'FontSize', fig5_label_fontsize);
grid on;

subplot(2,2,2);
for i = 1:num_Px
    plot(results_Px(i).ts, results_Px(i).f, 'LineWidth', 2, ...
         'Color', colors(i,:));
    hold on;
end
yline(fs, 'k--', 'LineWidth', 1.5);
ylabel('$f$ [Hz]', 'FontSize', fig5_label_fontsize);
xlabel('Time [s]', 'FontSize', fig5_label_fontsize);
grid on;

subplot(2,2,3);
Px_plot = [results_Px.Px]/1e6;
settling_Px = arrayfun(@(x) x.metrics.settling_time, results_Px);
peak_Px = arrayfun(@(x) x.metrics.peak_power_overshoot_pct, results_Px);

yyaxis left
plot(Px_plot, settling_Px, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Settling Time [s]', 'FontSize', fig5_label_fontsize);
set(gca, 'YColor', 'b');

yyaxis right
plot(Px_plot, peak_Px, 's-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Overshoot [\%]', 'FontSize', fig5_label_fontsize);
xlabel('$P_x$ [MW]', 'FontSize', fig5_label_fontsize);
set(gca, 'YColor', 'r');
grid on;

subplot(2,2,4);
freq_nadir_Px = arrayfun(@(x) x.metrics.frequency_nadir, results_Px);
plot(Px_plot, freq_nadir_Px, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', 'm');
ylabel('Frequency Nadir [Hz]', 'FontSize', fig5_label_fontsize);
xlabel('$P_x$ [MW]', 'FontSize', fig5_label_fontsize);
yline(fs, 'k--', 'LineWidth', 1.5);
grid on;

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', fig5_label_fontsize);

%% FIGURE 6: Inertia (K) Parameter Sweep
figure('Position', [350, 350, 1200, 800], 'Color', 'white');

fig6_label_fontsize = 30;

subplot(2,2,1);
for i = 1:num_K
    plot(results_K(i).ts, results_K(i).P/1e6, 'LineWidth', 2, ...
         'Color', colors(i,:), 'DisplayName', ...
         sprintf('$K=%.1f \\times K_{base}$', results_K(i).K_mult));
    hold on;
end
ylabel('$P$ [MW]', 'FontSize', fig6_label_fontsize);
%title('Effect of Generator Inertia', 'FontSize', 18);
xlabel('Time [s]', 'FontSize', fig6_label_fontsize);
legend('Location', 'best', 'FontSize', fig6_label_fontsize);
grid on;

subplot(2,2,2);
for i = 1:num_K
    plot(results_K(i).ts, results_K(i).f, 'LineWidth', 2, ...
         'Color', colors(i,:));
    hold on;
end
yline(fs, 'k--', 'LineWidth', 1.5);
ylabel('$f$ [Hz]', 'FontSize', fig6_label_fontsize);
xlabel('Time [s]', 'FontSize', fig6_label_fontsize);
grid on;

subplot(2,2,3);
K_plot = [results_K.K_mult];
settling_K = arrayfun(@(x) x.metrics.settling_time, results_K);
peak_K = arrayfun(@(x) x.metrics.peak_power_overshoot_pct, results_K);

yyaxis left
plot(K_plot, settling_K, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Settling Time [s]', 'FontSize', 18);
set(gca, 'YColor', 'b');

yyaxis right
plot(K_plot, peak_K, 's-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Overshoot [\%]', 'FontSize', fig6_label_fontsize);
xlabel('$K$', 'FontSize', fig6_label_fontsize);
set(gca, 'YColor', 'r');
grid on;

subplot(2,2,4);
freq_nadir_K = arrayfun(@(x) x.metrics.frequency_nadir, results_K);
plot(K_plot, freq_nadir_K, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', 'm');
ylabel('Frequency Nadir [Hz]', 'FontSize', fig6_label_fontsize);
xlabel('$K$', 'FontSize', fig6_label_fontsize);
yline(fs, 'k--', 'LineWidth', 1.5);
grid on;

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', fig6_label_fontsize);

%% FIGURE 7: Reference Power (Pref) Parameter Sweep
figure('Position', [400, 400, 1200, 800], 'Color', 'white');

fig7_label_fontsize = 30;

subplot(2,2,1);
for i = 1:num_Pref
    plot(results_Pref(i).ts, results_Pref(i).P/1e6, 'LineWidth', 2, ...
         'Color', colors(i,:), 'DisplayName', ...
         sprintf('$P_{ref}=%.2f \\times \\overline{P_L}$', ...
         results_Pref(i).Pref_frac));
    hold on;
end
ylabel('$P$ [MW]', 'FontSize', fig7_label_fontsize);
%title('Effect of Generator Reference Power', 'FontSize', 16);
xlabel('Time [s]', 'FontSize', fig7_label_fontsize);
legend('Location', 'best', 'FontSize', fig7_label_fontsize);
grid on;

subplot(2,2,2);
for i = 1:num_Pref
    plot(results_Pref(i).ts, results_Pref(i).f, 'LineWidth', 2, ...
         'Color', colors(i,:));
    hold on;
end
yline(fs, 'k--', 'LineWidth', 1.5);
ylabel('$f$ [Hz]', 'FontSize', fig7_label_fontsize);
xlabel('Time [s]', 'FontSize', fig7_label_fontsize);
grid on;

subplot(2,2,3);
Pref_plot = [results_Pref.Pref_frac];
settling_Pref = arrayfun(@(x) x.metrics.settling_time, results_Pref);
peak_Pref = arrayfun(@(x) x.metrics.peak_power_overshoot_pct, results_Pref);

yyaxis left
plot(Pref_plot, settling_Pref, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Settling Time [s]', 'FontSize', fig7_label_fontsize);
set(gca, 'YColor', 'b');

yyaxis right
plot(Pref_plot, peak_Pref, 's-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Overshoot [\%]', 'FontSize', fig7_label_fontsize);
xlabel('$P_{ref}$ [\%]', 'FontSize', fig7_label_fontsize);
set(gca, 'YColor', 'r');
grid on;

subplot(2,2,4);
freq_nadir_Pref = arrayfun(@(x) x.metrics.frequency_nadir, results_Pref);
plot(Pref_plot, freq_nadir_Pref, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', 'm');
ylabel('Frequency Nadir [Hz]', 'FontSize', fig7_label_fontsize);
xlabel('$P_{ref}$ [\%]', 'FontSize', fig7_label_fontsize);
yline(fs, 'k--', 'LineWidth', 1.5);
grid on;

% Set font sizes for all subplots
allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', fig7_label_fontsize);

%% ========================================================================
%  PART 4: PRINT METRICS TABLE TO CONSOLE
%  ========================================================================

disp('========================================');
disp('QUANTITATIVE METRICS SUMMARY (Alpha Sweep)');
disp('========================================');
fprintf('\n');
fprintf('%10s | %15s | %15s | %15s | %15s | %15s\n', 'Alpha', ...
        'Settling Time', 'Peak Power', 'Overshoot %', 'Freq Nadir', 'Max Angle');
fprintf('%10s | %15s | %15s | %15s | %15s | %15s\n', '[1/s]', '[s]', '[MW]', ...
        '[%]', '[Hz]', '[deg]');
fprintf('%s\n', repmat('-', 1, 100));

for i = 1:num_alpha
    fprintf('%10.2f | %15.3f | %15.2f | %15.2f | %15.3f | %15.2f\n', ...
            results(i).alpha, ...
            results(i).metrics.settling_time, ...
            results(i).metrics.peak_power/1e6, ...
            results(i).metrics.peak_power_overshoot_pct, ...
            results(i).metrics.frequency_nadir, ...
            results(i).metrics.max_rotor_angle);
end

fprintf('\n');

%% ========================================================================
%  HELPER FUNCTION: Calculate Quantitative Metrics
%  ========================================================================

function metrics = calculate_metrics(ts, P, f, delta, Pref, fs, PLmax)
    % Calculate key performance metrics from simulation data
    
    % 1. Settling time (time to reach within 2% of final value for frequency)
    f_final = f(end);
    tolerance = 0.02 * abs(f_final - fs);  % 2% tolerance
    settled_indices = find(abs(f - f_final) <= tolerance);
    
    if ~isempty(settled_indices)
        % Find first time when it stays within tolerance
        settled_idx = settled_indices(1);
        for j = settled_idx:length(f)-1
            if abs(f(j) - f_final) > tolerance
                settled_idx = j+1;
            end
        end
        metrics.settling_time = ts(settled_idx);
    else
        metrics.settling_time = ts(end);  % Did not settle
    end
    
    % 2. Peak power output
    metrics.peak_power = max(P);
    
    % 3. Peak power overshoot (percentage above Pref)
    metrics.peak_power_overshoot_pct = ((metrics.peak_power - Pref) / Pref) * 100;
    
    % 4. Frequency nadir (minimum frequency)
    metrics.frequency_nadir = min(f);
    
    % 5. Maximum rotor angle deviation
    metrics.max_rotor_angle = max(abs(delta * 180/pi));
    
    % 6. Time to peak power
    [~, idx_peak] = max(P);
    metrics.time_to_peak = ts(idx_peak);
    
    % 7. Frequency deviation magnitude
    metrics.max_freq_deviation = abs(min(f) - fs);
    
    % 8. Number of oscillations (zero-crossings of frequency deviation)
    f_dev = f - f_final;
    zero_crossings = sum(diff(sign(f_dev)) ~= 0);
    metrics.num_oscillations = zero_crossings / 2;  % Full cycles
end
