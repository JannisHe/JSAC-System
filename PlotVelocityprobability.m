% This MATLAB code loads the position errors for a
% Joint Sensing and Communication (JSAC) system.
% Afterwards it plots the Probability of a correct
% Detection vs SNR values for a chosen ECA algorithm.
%
% Author: Jannis Held

% Load the data
ECA = load('radar_results_ECA.mat');  
ECA = ECA.velocity_errors;
ECA_B = load('radar_results_ECA-B.mat');  
ECA_B = ECA_B.velocity_errors;
ECA_C = load('radar_results_ECA-C.mat');  
ECA_C = ECA_C.velocity_errors;

SNR_values = [-60:2:20];

vel_error_threshold = 3* 3.6142;  % Set to 3x the resolution (3.6142 m/s)

prob_correct_vel_ECA = sum(ECA <= vel_error_threshold, 1) / length(ECA) * 100;
prob_correct_vel_ECA_B = sum(ECA_B <= vel_error_threshold, 1) / length(ECA) * 100;
prob_correct_vel_ECA_C = sum(ECA_C <= vel_error_threshold, 1) / length(ECA) * 100;

figure;
plot(SNR_values, prob_correct_vel_ECA, '-o', 'LineWidth', 1.5,'MarkerSize', 3);
hold on;
plot(SNR_values, prob_correct_vel_ECA_B, '-x', 'LineWidth', 1.5,'MarkerSize', 3);
plot(SNR_values, prob_correct_vel_ECA_C, '-s', 'LineWidth', 1.5,'MarkerSize', 3);
hold off;

%  plot
xlabel('SNR (dB)');
ylabel('Probability of Correct Velocity Detection (%)');
title('Probability of Correct Velocity Detection vs. SNR');
legend('ECA', 'ECA-B', 'ECA-C', 'Location', 'SouthEast');
grid on;


% Apply moving average smoothing
window_size = 3;  % Define the window size for moving average smoothing
prob_correct_vel_ECA_smooth = movmean(prob_correct_vel_ECA, window_size);
prob_correct_vel_ECA_B_smooth = movmean(prob_correct_vel_ECA_B, window_size);
prob_correct_vel_ECA_C_smooth = movmean(prob_correct_vel_ECA_C, window_size);

% Plot the smoothed probability of correct velocity detection vs. SNR
figure;
plot(SNR_values, prob_correct_vel_ECA_smooth, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
hold on;
plot(SNR_values, prob_correct_vel_ECA_B_smooth, '-x', 'LineWidth', 1.5, 'MarkerSize', 3);
plot(SNR_values, prob_correct_vel_ECA_C_smooth, '-s', 'LineWidth', 1.5, 'MarkerSize', 3);
hold off;

% Customize plot
xlabel('SNR (dB)');
ylabel('Probability of Correct Velocity Detection (%)');
title('Probability of Correct Velocity Detection vs. SNR');
legend('ECA', 'ECA-B', 'ECA-C', 'Location', 'SouthEast');
grid on;