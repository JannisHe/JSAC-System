% This MATLAB code loads the velocity errors for a
% Joint Sensing and Communication (JSAC) system.
% Afterwards it plots the Probability of a correct
% Detection vs SNR values for a chosen ECA algorithm.
%
% Author: Jannis Held

% Load the data
ECA = load('radar_results_ECA.mat');  
ECA = ECA.range_error_matrix;
ECA_B = load('radar_results_ECA-B.mat');  
ECA_B = ECA_B.range_error_matrix;
ECA_C = load('radar_results_ECA-C.mat');  
ECA_C = ECA_C.range_error_matrix;

SNR_values = [-60:2:20];

range_error_threshold = 4.5;  % Set to 3x the resolution (1.5 meters)

prob_correct_detection_ECA = sum(ECA <= range_error_threshold, 1) / length(ECA) * 100;
prob_correct_detection_ECA_B = sum(ECA_B <= range_error_threshold, 1) / length(ECA) * 100;
prob_correct_detection_ECA_C = sum(ECA_C <= range_error_threshold, 1) / length(ECA) * 100;

figure;
plot(SNR_values, prob_correct_detection_ECA, '-o', 'LineWidth', 1.5,'MarkerSize', 3);
hold on;
plot(SNR_values, prob_correct_detection_ECA_B, '-x', 'LineWidth', 1.5,'MarkerSize', 3);
plot(SNR_values, prob_correct_detection_ECA_C, '-s', 'LineWidth', 1.5,'MarkerSize', 3);
hold off;

xlabel('SNR (dB)');
ylabel('Probability of Correct Detection (%)');
title('Probability of Correct Detection vs. SNR');
legend('ECA', 'ECA-B', 'ECA-C', 'Location', 'SouthEast');
grid on;


% Define the window size for moving average
window_size = 3;  % You can adjust this to control the level of smoothing

% Apply the moving average filter to smooth the detection probabilities
prob_correct_detection_ECA_smooth = movmean(prob_correct_detection_ECA, window_size);
prob_correct_detection_ECA_B_smooth = movmean(prob_correct_detection_ECA_B, window_size);
prob_correct_detection_ECA_C_smooth = movmean(prob_correct_detection_ECA_C, window_size);

% Plot the smoothed probabilities of correct detection vs. SNR for each algorithm
figure;
plot(SNR_values, prob_correct_detection_ECA_smooth, '-o', 'LineWidth', 1.5, 'MarkerSize', 3);
hold on;
plot(SNR_values, prob_correct_detection_ECA_B_smooth, '-x', 'LineWidth', 1.5, 'MarkerSize', 3);
plot(SNR_values, prob_correct_detection_ECA_C_smooth, '-s', 'LineWidth', 1.5, 'MarkerSize', 3);
hold off;

% Customize plot
xlabel('SNR (dB)');
ylabel('Probability of Correct Detection (%)');
title('Probability of Correct Detection vs. SNR');
legend('ECA', 'ECA-B', 'ECA-C', 'Location', 'SouthEast');
grid on;


