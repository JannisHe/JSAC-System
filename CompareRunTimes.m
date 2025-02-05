% This MATLAB code simulates a number of iterations for a 
% Joint Sensing and Communication (JSAC) system. The runtimes will be saved
% for different clutter suppression algorithms (ECA, ECA-B, ECA-B Parallel, and ECA-C)
% and visualized in a normalized cumulative runtime plot.
%
% Author: Jannis Held
clear;
num_iterations = 100;  % Number of iterations
algorithms = {'ECA', 'ECA-B','ECA-B Parallel', 'ECA-C'};  % List of ECA algorithms

SNR_dB = 10;                              % Desired SNR in dB
fc = 28e9;                                % Carrier frequency (Hz)
B = 100e6;                                % Bandwidth (Hz)
fs = B;                                   % Sampling frequency
Nsc = 1024;                               % Number of subcarriers

Pt_dBm = 23;                              % Transmit power in dBm
Pt = 10^-3 * 10^(Pt_dBm/10);              % Transmit power in Watts
Gtx = 30;                                 % Transmitter antenna gain (dB)
Grx = 20;                                 % Receiver antenna gain (dB)
NF = 2.9;                                 % Noise figure (dB)
Tref = 290;                               % Reference temperature (Kelvin)

% Platform setup for node (transceiver)
NodePos = [0; 0; 0];                      % Node position in 3D space
NodeVel = [0; 0; 0];                      % Node velocity (stationary in this case)
NodePlatform = phased.Platform('InitialPosition', NodePos, 'Velocity', NodeVel);

df = B / Nsc;                             % Subcarrier spacing
Tsym = 1 / df;                            % OFDM symbol duration (without CP)

Rmax = 200;                               % Maximum target range (meters)
vrelmax = 60;                             % Maximum relative velocity (m/s)

% Target setup with positions and velocities
TargetPos = [110 65 80; -10 5 -4; 0 0 0]; % Positions of three targets
TargetRanges = sqrt(sum(TargetPos.^2, 1)); % Ranges for each target
TargetVel = [-15 40 -32; 0 0 0; 0 0 0];   % Target velocities (m/s)
TargetPlatfrom = phased.Platform('InitialPosition', TargetPos, 'Velocity', TargetVel);

% Set up the transmitter and radiator
Transmitter = phased.Transmitter('Gain', Gtx, 'PeakPower', Pt);
Ant = phased.IsotropicAntennaElement;     % Isotropic antenna
Radiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);

% Set up the collector and receiver
Collector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
Receiver = phased.ReceiverPreamp('SampleRate', fs, 'Gain', Grx, 'NoiseFigure', NF, 'ReferenceTemperature', Tref);

% Set up the radar channel for free space
Channel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);

% Target Radar Cross Section (RCS) values for three targets
TargetRCS = [1 0.5 0.7]; 
Target = phased.RadarTarget('Model', 'Swerling1', 'MeanRCS', TargetRCS, 'OperatingFrequency', fc);

% Calculate cyclic prefix (CP) length
Tcp = range2time(Rmax);                   % Duration of the cyclic prefix based on max range
Ncp = ceil(fs * Tcp);                     % Number of samples in the CP
Tcp = Ncp / fs;                           

% OFDM symbol and data setup
Tofdm = Tsym + Tcp;                       % Total OFDM symbol duration (with CP)
Nofdm = Nsc + Ncp;                        % Total number of samples in one OFDM symbol
nullIdx = [1:9 (Nsc/2+1) (Nsc-8:Nsc)]';   % Guard bands and DC subcarrier indices
Nscd = Nsc - length(nullIdx);             % Number of data subcarriers
bps = 6;                                  % Bits per QAM symbol
K = 2^bps;                                % QAM modulation order
Mofdm = 128;                              % Number of OFDM symbols

% Generate random data bits for transmission
BitsTx = randi([0,1], [Nscd * bps Mofdm]);
SymbolsTx = qammod(BitsTx, K, 'InputType', 'bit', 'UnitAveragePower', true);  % QAM modulation

% Generate OFDM signal
OFDM_Signal = ofdmmod(SymbolsTx, Nsc, Ncp, nullIdx); % OFDM modulation
Delta_f = 1 / (Tofdm * Mofdm);                      % Subcarrier spacing in Hz
Tx_OFDM = reshape(OFDM_Signal, Nofdm, Mofdm);       
Tx_OFDM = Tx_OFDM / rms(Tx_OFDM, 'all');            % Normalize transmit signal power
RX_OFDM = zeros(size(Tx_OFDM));                     

reset(NodePlatform);                                % Reset the transceiver platform
reset(TargetPlatfrom);                              % Reset the target platform

% Add clutter if enabled
addClutter = 1;
if addClutter
    % Define clutter positions and Radar Cross Sections (RCS)
    ClutterPos = [50 30 0; 80 -20 0; 120 30 0; 20 -40 0]';  % Clutter positions 
    ClutterRanges = sqrt(sum(ClutterPos.^2, 1));            % Clutter ranges
    ClutterRCS = [50 30 100 80];                            % Clutter RCS values
    
    % Set up the clutter radar targets and channels
    Clutter = phased.RadarTarget('Model', 'Nonfluctuating', 'MeanRCS', ClutterRCS, 'OperatingFrequency', fc);
    ClutterRadiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);
    ClutterChannel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);
    ClutterCollector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
end

% Main loop to transmit and receive OFDM signals
for m = 1:Mofdm
    % Update node and target positions for each OFDM symbol
    [NodePos, NodeVel] = NodePlatform(Tofdm);
    [TargetPos, TargetVel] = TargetPlatfrom(Tofdm);
    
    % Calculate the range and angle of the targets relative to the node
    [TargetRange, TargetAngle] = rangeangle(TargetPos, NodePos);
    
    % Transmit OFDM symbol
    TxSignal = Transmitter(Tx_OFDM(:, m));
    radiatedSignal = Radiator(TxSignal, TargetAngle);  % Radiate signal towards target
    
    % Propagate signal through free space channel and receive target reflection
    ChannelSignal = Channel(radiatedSignal, NodePos, TargetPos, NodeVel, TargetVel);
    TargetReflection = Target(ChannelSignal, false);   % Target reflection
    
    % Collect the reflected signal at the receiver
    RxSignal = Collector(TargetReflection, TargetAngle);
    
    % Add clutter reflections if enabled
    if addClutter
        [ClutterRange, ClutterAng] = rangeangle(ClutterPos, NodePos);
        RadiatedClutterSig = ClutterRadiator(TxSignal, ClutterAng);  % Radiate signal towards clutter
        ChannelClutter = ClutterChannel(RadiatedClutterSig, NodePos, ClutterPos, NodeVel, zeros(size(ClutterPos)));
        ClutterSignal = Clutter(ChannelClutter);                      % Clutter reflection
        RxClutter = ClutterCollector(ClutterSignal, ClutterAng);      % Collect clutter signal
        RxSignal = RxSignal + RxClutter;                             % Add clutter to received signal
    end
    
    RX_OFDM(:, m) = Receiver(RxSignal);  % Pass received signal through receiver
end

% Add Gaussian noise to the received signal based on the desired SNR
signal_power = rms(TargetReflection(:))^2;  % Calculate signal power
noise_power = signal_power / (10^(SNR_dB / 10));  % Calculate noise power based on SNR
noise = sqrt(noise_power) * (randn(size(RX_OFDM)) + 1j * randn(size(RX_OFDM))) / sqrt(2);  % Generate complex Gaussian noise
RxSig = RX_OFDM + noise;  % Add noise to the received OFDM signal

RxSigArray = reshape(RxSig, Nofdm * Mofdm, 1);  % Reshape received signal into a vector

timing_results = zeros(num_iterations, length(algorithms));

parpool('local');

% Run timing tests for each algorithm
for alg_idx = 1:length(algorithms)
    ECAAlgorithm = algorithms{alg_idx}; 
    fprintf('Running timing test for %s...\n', ECAAlgorithm);
    
    for iter = 1:num_iterations
        tic; 
        RxSigArray = ClutterRemoval(ECAAlgorithm, Tx_OFDM, Nofdm, Mofdm, df, Tsym, RxSigArray, Nsc, Ncp);
        timing_results(iter, alg_idx) = toc; 
    end
end

% Shut down the parallel pool after timing tests
delete(gcp('nocreate'));  

% Cumulative times and normalize to ECA runtime
cumulative_times = cumsum(timing_results); 
normalized_times = cumulative_times ./ max(cumulative_times(:, 1));  % Normalize using ECA runtime

% Plot results
figure;
plot(1:num_iterations, normalized_times(:, 1), 'r-', 'LineWidth', 2); 
hold on;
plot(1:num_iterations, normalized_times(:, 2), 'b--', 'LineWidth', 2); 
plot(1:num_iterations, normalized_times(:, 3), 'b-.', 'LineWidth', 2);  
plot(1:num_iterations, normalized_times(:, 4), 'g:', 'LineWidth', 2);   

hold off;

xlabel('Number of Iterations');
ylabel('Normalized Cumulative Runtime');
title('Normalized Runtime of ECA, ECA-B, ECA-B Parallel, and ECA-C');
legend('ECA', 'ECA-B','ECA-B Parallel', 'ECA-C', 'Location', 'NorthWest');

xlim([1 num_iterations]);
ylim([0 1]);
grid on;
