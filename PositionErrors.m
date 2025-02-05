% This MATLAB code simulates the position errors for a
% Joint Sensing and Communication (JSAC) system. It includes options to choose 
% between different clutter suppression algorithms (ECA, ECA-B, ECA-B Parallel, or ECA-C).
% The results include position errors for varying SNR values and are saved for later analysis.
%
% Author: Jannis Held
clear;

parpool('local');

num_iterations = 2;  % Number of iterations for error calculation
addClutter = 1;  % Set to true to add clutter
num_clutter = 4;  % Number of static clutter objects
ECAAlgorithm = 'ECA';  % Choose Algorithm: 'ECA', 'ECA-B','ECA-B Parallel', 'ECA-C'
SNR_values = [-60:2:20];  % Array of SNR values to loop over
fc = 28e9;  % Carrier frequency (Hz)
B = 100e6;  % Bandwidth (Hz)
fs = B;  % Sampling frequency

% Transmitter parameters
Pt_dBm = 23;
Pt = 10^-3 * 10^(Pt_dBm / 10);  % Transmit power in Watts
Gtx = 30;  % Transmitter antenna gain (dB)

% Radar and system parameters
Rmax = 200;  % Maximum range of interest (m)
vrelmax = 60;  % Maximum relative velocity (m/s)
Grx = 20;  % Receiver antenna gain (dB)
NF = 2.9;  % Noise figure (dB)
Tref = 290;  % Reference temperature (Kelvin)
fdmax = 2 * speed2dop(vrelmax, physconst('Lightspeed') / fc);  % Maximum Doppler shift

% OFDM parameters
Nsc = 1024;  % Number of subcarriers
df = B / Nsc;  % Subcarrier spacing (Hz)
Tsym = 1 / df;  % OFDM symbol duration (seconds)

Tcp = range2time(Rmax);  % Calculate cyclic prefix (CP) duration based on max range
Ncp = ceil(fs * Tcp);  % Number of samples in the CP
Tofdm = Tsym + Tcp;  % Total OFDM symbol duration (with CP)
Nofdm = Nsc + Ncp;  % Total number of samples in one OFDM symbol

% OFDM subcarrier allocation (guard bands and data subcarriers)
nullIdx = [1:9 (Nsc/2+1) (Nsc-8:Nsc)]';  % Guard bands and DC subcarrier
Nscd = Nsc - length(nullIdx);  % Number of data subcarriers

bps = 6;  % Bits per QAM symbol
K = 2^bps;  % QAM modulation order
Mofdm = 128;  % Number of OFDM symbols transmitted

range_error_matrix = zeros(num_iterations, length(SNR_values));
velocity_error_matrix = zeros(num_iterations, length(SNR_values));

BitsTx = randi([0, 1], [Nscd * bps Mofdm]);  % Bits
SymbolsTx = qammod(BitsTx, K, 'InputType', 'bit', 'UnitAveragePower', true);  % QAM modulated data
OFDM_Signal = ofdmmod(SymbolsTx, Nsc, Ncp, nullIdx);  % OFDM modulated signal
Tx_OFDM = reshape(OFDM_Signal, Nofdm, Mofdm); 
Tx_OFDM = Tx_OFDM / rms(Tx_OFDM, 'all');  % Normalize OFDM transmission signal power

% Target and clutter setup
TargetRCS = 1;  % Radar Cross Section (RCS) of the target
ClutterPos = [50 30 0; 80 -20 0; 150 60 0; 20 -40 0]';  % Static clutter positions
ClutterRCS = [30 20 50 40];  % RCS values for the clutter objects

% Setup radar target model
Target = phased.RadarTarget('Model', 'Swerling1', 'MeanRCS', TargetRCS, 'OperatingFrequency', fc);

% Setup the transmitter, radiator, collector, and receiver
Transmitter = phased.Transmitter('Gain', Gtx, 'PeakPower', Pt);
Ant = phased.IsotropicAntennaElement;
Radiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);
Collector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
Receiver = phased.ReceiverPreamp('SampleRate', fs, 'Gain', Grx, 'NoiseFigure', NF, 'ReferenceTemperature', Tref);
Channel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);

% Setup clutter targets if clutter is enabled
if addClutter
    Clutter = phased.RadarTarget('Model', 'Nonfluctuating', 'MeanRCS', ClutterRCS, 'OperatingFrequency', fc);
    ClutterRadiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);
    ClutterChannel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);
    ClutterCollector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
end

true_positions = zeros(num_iterations, 3, length(SNR_values));
true_ranges = zeros(num_iterations, length(SNR_values));
true_velocities = zeros(num_iterations, 3, length(SNR_values));
estimated_ranges = zeros(num_iterations, length(SNR_values));
estimated_velocities = zeros(num_iterations, length(SNR_values));
velocity_errors = zeros(num_iterations, length(SNR_values));

% Main parallel loop to calculate errors for different SNR values
parfor snr_idx = 1:length(SNR_values)
    SNR_dB = SNR_values(snr_idx);  
    
    for iter = 1:num_iterations
        % Generate random target position and velocity for this iteration
        TargetPos = [120 * rand(1, 1); 25 * (2 * rand(1, 1) - 1); 0];  % Random target position
        TargetVel = [(2 * rand(1, 1) - 1) * (1 + (40 - 1) * rand); 0; 0];  % Random target velocity
        
        % Save true position and velocity
        true_positions(iter, :, snr_idx) = TargetPos';
        true_velocities(iter, :, snr_idx) = TargetVel';
        
        % Setup target platform
        TargetPlatform = phased.Platform('InitialPosition', TargetPos, 'Velocity', TargetVel);
        RX_OFDM = zeros(size(Tx_OFDM));  % Preallocate received OFDM signal
        
        % Transmit and receive OFDM symbols for each iteration
        for m = 1:Mofdm
            % Update target position and velocity
            [TargetPos, TargetVel] = TargetPlatform(Tofdm);
            [TargetRange, TargetAngle] = rangeangle(TargetPos, [0; 0; 0]);  % Range and angle of the target
            
            % Transmit OFDM signal
            TxSignal = Transmitter(Tx_OFDM(:, m));
            radiatedSignal = Radiator(TxSignal, TargetAngle);  % Radiate towards target
            
            % Propagate the signal through free space and receive the reflection
            ChannelSignal = Channel(radiatedSignal, [0; 0; 0], TargetPos, [0; 0; 0], TargetVel);
            TargetReflection = Target(ChannelSignal, false);
            RxSignal = Collector(TargetReflection, TargetAngle);  % Collect reflected signal
            
            % Add clutter reflections if clutter is enabled
            if addClutter
                [ClutterRange, ClutterAng] = rangeangle(ClutterPos, [0; 0; 0]);  % Calculate clutter angles
                RadiatedClutterSig = ClutterRadiator(TxSignal, ClutterAng);  % Radiate signal towards clutter
                ChannelClutter = ClutterChannel(RadiatedClutterSig, [0; 0; 0], ClutterPos, [0; 0; 0], zeros(size(ClutterPos)));
                ClutterSignal = Clutter(ChannelClutter);  % Get clutter reflections
                RxClutter = ClutterCollector(ClutterSignal, ClutterAng);  % Collect clutter signal
                RxSignal = RxSignal + RxClutter;  % Add clutter to received signal
            end
            
            RX_OFDM(:, m) = Receiver(RxSignal); 
        end
        
        % Add Gaussian noise to the received signal based on the SNR
        signal_power_target = rms(TargetReflection(:))^2;  % Calculate signal power
        noise_power = signal_power_target / (10^(SNR_dB / 10));  % Calculate noise power based on SNR
        noise = sqrt(noise_power) * (randn(size(RxSignal)) + 1j * randn(size(RxSignal))) / sqrt(2);  % Generate Gaussian noise
        RxSig = RX_OFDM + noise;  % Add noise to the received signal
        
        % Reshape received signal for processing
        RxSig = reshape(RxSig, Nofdm * Mofdm, 1);
        
        % Apply the selected clutter suppression algorithm
        RxSig = ClutterRemoval(ECAAlgorithm, Tx_OFDM, Nofdm, Mofdm, df, Tsym, RxSig, Nsc, Ncp);
        
        % Demodulate the received OFDM signal
        RxDemod = ofdmdemod(RxSig, Nsc, Ncp, Ncp, nullIdx);
        Zofdm = RxDemod ./ SymbolsTx;  % Calculate OFDM symbol ratios
        
        % Compute Range-Doppler response
        rdr = phased.RangeDopplerResponse('RangeMethod', 'FFT', 'SampleRate', fs, ...
            'SweepSlope', -B / Tsym, 'OperatingFrequency', fc, 'DopplerOutput', 'Speed', ...
            'PRFSource', 'Property', 'PRF', 1 / Tofdm);
        [resp, rng_grid, vel_grid] = rdr(Zofdm);  % Get response, range, and velocity grids
        
        % Find the peak response (estimate of range and velocity)
        [~, idx] = max(abs(resp(:)));
        [rng_idx, vel_idx] = ind2sub(size(resp), idx);
        rng_est = rng_grid(rng_idx);  % Estimated range
        vel_est = vel_grid(vel_idx);  % Estimated velocity
        
        % Save estimated range and calculate range error
        estimated_ranges(iter, snr_idx) = rng_est;
        true_range = norm(TargetPos);  % True target range
        range_error_matrix(iter, snr_idx) = abs(rng_est - true_range);
        
        % Save estimated velocity and calculate velocity error
        true_velocity = TargetVel(1);
        LOS = TargetPos / norm(TargetPos);  
        true_vel_Los = true_velocity * LOS(1);  % True velocity along line of sight
        estimated_velocities(iter, snr_idx) = vel_est;
        velocity_error = abs(vel_est - true_vel_Los);  % Velocity error
        
        % Store velocity error for analysis
        velocity_errors(iter, snr_idx) = velocity_error;
        
    end
end

% Save the results to a file
save(['radar_results_' ECAAlgorithm '.mat'], 'true_positions', 'true_ranges', 'true_velocities', ...
   'estimated_ranges', 'estimated_velocities', 'velocity_errors', 'SNR_values', ...
    'range_error_matrix', 'ECAAlgorithm');

% Shut down the parallel pool after processing
delete(gcp('nocreate'));  
toc
