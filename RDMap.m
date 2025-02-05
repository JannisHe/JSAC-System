% This MATLAB code simulates the generation of a range-Doppler map for a 
% Joint Sensing and Communication (JSAC) system. It includes options to add 
% clutter to the environment and allows the user to choose different clutter 
% suppression algorithms (ECA, ECA-B, ECA-B Parallel, or ECA-C)
%
% Author: Jannis Held

clear
rng('default');

addClutter = 1;                         % Set to true to add clutter
ECAAlgorithm = 'ECA-C';                 % Choose Algorithm: 'ECA', 'ECA-B','ECA-B Parallel', 'ECA-C'
SNR_dB = 10;                            % Desired SNR in dB
fc = 28e9;                              % Carrier frequency (Hz)
B = 100e6;                              % Bandwidth (Hz)
fs = B;

Pt_dBm = 23;
Pt = 10^-3*10^(Pt_dBm/10);               %Transmit Power
Gtx = 30;                               % Tx antenna gain (dB)

NodePos = [0; 0; 0];                     % Transceiver position
NodeVel = [0; 0; 0];                     % Transceiver velocity

NodePlatform = phased.Platform('InitialPosition', NodePos, 'Velocity', NodeVel);

Grx = 20;                               % Rx antenna gain (dB)
NF = 2.9;                               % Noise Figure (dB)
Tref = 290;                             % Reference temperature (K)

Rmax = 200;                             % Maximum range
vrelmax = 60;                           % Maximum relative velocity

TargetPos = [110 65 80; -10 5 -4; 0 0 0];  % Target positions
TargetRanges = sqrt(sum(TargetPos.^2, 1));  % Target Ranges
TargetVel = [-15 40 -32; 0 0 0; 0 0 0];    % Target velocities
TargetPlatform = phased.Platform('InitialPosition', TargetPos, 'Velocity', TargetVel);

TargetRCS = [1 0.5 0.7];  % RCS for three targets
Target = phased.RadarTarget('Model', 'Swerling1', 'MeanRCS', TargetRCS, 'OperatingFrequency', fc);

%Transmitter
Transmitter = phased.Transmitter('Gain', Gtx, 'PeakPower', Pt);
%Antenna
Ant = phased.IsotropicAntennaElement;
%Radiator
Radiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);

% Collector and Receiver
Collector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
Receiver = phased.ReceiverPreamp('SampleRate', fs, 'Gain', Grx, 'NoiseFigure', NF, 'ReferenceTemperature', Tref);

% Channel
Channel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);

%%
fdmax = 2 * speed2dop(vrelmax, physconst('Lightspeed')/fc);

Nsc = 1024;                              % Number of subcarriers
df = B/Nsc;                              % SCS
Tsym = 1/df;                             % OFDM symbol duration

Tcp = range2time(Rmax);                  % Duration of the cyclic prefix (CP)
Ncp = ceil(fs*Tcp);                      % Length of the CP in samples
Tcp = Ncp/fs;                           

Tofdm = Tsym + Tcp;                      % Total OFDM symbol duration with CP
Nofdm = Nsc + Ncp;                       % Number of samples in one OFDM symbol

nullIdx = [1:9 (Nsc/2+1) (Nsc-8:Nsc)]';  % Guard bands and DC subcarrier
Nscd = Nsc-length(nullIdx);              % Number of data subcarriers

bps = 6;                                 % Bits per Symbol
K = 2^bps;                               % Modulation order
Mofdm = 128;                             % Number of transmitted OFDM symbols
BitsTx = randi([0,1], [Nscd*bps Mofdm]);
SymbolsTx = qammod(BitsTx, K, 'InputType', 'bit', 'UnitAveragePower', true);
OFDM_Signal = ofdmmod(SymbolsTx, Nsc, Ncp, nullIdx);

Delta_f = 1/(Tofdm*Mofdm);

Tx_OFDM = reshape(OFDM_Signal, Nofdm, Mofdm);

%Normalize Power
Tx_OFDM = Tx_OFDM / rms(Tx_OFDM, 'all');

RX_OFDM = zeros(size(Tx_OFDM));

reset(NodePlatform);
reset(TargetPlatform);


% Clutter: Define clutter positions and Radar Cross Sections (RCS)
if addClutter
    ClutterPos = [50 30 0; 80 -20 0; 120 30 0; 20 -40 0]';  % Clutter positions
    ClutterRanges = sqrt(sum(ClutterPos.^2, 1));
    ClutterRCS = [50 30 100 80];  % RCS of CLutter
    Clutter = phased.RadarTarget('Model', 'Nonfluctuating', 'MeanRCS', ClutterRCS, 'OperatingFrequency', fc);
    ClutterRadiator = phased.Radiator('Sensor', Ant, 'OperatingFrequency', fc);
    ClutterChannel = phased.FreeSpace('SampleRate', fs, 'TwoWayPropagation', true, 'OperatingFrequency', fc);
    ClutterCollector = phased.Collector('Sensor', Ant, 'OperatingFrequency', fc);
end

% Transmission Loop: Transmitting and receiving signals for M OFDM symbols
for m = 1:Mofdm

    [NodePos, NodeVel] = NodePlatform(Tofdm);
    [TargetPos, TargetVel] = TargetPlatform(Tofdm);
    
    [TargetRange, TargetAngle] = rangeangle(TargetPos, NodePos);
    
    TxSignal = Transmitter(Tx_OFDM(:, m));
    
    radiatedSignal = Radiator(TxSignal, TargetAngle);
    
    ChannelSignal = Channel(radiatedSignal, NodePos, TargetPos, NodeVel, TargetVel);
    
    TargetReflection = Target(ChannelSignal, false);
    
    RxSignal = Collector(TargetReflection, TargetAngle);
    
    % Add clutter
    if addClutter
        [ClutterRange, ClutterAng] = rangeangle(ClutterPos, NodePos);
        RadiatedClutterSig = ClutterRadiator(TxSignal, ClutterAng);
        ChannelClutter = ClutterChannel(RadiatedClutterSig, NodePos, ClutterPos, NodeVel, zeros(size(ClutterPos)));
        ClutterSignal = Clutter(ChannelClutter);
        RxClutter = ClutterCollector(ClutterSignal, ClutterAng);
        RxSignal = RxSignal + RxClutter;
    end
    
        RX_OFDM(:, m) = Receiver(RxSignal);
end

% Adding Gaussian noise based on desired SNR (signal-to-noise ratio)
signal_power = rms(TargetReflection(:))^2;  % Signal power
noise_power = signal_power / (10^(SNR_dB/10));  % Noise Power based on SNR
noise = sqrt(noise_power) * (randn(size(RX_OFDM)) + 1j*randn(size(RX_OFDM))) / sqrt(2);  % Gaussian Noise
RxSig = RX_OFDM + noise;  

RxSigArray = reshape(RxSig, Nofdm*Mofdm, 1);

%%

%ECA Algorithm
[RxSigArray] = ClutterRemoval(ECAAlgorithm,Tx_OFDM,Nofdm,Mofdm,Delta_f,Tsym,RxSigArray,Nsc,Ncp);


%%
% Demodulate the OFDM Signal
RxDemod = ofdmdemod(RxSigArray, Nsc, Ncp, Ncp, nullIdx);

Zofdm = RxDemod./SymbolsTx;

% Range Doppler response
rdr = phased.RangeDopplerResponse('RangeMethod', 'FFT', 'SampleRate', fs, 'SweepSlope', -B/Tsym,...
    'DopplerOutput', 'Speed', 'OperatingFrequency', fc, 'PRFSource', 'Property', 'PRF', 1/Tofdm, ...
    'ReferenceRangeCentered', true);

%Plot range-Doppler Map
figure;
plotResponse(rdr, Zofdm, 'Unit', 'db');
xlim([-vrelmax vrelmax]);
ylim([0 Rmax]);
colormap(parula);  
colorbar; 

