function [s_ECA] = ClutterRemoval(ECAAlgorithm,s_R,Nofdm,Mofdm,Delta_f,Tsym,s_surv,Nsc,Ncp)
% ClutterRemoval: Applies the selected ECA clutter removal algorithm to the input signal
% 
% Inputs:
% - ECAAlgorithm: String specifying the algorithm ('ECA', 'ECA-B', 'ECA-B Parallel', 'ECA-C')
% - s_R: Reference signal used for clutter removal
% - Nofdm: Number of OFDM samples in one symbol
% - Mofdm: Number of OFDM symbols
% - Delta_f: Subcarrier spacing
% - Tsym: OFDM symbol duration
% - s_surv: Received surveillance signal (with clutter)
% - Nsc: Number of subcarriers
% - Ncp: Cyclic prefix length
% 
% Output:
% - s_ECA: Output signal after clutter removal
%
% Author: Jannis Held
switch ECAAlgorithm
    case 'None'
        s_ECA = s_surv;
    case 'ECA'
        [s_ECA] = ECA(s_R,Nofdm,Mofdm,Delta_f,Tsym,s_surv);
    case 'ECA-B'
        batch_size = 37056; %Batch Size
        [s_ECA] = ECA_B(s_R, Nofdm, Mofdm, Delta_f, Tsym, s_surv, batch_size);
        
    case 'ECA-B Parallel'
        batch_size = 37056; %Batch Size
        [s_ECA] = ECA_BParallel(s_R, Nofdm, Mofdm, Delta_f, Tsym, s_surv, batch_size);
        
    case 'ECA-C'
        s_surv = reshape(s_surv, Nofdm, Mofdm);     %Reshape to Matrix
        s_surv_noCP = s_surv(Ncp+1:end, :);      % Remove cyclic prefix
        s_surv_noCP = reshape(s_surv_noCP, Nsc * Mofdm, 1);  % Reshape back to Vector
        
        [s_ECA_C] = ECA_C(s_R(Ncp+1:end,:), Nsc, Mofdm, s_surv_noCP);
        s_ECA_C_noCP = reshape(s_ECA_C, Nsc, Mofdm);
        s_ECA_C_CP = [s_surv(1:Ncp, :); s_ECA_C_noCP];  % Add back the original CP for Demodulation
        s_ECA = reshape(s_ECA_C_CP, Nofdm * Mofdm, 1);  % Reshape to match the original size
end

end

