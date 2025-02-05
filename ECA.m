function [s_ECA] = ECA(xofdm,Nofdm,Mofdm,Delta_f,Tsym,s_surv)
% ECA: Extensive Cancellation Algorithm for clutter suppression
%
% Inputs:
% - xofdm: Reference OFDM signal
% - Nofdm: Number of OFDM samples in one symbol
% - Mofdm: Number of OFDM symbols
% - Delta_f: Subcarrier spacing
% - Tsym: OFDM symbol duration
% - s_surv: Received surveillance signal (with clutter)
%
% Outputs:
% - s_ECA: Cleaned surveillance signal after clutter suppression
%
% Author: Jannis Held

s_ref = reshape(xofdm, Nofdm*Mofdm, 1);
N = length(s_ref);
K = 200;  % Number of range bins (delays)
P = 0;   % Number of Doppler bins (for Doppler shifts)

M = (2*P + 1) * K;  % Total number of columns in X: K range bins, 2P + 1 Doppler shifts
X = zeros(N, M);   
  

col = 1;
for p = -P:P
    %Doppler Shift
    doppler_shift = exp(1j * 2 * pi * p * Delta_f * (0:N-1)' * Tsym);
    
    for k = 0:K-1
        %Range Delay
        delayed_signal = [zeros(k, 1); s_ref(1:N-k)];       
        shifted_signal = delayed_signal .* doppler_shift;   %Apply Doppler Shift     
        X(:, col) = shifted_signal;        
        col = col + 1;
    end
end

alpha = inv(X'*X)*X'*s_surv;
s_ECA = s_surv - X*alpha;   %Substract Clutter

end

