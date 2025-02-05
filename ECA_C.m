function [s_ECA_C] = ECA_C(ref_signal, Nofdm, Mofdm, received_signal)
% ECA_C: Implements the ECA by Subcarrier (ECA-C) for clutter suppression
%
% Inputs:
% - ref_signal: Transmitted OFDM signal (reference signal)
% - Nofdm: Number of samples in one OFDM symbol
% - Mofdm: Number of OFDM symbols
% - received_signal: Surveillance signal (received signal with clutter)
%
% Output:
% - s_ECA_C: Cleaned surveillance signal after clutter suppression using ECA-C
%
% Author: Jannis Held

s_ref = reshape(ref_signal, Nofdm * Mofdm, 1);
s_surv = reshape(received_signal, Nofdm * Mofdm, 1); 

% Perform FFT to move both signals into the subcarrier domain
s_ref_SC = fft(reshape(s_ref, Nofdm, Mofdm), Nofdm, 1);  
s_surv_SC = fft(reshape(s_surv, Nofdm, Mofdm), Nofdm, 1); 

s_ECA_C = zeros(size(s_surv_SC));

% Process each subcarrier independently
for k = 1:Nofdm
    % Stack subcarrier signals across multiple OFDM symbols
    Ck = s_ref_SC(k, :).'; 
    Yk = s_surv_SC(k, :).'; 
    
    % Estimate clutter subspace and perform orthogonal projection
    lambda = 1e-6;  % Small regularization term
    alpha_k = (Ck' * Yk) / (Ck' * Ck + lambda);
    
    s_ECA_C(k, :) = (Yk - Ck * alpha_k).';
end

% Convert the cleaned subcarrier signals back to the time domain using IFFT
s_ECA_C = real(ifft(s_ECA_C, Nofdm, 1));

s_ECA_C = reshape(s_ECA_C, Nofdm * Mofdm, 1); 

end
