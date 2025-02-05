function [s_ECA_B] = ECA_B(xofdm, Nofdm, Mofdm, Delta_f, Tsym, s_surv, batch_size)
% ECA_B: Implements the ECA-Batches (ECA-B) algorithm for clutter suppression
%
% Inputs:
% - xofdm: Transmitted OFDM signal (reference signal)
% - Nofdm: Number of samples in one OFDM symbol
% - Mofdm: Number of OFDM symbols
% - Delta_f: Subcarrier spacing
% - Tsym: OFDM symbol duration
% - s_surv: Surveillance signal (received signal with clutter)
% - batch_size: Number of samples in each batch for the batch-wise processing
%
% Output:
% - s_ECA_B: Output signal after clutter suppression using the ECA-B algorithm
%
% Author: Jannis Held

s_ref = reshape(xofdm, Nofdm * Mofdm, 1);
N = length(s_ref);  % Total number of samples in the reference signal
K = 200;            % Number of range bins (delays)
P = 0;              % Number of Doppler bins (for Doppler shifts)

% Number of batches to process
num_batches = ceil(N / batch_size);

s_ECA_B = zeros(size(s_surv));

% Process each batch separately
for batch_idx = 1:num_batches
    % Define the start and end of the current batch
    batch_start = (batch_idx - 1) * batch_size + 1;
    batch_end = min(batch_idx * batch_size, N);
    
    % Extract the current batch of the surveillance and reference signals
    s_surv_batch = s_surv(batch_start:batch_end);
    s_ref_batch = s_ref(batch_start:batch_end);
    N_batch = length(s_ref_batch);  % Number of samples in the current batch
    
    % Build matrix X for the current batch
    M = (2 * P + 1) * K;  % Number of columns in X (delayed and Doppler shifted versions)
    X_batch = zeros(N_batch, M);
    
    col = 1;
    for p = -P:P
        doppler_shift = exp(1j * 2 * pi * p * Delta_f * (0:N_batch-1)' * Tsym);  % Doppler Shift
        for k = 0:K-1
            if k < N_batch
                delayed_signal = [zeros(k, 1); s_ref_batch(1:N_batch-k)];  % Apply delay
            else
                delayed_signal = zeros(N_batch, 1);  % If delay exceeds batch size, zero-fill
            end
            shifted_signal = delayed_signal .* doppler_shift;  % Apply Doppler shift
            X_batch(:, col) = shifted_signal;
            col = col + 1;
        end
    end
    
    % Apply ECA for the current batch
    
    lambda = 1e-2;  % Small regularization parameter
    alpha_batch = inv(X_batch' * X_batch + lambda * eye(size(X_batch, 2))) * X_batch' * s_surv_batch;
    s_ECA_batch = s_surv_batch - X_batch * alpha_batch;
    
    % Combine the processed batch into the final signal
    s_ECA_B(batch_start:batch_end) = s_ECA_batch;
end

end
