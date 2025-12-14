clc; clear; close all;
%% Requirmnet 1
% Parameters
T = 1; % Symbol duration in seconds
samples_per_symbol = 5;
Ts = 0.2; % sampling duartion in seconds
Num_bits = 10;
% Define normalized pulse
p = [5 4 3 2 1]/sqrt(55);
% Generate random binary bits (0s and 1s)
bits = randi([0 1], 1, Num_bits);
% Convert to +1/-1
signal = 2*bits - 1; % Converts 0->-1 and 1->+1
% Upsample the signal (insert 4 zeros between samples)
upsampled_signal = upsample(signal, samples_per_symbol);
% Convolve with pulse shaping function
y = conv(upsampled_signal, p, 'full');

figure;
stem(upsampled_signal, 'LineWidth', 1.5);
title("Upsmapling sequence every 200ms");
xlabel("time (ms)");
ylabel("Amp");
figure;
plot(y, 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1); % Add a dashed black line at y = 0
title("Convolution with pulse shape");
xlabel("time (ms)");
ylabel("Amp");
% Matched filter
h_matched = fliplr(p);
% Alternative Simple hold filter (all ones, normalized)
h_alt = ones(1,5)/sqrt(5); % Simple hold filter with 5 samples, normalized
% Filter the signal
y_matched = conv(y, h_matched, 'full');
y_alt = conv(y, h_alt, 'full');
% Correlator filter
y_reshaped = reshape(y(:, 1:Num_bits * samples_per_symbol), samples_per_symbol,
Num_bits); % size(y_reshaped) = [5, 10], 5 samples, 10 bits
y_correlator = zeros(samples_per_symbol * Num_bits, 1); % Output size: 5 * 10 = 50
% Loop over each bit
for n = 1:Num_bits
% Loop over each correlation output for the current bit
for i = 1:samples_per_symbol
sum_result = 0;
% Compute correlation: sum p(m) at each sample * y_reshaped(m, n) for each
sample within a bit in output y
for m = 1:i
sum_result = sum_result + p(m) * y_reshaped(m, n);
end
% Store result in y_correlator (dumb summation at T) for each bit to then
start over
y_correlator((n-1)*samples_per_symbol + i) = sum_result;
end
end
% Sampling points
sampling_indices = (1:Num_bits) * samples_per_symbol;
% Extract matched and alt samples at symbol instants
matched_samples = y_matched(sampling_indices);
alt_samples = y_alt(sampling_indices);
% Time vector adjusting
time_vector_matched = (0:length(y_matched)-1) * Ts + 0.2; % added 0.2 filter delay to
overcome the transition from 50ms domain to 10sec domain
time_vector_alt = (0:length(y_alt)-1) * Ts + 0.2;

time_vector_corr = (0:length(y_correlator)-1) * Ts + 0.2;
% Plot: Matched Filter vs Simple Hold Filter
figure;
subplot(2,1,1);
plot(time_vector_matched , y_matched, 'b-', 'LineWidth', 1.5); hold on;
stem(sampling_indices * Ts, matched_samples, 'bo', 'LineWidth', 1);
title('Matched Filter Output');
xlabel('Time (s)');
xlim([0,10]);
ylabel('Amplitude');
legend('Matched Filter','Samples');
grid on;
subplot(2,1,2);
plot(time_vector_alt, y_alt, 'r-', 'LineWidth', 1.5); hold on;
stem(sampling_indices * Ts, alt_samples, 'ro', 'LineWidth', 1);
title('Simple Hold Filter Output');
xlabel('Time (s)');
xlim([0,10]);
ylabel('Amplitude');
legend('Hold Filter','Samples');
grid on;
set(gcf, 'Position', [100 100 800 600]);
sgtitle('Filter Outputs Comparison');
% Plot: Matched Filter vs Correlator
figure;
plot(time_vector_matched,y_matched, 'm-', 'LineWidth', 1.5); hold on;
plot(time_vector_corr,y_correlator, 'Color', '[0 0.5 0.5]', 'LineWidth', 1.5);
stem(sampling_indices * Ts, matched_samples, 'r--o', 'LineWidth', 1);
title('Matched Filter vs Correlator Output');
xlabel('Time (s)');
xlim([0,10]);
ylabel('Amplitude');
legend('Matched Filter', 'Correlator', 'Samples');
grid on;
set(gcf, 'Position', [100 100 800 600]);
%% Requirment 2
% Parameters
Num_bits = 10000;
EbN0_dB = -2:1:5;
% Generate random binary bits (0s and 1s)
bits = randi([0 1], 1, Num_bits);
% Convert to +1/-1
signal = 2*bits - 1; % Converts 0->-1 and 1->+1
% Upsample the signal (insert 4 zeros between samples)
upsampled_signal = upsample(signal, samples_per_symbol);
% Convolve with pulse shaping function
y = conv(upsampled_signal, p, 'full');

% Adding Noise and Calculating BER
% Initialize BER arrays for matched filter, hold filters, and theoretical
BER_matched = zeros(size(EbN0_dB));
BER_hold = zeros(size(EbN0_dB));
BER_theo = zeros(size(EbN0_dB));
% Matched filter
h_matched = fliplr(p);
h_hold = 5 * ones(1,5); % Simple hold filter with 5 samples, un-normalized
h_hold = h_hold/(5*sqrt(5)); % Simple hold filter with 5 samples, normalized
% Sampling points
sampling_indices = (1:Num_bits) * samples_per_symbol;
% Loop over Eb/N0 values
for idx = 1:length(EbN0_dB)
EbN0 = 10^(EbN0_dB(idx)/10); % Convert Eb/N0 fr%om dB to linear scale
N0 = 1 / EbN0; % Since pulse p is normalized, Eb = 1
% Generate Gaussian noise with appropriate variance
noise = sqrt(N0/2) * randn(size(y)); % Generate AWGN noise
rx_signal = y + noise; % Received signal after adding noise
% Matched filter output
rx_matched = conv(rx_signal, h_matched, 'full');
matched_samples = rx_matched(sampling_indices);
% Hold filter output
rx_hold = conv(rx_signal, h_hold, 'full');
hold_samples = rx_hold(sampling_indices);
% Decision making (threshold = 0)
matched_decisions = matched_samples > 0;
hold_decisions = hold_samples > 0;
% Calculate BER for matched filter
BER_matched(idx) = sum(matched_decisions ~= bits) / Num_bits;
% Calculate BER for hold filter
BER_hold(idx) = sum(hold_decisions ~= bits) / Num_bits;
% Calculate Theoretical BER value
BER_theo(idx) = 0.5 * erfc(sqrt(EbN0));
end
% Plotting the BER vs SNR
figure;
semilogy(EbN0_dB, BER_matched, 'bo-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, BER_hold, 'g-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, BER_theo, 'k--', 'LineWidth', 1.5); % Theoretical BER (dashed line)
title('BER with Matched & Simple Hold Filters');
xlabel('SNR in dB');
ylabel('Bit Error Rate');
legend('BER using matched filter', 'BER using hold filter', 'Theoretical BER');
grid on;

%% Requirmnet 3
% ISI && Raised cosine
data_length = 100; % Length of the data (bits)
data = randi([0, 1], 1, data_length);
Data_sequence_maped = 2*data - 1; % maping the output to 1 or -1 %
Data_sequence_UpSample = upsample(Data_sequence_maped, samples_per_symbol);
i=1;
for R = [0, 1]
for delay = [2, 8]
h = rcosdesign(R, delay, samples_per_symbol ,"sqrt");
% Apply square root raised cosine first
A_Tx = filter(h ,1 , Data_sequence_UpSample);
% ideal channel free noise mean H(s)=1not make any change
B_Rx = filter(h ,1 , A_Tx); % Apply square root raised cosine first
%The outcome was applying raised cosine to the data.
% draw rsosine filter
figure(i);
i = i + 7;
plot(h);
title("Square Root Raised Cosine Filter (R: " + num2str(R) + ", Delay: "
+ num2str(delay) + ")");
xlabel("Time Samples");
ylabel("Amplitude");
ylim([min(h)-0.01 max(h)+0.01]);
% draw Eye diagram
eye_fig_TX=eyediagram(A_Tx(delay * samples_per_symbol + 1:end - delay *
samples_per_symbol)', 2 * samples_per_symbol);
set(eye_fig_TX,'Name',"A_Tx_eyediagram for R :" + R + " Delay: " +
delay);
eye_fig_RX=eyediagram(B_Rx(delay * samples_per_symbol +1:end - delay *
samples_per_symbol)', 2 * samples_per_symbol);
set(eye_fig_RX,'Name',"B_Rx_eyediagram for R :" + R + " Delay: " +
delay);
end
end