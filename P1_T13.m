% Initialize Parameters
num_waveforms = 500;        
Amplitude = 4;               
samples_per_bit = 7; 
half_bit_duration = 4;
Samples_length = 700;
num_bits_with_delay = 101; 

% Matrix to store waveforms
ensemble_unipolar_nrz = zeros(num_waveforms, Samples_length);
ensemble_polar_nrz = zeros(num_waveforms, Samples_length); 
ensemble_polar_rz = zeros(num_waveforms, Samples_length);

% Generate Waveforms Line code and delay
for i = 1:num_waveforms
    data = randi([0 1], 1, num_bits_with_delay); % choose random number 0 or 1
    
    % Unipolar NRZ Signal
    unipolar_signal_nrz = data * Amplitude;
    unipolar_nrz_waveform = repmat(unipolar_signal_nrz, samples_per_bit, 1);
    unipolar_nrz_waveform = reshape(unipolar_nrz_waveform, 1, []);

    % polar NRZ Signal
    polar_signal = ((2 * data) - 1) * Amplitude;
    polar_nrz_waveform = repmat(polar_signal, samples_per_bit, 1);
    polar_nrz_waveform = reshape(polar_nrz_waveform, 1, []);

    % polar RZ Signal
    
    polar_rz_waveform = zeros(samples_per_bit, num_bits_with_delay); % Pre-allocation to avoid warning and set zeros for last 3 samples
    polar_rz_waveform(1:half_bit_duration, :) = repmat(polar_signal, half_bit_duration, 1); % first 4 samples are polar signal , last 3 samples returns to zero 
    polar_rz_waveform = reshape(polar_rz_waveform, 1, []);

    random_delay = randi([1 7]); % choose random number from 1 to 7
    % Shift waveforms to add random delay
    delayed_unipolar_nrz = unipolar_nrz_waveform(1 + random_delay : end); 
    delayed_polar_nrz = polar_nrz_waveform(1 + random_delay : end); 
    delayed_polar_rz = polar_rz_waveform(1 + random_delay : end);

    % Trim excess samples more than 700
    if length(delayed_unipolar_nrz) > Samples_length
       delayed_unipolar_nrz = delayed_unipolar_nrz(1:Samples_length);  
    end

    if length(delayed_polar_nrz) > Samples_length
       delayed_polar_nrz = delayed_polar_nrz(1:Samples_length);  
    end

    if length(delayed_polar_rz) > Samples_length
       delayed_polar_rz = delayed_polar_rz(1:Samples_length);  
    end

    ensemble_unipolar_nrz(i, :) = delayed_unipolar_nrz;
    ensemble_polar_nrz(i, :) = delayed_polar_nrz;
    ensemble_polar_rz(i, :) = delayed_polar_rz;
end
% User input to choose waveform type (unchanged)
disp('Choose a linecode type to view:');
disp('1. Unipolar NRZ');
disp('2. Polar NRZ');
disp('3. Polar RZ');
choice = input('Enter the number corresponding to your choice (1-3): ');

switch choice
    case 1
        plotWaveforms(ensemble_unipolar_nrz, 'Unipolar NRZ', Amplitude, samples_per_bit);
        plotStats(ensemble_unipolar_nrz, 'Unipolar NRZ', Amplitude, samples_per_bit, num_waveforms, Samples_length);
    case 2
        plotWaveforms(ensemble_polar_nrz, 'Polar NRZ', Amplitude, samples_per_bit);
        plotStats(ensemble_polar_nrz, 'Polar NRZ', Amplitude, samples_per_bit, num_waveforms, Samples_length);
    case 3
        plotWaveforms(ensemble_polar_rz, 'Polar RZ', Amplitude, samples_per_bit);
        plotStats(ensemble_polar_rz, 'Polar RZ', Amplitude, samples_per_bit, num_waveforms, Samples_length);
    otherwise
        disp('Invalid choice! Please run again and enter a number between 1 and 3.');
end
% Function to plot waveforms
function plotWaveforms(ensemble, type, Amplitude, samples_per_bit)
    view_time = 10 * samples_per_bit + 1; % view excatly 70 samples can extend up to 700 
    num_plots = 3;  

    figure;
    for i = 1:num_plots
        subplot(num_plots, 1, i);
        stairs(0:view_time-1, ensemble(i, 1:view_time), 'LineWidth', 1.5);
        xlabel('Sample Index');  
        ylabel('Voltage');  
        title([type, ' Waveform ', num2str(i)]);  
        if strcmp(type, 'Unipolar NRZ') % comapre type of line code chose by user to "unipolar" to decided amplitude  
            ylim([-1, Amplitude + 1]);  
        else
            ylim([-Amplitude - 1, Amplitude + 1]);  
        end
        grid on;
    end
end

% Function to calculate statistical mean (ensemble mean at each sample point)
function stat_mean = calc_statistical_mean(ensemble, num_waveforms)
    stat_mean=sum(ensemble, 1)/num_waveforms;
end

% Function to calculate time mean for one and all waveforms
function time_means = calc_time_means(ensemble, Samples_length)
  % for 1 waveform
  first_waveform = ensemble(1, :);
  time_mean = sum(first_waveform) / Samples_length;
  disp("The time mean of the 1st first waveform " + time_mean);

  %for average 
  time_means = sum(ensemble, 2)/Samples_length;
end

% Function to calculate ensemble autocorrelation Rx(tau)
function Rx = calc_ensemble_autocorr(ensemble, num_waveforms, Samples_length, max_tau)
    Rx_pos = zeros(1, max_tau + 1);
    for tau = 0:max_tau
        sum_val = 0;
        for t = 1:(Samples_length - tau)
            for i = 1:num_waveforms    % we can remove this loop and use sum_val = sum_val + sum(ensemble(:, t) .* ensemble(:, t+tau));
                sum_val = sum_val + ensemble(i, t) * ensemble(i, t + tau); 
            end
        end
        Rx_pos(tau + 1) = sum_val / (num_waveforms * (Samples_length - tau));
    end
    Rx = [flip(Rx_pos(2:end)), Rx_pos]; % Mirror for negative lags
end

% Function to calculate statistical autocorrelation at a specific time t1
function Rx_stat_t = calc_stat_autocorr_at_t(ensemble, num_waveforms, Samples_length, max_tau, t1)
    Rx_stat_t_pos = zeros(1, max_tau + 1);
    for tau = 0:max_tau
        if t1 + tau <= Samples_length  % Check bounds for positive tau
            sum_val = 0;
            for i = 1:num_waveforms
                sum_val = sum_val + ensemble(i, t1) * ensemble(i, t1 + tau);
            end
            Rx_stat_t_pos(tau + 1) = sum_val / num_waveforms;
        else
            Rx_stat_t_pos(tau + 1) = 0; % Out of bounds
        end
    end
    Rx_stat_t = [flip(Rx_stat_t_pos(2:end)), Rx_stat_t_pos]; % Mirror for negative lags
end

% Function to calculate time autocorrelation Rx(tau) for one waveform
function Rx_time = calc_time_autocorr(signal, sample_length,max_tau)
    shifted_matrix = zeros(sample_length, sample_length);
    % Create circularly shifted versions of the input signal
    for i = 1:sample_length
        shifted_matrix(i, :) = [signal(i:end), signal(1:i-1)];
    end
    % Compute time-domain autocorrelation
    autocorr_positive = zeros(1, max_tau + 1);  % Initialize positive-lag autocorrelation vector
    for t = 1:sample_length
        autocorr_positive(t) = sum(shifted_matrix(1,:) .* shifted_matrix(t,:)) / sample_length;
    end

    % Mirror the autocorrelation to include negative lags (symmetric autocorrelation)
    Rx_time = [flip(autocorr_positive(2:end))'; autocorr_positive'];
end

% Function to plot all statistics
function plotStats(ensemble, type, Amplitude, samples_per_bit, num_waveforms, Samples_length)
    max_tau = 699;
    tau = -max_tau:max_tau;

    % Statistical mean
    stat_mean = calc_statistical_mean(ensemble, num_waveforms);
    
    % Time means for all waveforms
    time_means = calc_time_means(ensemble, Samples_length);

    % Ensemble autocorrelation
    Rx = calc_ensemble_autocorr(ensemble, num_waveforms, Samples_length, max_tau);

    % Time autocorrelation (single waveform)
    Rx_time = calc_time_autocorr(ensemble(1, :), Samples_length, max_tau);

    % Statistical autocorrelation at t1 = 1 and t1 = 8 (MATLAB indices start at 1)
    Rx_stat_t1 = calc_stat_autocorr_at_t(ensemble, num_waveforms, Samples_length, max_tau, 1);
    Rx_stat_t8 = calc_stat_autocorr_at_t(ensemble, num_waveforms, Samples_length, max_tau, 8);

    % Theoretical values
    if strcmp(type, 'Unipolar NRZ')
        theo_mean = Amplitude / 2;
        theo_Rx = (Amplitude^2 / 4) * (1 - abs(tau) / samples_per_bit) .* (abs(tau) <= samples_per_bit) + (Amplitude/2)^2;
    elseif strcmp(type, 'Polar NRZ')
        theo_mean = 0;
        theo_Rx = Amplitude^2 * (1 - abs(tau) / samples_per_bit) .* (abs(tau) <= samples_per_bit);
    else % Polar RZ
        theo_mean = 0;
        theo_Rx = (Amplitude^2 * 4 / samples_per_bit) * (1 - abs(tau) / 4) .* (abs(tau) <= 4);
    end
    % PSD from ensemble autocorrelation
    Ts = 0.01;
    Tb = samples_per_bit * Ts; % 0.07 s
    fs = 1 / Ts; % 100 Hz
    N = length(Rx); % 1399 with max_tau = 699
    f = (-N/2:N/2-1) * (fs / N); % Frequency axis in Hz
    delta_f = zeros(size(f));
    [~,idx0]= min(abs(f));        % Find index where f ≈ 0
    psd = abs(fftshift(fft(Rx))) / fs; % No padding, use full Rx
    if strcmp(type, 'Unipolar NRZ')
        tho_psd = (Amplitude^2 / 4) * Tb * (sinc(f * Tb)).^2; % Main term
        delta_f(idx0) = (Amplitude^2 )/4; 
        tho_psd = tho_psd + delta_f;
    elseif strcmp(type, 'Polar NRZ')
        tho_psd = Amplitude^2 * Tb * (sinc(f * Tb)).^2;
        tho_psd = tho_psd + delta_f;
    else % Polar RZ
        Tp = 4 * Ts; % 0.04 s, pulse width
        tho_psd = (Amplitude^2 *4/7) * 4/7 *Tb * (sinc(f * Tp)).^2; % Duty cycle adjusted
        tho_psd = tho_psd + delta_f;
    end

    % Plotting
    figure('Color', 'w');
    plot(1:Samples_length, stat_mean, 'b', 'LineWidth', 1.5);
    hold on;
    plot(1:Samples_length, theo_mean * ones(1, Samples_length), 'r', 'LineWidth', 1.5);
    xlabel('Sample Index'); ylabel('Statistical Mean'); 
    ylim([-Amplitude - 1, Amplitude + 1]);
    title([type, ' Statistical Mean (Ensemble)']);
    legend('Simulated', 'Theoretical'); 
    grid on;

    figure('Color', 'w');
    plot(1:num_waveforms, time_means, 'b', 'LineWidth', 1.5);
    hold on;
    plot(1:num_waveforms, theo_mean * ones(1, num_waveforms), 'r', 'LineWidth', 1.5);
    xlabel('Waveform Index'); ylabel('Time Mean');
    ylim([-Amplitude - 1, Amplitude + 1]);
    title([type, ' Time Mean Across All Waveforms']);
    legend('Simulated', 'Theoretical'); 
    grid on;

    figure('Color', 'w');
    plot(tau, Rx, 'b', 'LineWidth', 1.5);
    hold on;
    plot(tau, theo_Rx, 'r--', 'LineWidth', 1.5);
    xlabel('Lag (\tau)'); ylabel('Rx(\tau)'); 
    title([type, ' Ensemble Autocorrelation']);
    legend('Simulated', 'Theoretical'); 
    grid on;
    xlim([-700, 700]);

    figure('Color', 'w');
    plot(tau, Rx_time, 'b', 'LineWidth', 1.5);
    hold on;
    plot(tau, theo_Rx, 'r--', 'LineWidth', 1.5);
    xlabel('Lag (\tau)'); ylabel('Rx(\tau)'); 
    title([type, ' Time Autocorrelation (Waveform 1)']);
    legend('Simulated', 'Theoretical'); grid on;
    xlim([-700, 700]);

    figure('Color', 'w');
    plot(tau, Rx_stat_t1, 'b', 'LineWidth', 1.5);
    hold on;
    plot(tau, Rx_stat_t8, 'r--', 'LineWidth', 1.5);
    xlabel('Lag (\tau)'); ylabel('Rx(\tau)'); 
    title([type, ' Statistical Autocorrelation at t1=1 and t1=8']);
    legend('t1=1', 't1=8'); grid on;
    xlim([-100, 100]);

    % PSD Plot
    figure('Color', 'w');
    plot(f, psd, 'b', 'LineWidth', 1);
    hold on;
    plot(f, tho_psd, 'r--', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density');
    legend('Simulated', 'Theoretical');
    title([type, ' Power Spectral Density']);
    xlim([-50, 50]);
    ylim([-0.5,1.5])
    grid on;
end
