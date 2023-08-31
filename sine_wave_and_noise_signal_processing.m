% Generate a complex signal (sum of sine waves + noise)
fs = 1000; % Sampling frequency
t = 0:1/fs:5; % Time vector
frequencies = [50, 120, 250]; % Frequencies of the sine waves
amplitudes = [1, 0.5, 0.3]; % Amplitudes of the sine waves
phases = [0, 2/3*pi, 13/6*pi]; % Phase offsets for each sine wave
num_sines = length(frequencies);
signal = zeros(size(t));

for i = 1:num_sines
    signal = signal + amplitudes(i) * sin(2*pi*frequencies(i) * t + phases(i));
end

% Add noise
noise_amplitude = 0.5;
signal = signal + noise_amplitude * randn(size(t));

% Preprocessing: Remove DC offset
signal = signal - mean(signal);

% Design a low-pass Butterworth filter
cutoff_frequency = 150; % Cutoff frequency for the filter
order = 4; % Filter order
[b, a] = butter(order, cutoff_frequency/(0.5*fs), 'low');

% Apply the filter to the signal
filtered_signal = filtfilt(b, a, signal);

% Plot the original and filtered signals
figure;
subplot(2,1,1);
plot(t, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, filtered_signal);
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency analysis
N = length(filtered_signal);
frequency_spectrum = abs(fft(filtered_signal));
frequencies = linspace(0, fs/2, N/2 + 1);

% Calculate statistical features
mean_value = mean(filtered_signal);
std_deviation = std(filtered_signal);
skewness_value = skewness(filtered_signal);
kurtosis_value = kurtosis(filtered_signal);

% Analyze time-domain characteristics
[~, peak_index] = max(filtered_signal);
rise_time = t(peak_index) - t(find(filtered_signal > 0.1*max(filtered_signal), 1, 'first'));
decay_time = t(find(filtered_signal < 0.1*max(filtered_signal), 1, 'last')) - t(peak_index);

% Display results
disp('Statistical Features:');
disp(['Mean: ', num2str(mean_value)]);
disp(['Standard Deviation: ', num2str(std_deviation)]);
disp(['Skewness: ', num2str(skewness_value)]);
disp(['Kurtosis: ', num2str(kurtosis_value)]);
disp('Time-Domain Characteristics:');
disp(['Rise Time: ', num2str(rise_time)]);
disp(['Decay Time: ', num2str(decay_time)]);

% Plot frequency spectrum
figure;
plot(frequencies, 2/N * frequency_spectrum(1:N/2 + 1));
title('Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Spectrogram
figure;
spectrogram(filtered_signal, 256, 250, [], fs, 'yaxis');
title('Spectrogram');
