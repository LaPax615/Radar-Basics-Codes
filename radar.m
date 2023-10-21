close all;
clear;
clc;

B = 200e6;  % Frequency range
T = 0.4e-3;  % Chirp duration
fs = 512e6;  % Sampling frequency
fc = 100;    % Carrier frequency
beta = B / T;
K = 3;

t = 0:1/fs:(T - 1/fs);  % Time vector one chirp
fi = beta * t;
y = cos(2 * pi * fc * t + beta * pi * t.^2);


% Méthode pour créer K chirps
t = 0:1/fs:(K*T - 1/fs); % Time vector K chirps
fi = repmat({fi}, 1, K); fi = cell2mat(fi); % Duplication and translation 
y = repmat({y}, 1, K); y = cell2mat(y);


N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  % Frequency vector

Y = fft(y);
Spectrum = abs(fftshift(Y));
threshold = max(Spectrum) / sqrt(2);  % Threshold at -3dB
indices = find(Spectrum > threshold);  % Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  % Compute the bandwidth
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);

figure;
subplot(3, 1, 2);  % Create three subplots vertically, this is the first one
plot(t, y);
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(3, 1, 1);  % This is the second subplot
plot(t, fi);
title('fi = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on

subplot(3, 1, 3);  % This is the third subplot
plot(frequencies, abs(fftshift(Y)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on