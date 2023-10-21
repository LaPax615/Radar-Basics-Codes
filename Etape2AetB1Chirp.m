close all;
clear;
clc;

B = 200e6;  % Frequency range
T = 0.4e-3;  % Chirp duration
fs = 512e6;  % Sampling frequency
fc = 24e2;    % Carrier frequency
beta = B / T;   % Modulation index

duration = T;  % Signal duration
t = 0:1/fs:(duration - 1/fs);  % Time vector for the sim
N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  % Frequency vector

fi = beta * t;
y = cos(2 * pi * fc * t + beta * pi * t.^2);

delay = 1.6e-4;  % Delay in seconds - should eclater le fourier if it works
Attenuation = 0.1;  % Attenuation factor minimizes/maximizes
y_channel = Attenuation * y .* (heaviside(t - delay));



Y = fft(y);
Z = fft(y_channel);
Spectrum = abs(fftshift(Y));
threshold = max(Spectrum) / sqrt(2);  % Threshold at -3dB
indices = find(Spectrum > threshold);  % Find indices above the threshold
bandwidth = frequencies(indices(end)) - frequencies(indices(1));  % Compute the bandwidth

fprintf('The bandwidth of the signal is %f Hz.\n', bandwidth);



figure;
subplot(4, 1, 1);  % Triangles thingys
plot(t, fi);
title('fi')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on


subplot(4, 1, 2);  %Time signal
plot(t, y);
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on


subplot(4, 1, 3);  %Fourier before channel
plot(frequencies, abs(fftshift(Y)))
title('Frequency Domain Spectrum Before Channel')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on

subplot(4, 1, 4);  % Fourier after channel
plot(frequencies, abs(fftshift(Z)))
title('Frequency Domain Spectrum After Channel')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on
