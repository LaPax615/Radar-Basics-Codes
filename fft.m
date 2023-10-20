close all; clear; clc;

B = 200e6; %frequency range
T = 0.4e-3; %chirp duration
fs = 512e6; %sampling frequency
fc = 100;  %carrier frequency
beta = B/T;
 
duration = T; %signal duration
t = 0:1/fs:(duration - 1/fs); %time vector
N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  %frequencies vector

y = cos(2*pi*fc*t + beta*pi*t.^2);

figure;
subplot(2, 1, 1);  % Create two subplots vertically, this is the first one
plot(t, y);
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

Y = fft(y);

% Create a new figure for the frequency domain plot
subplot(2, 1, 2);  % This is the second subplot
plot(frequencies, abs(fftshift(Y)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on