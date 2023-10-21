close all;
clear;
clc;

B = 200e2;   % Frequency range
T = 0.4e-3*10;  % Chirp duration
fs = 512e6;  % Sampling frequency
fc = 100;    % Carrier frequency
K = 4;       % Number of chirps
beta = B/T;


t_prime = 0:1/fs:(T - 1/fs);  % Time vector one chirp
t = 0:1/fs:(K*T - 1/fs); % Time vector K chirps

fi0 = beta * t_prime;
phi0 = beta * pi * t_prime.^2;
y0 = cos(2 * pi * fc * t_prime + phi0);

phi = zeros(K, length(phi0));
for i=1:K
    k = i-1;
    fprintf('k = %d\n', k);
    phi(i,:) = phi0 + k*pi*beta*T^2; 
end

fi = duplicate(fi0, K);
phi = merge_mat(phi);
y = cos(2 * pi * fc * t + phi);

N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  % Frequency vector

Y = fft(y);
Spectrum = abs(fftshift(Y));
threshold = max(Spectrum) / sqrt(2);  % Threshold at -3dB
indices = find(Spectrum > threshold);  % Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  % Compute the bandwidth
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);

figure;
subplot(4, 1, 1);  % This is the second subplot
plot(t, fi);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on

subplot(4, 1, 2);  % This is the second subplot
plot(t, phi);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
grid on

subplot(4, 1, 3);  % Create three subplots vertically, this is the first one
plot(t, y);
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on


subplot(4, 1, 4);  % This is the third subplot
plot(frequencies, abs(fftshift(Y)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on