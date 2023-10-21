close all; clear; clc;
%y0 = cos(2 * pi * fc * t_prime + phi0);
B = 200e3;   % Frequency range
T = 0.4e-3;  % Chirp duration
fs = 512e6;  % Sampling frequency
fc = 100;    % Carrier frequency
K = 2;       % Number of chirps
beta = B/T;
delay = 0.4e-3/15;
t_prime = 0:1/fs:(T - 1/fs);  % Time vector one chirp
t = 0:1/fs:(K*T - 1/fs); % Time vector K chirps

fi0 = @(t_prime) beta * t_prime;

fi = duplicate(fi0(t_prime), K);
fi_received = duplicate(fi0(t_prime-delay), K);


phi_0 = beta * pi * t_prime.^2;
phi = generate_phases(phi_0, K, beta, T);

phi_01 = beta * pi * (t_prime-delay).^2;
phi_1 = generate_phases(phi_01, K, beta, T-delay);

transmitted_signal = @(t) cos(2 * pi * fc * t + phi);
received_signal = @(t) cos(2 * pi * fc * (t-delay) + phi_1);

N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  % Frequency vector

Y = fft(transmitted_signal(t));
Spectrum = abs(fftshift(Y));
threshold = max(Spectrum) / sqrt(2);  % Threshold at -3dB
indices = find(Spectrum > threshold);  % Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  % Compute the bandwidth
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);



figure;
subplot(4, 1, 1);  % This is the second subplot
plot(t, fi, t+delay, fi_received);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on

subplot(4, 1, 2);  % This is the second subplot
plot(t, phi, t+delay, phi_1);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
grid on

subplot(4, 1, 3);  % Create three subplots vertically, this is the first one
plot(t, transmitted_signal(t), t+delay, received_signal(t));
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

