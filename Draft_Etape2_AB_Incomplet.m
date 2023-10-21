close all;
clear;
clc;

B = 200e6;   % Frequency range
T = 0.4e-3;  % Chirp duration
fs = 512e6;  % Sampling frequency
fc = 24e9;    % Carrier frequency
K = 10;       % Number of chirps
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

%fi = repmat(fi0, 1, K);
%phi = reshape(phi, 1, []);


fi = duplicate(fi0, K);
phi = merge_mat(phi);
y = cos(2 * pi * fc * t + phi);



delay = 1.6e-3;  % Delay in seconds - should affect the Fourier transform
%delay = 0;
Attenuation = 0.1;  % Attenuation factor
z = Attenuation * y .* (heaviside(t - delay) - heaviside(t - delay - T));

N_samples = length(t);
frequencies = (-N_samples/2:N_samples/2-1) * fs / N_samples;  % Frequency vector

Y = fftshift(fft(y))/N_samples;
Z = fftshift(fft(z))/N_samples;
Spectrum_Y = abs(Y);
Spectrum_Z = abs(Z);
threshold = max(Spectrum_Y) / sqrt(2);  % Threshold at -3dB
indices = find(Spectrum_Y > threshold);  % Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  % Compute the bandwidth
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);

figure;

subplot(6, 1, 1);  % FIT
plot(t, fi);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on

subplot(6, 1, 2);  % PHI
plot(t, phi);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
grid on

subplot(6, 1, 3);  % Transmitted
plot(t, y);
title('Transmitted Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(6, 1, 4);  % Received
plot(t, z);
title('Received Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(6, 1, 5);  % Fourier Transmitted
plot(frequencies, Spectrum_Y)
title('Frequency Domain Spectrum of Transmitted Signal')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on

subplot(6, 1, 6);  % Fourier Received
plot(frequencies, Spectrum_Z)
title('Frequency Domain Spectrum of Received Signal')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on
