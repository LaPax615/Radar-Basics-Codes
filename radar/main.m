close all; clear; clc;
c = 3e8; %speed of light
R_0 = 1e5; %initial range
B = 200e6;   %Frequency range
T = 0.4e-3;  %Chirp duration
fs = 512e6;  %Sampling frequency
fc = 100;    %Carrier frequency
K = 1;       %Number of chirps
beta = B/T;
alpha = 1;
t_prime = 0:1/fs:(T - 1/fs);  %Time vector one chirp
t = 0:1/fs:(K*T - 1/fs); %Time vector K chirps
delay = @(t) T/10 + t*0;
%delay = @(t) (2/c)*(R_0 + v*t);


fi0 = @(t_prime) beta * t_prime;
fi_emitted = duplicate(fi0(t_prime), K);
fi_received = duplicate(fi0(t_prime-delay(t_prime)), K);

phi_emitted = 2*pi*cumtrapz(t, fi_emitted); %Integrale
phi_received = 2*pi*cumtrapz(t, fi_received);

transmitted_signal = @(t) cos(2 * pi * fc * t + phi_emitted);
received_signal = @(t) alpha * cos(2 * pi * fc * (t-delay(t)) + phi_received);
video_signal =  @(t) (alpha/2) * exp(-1j * 2 * pi * fc * delay(t)) .* exp(1j * (phi_received - phi_emitted));


N_samples = length(t);
frequencies = fs * (-N_samples/2:N_samples/2-1) / N_samples;  %Frequency vector

S = fft(transmitted_signal(t));
R_0 = fft(received_signal(t));
X = fft(video_signal(t));

Spectrum = abs(fftshift(S));
threshold = max(Spectrum) / sqrt(2);  %Threshold at -3dB
indices = find(Spectrum > threshold);  %Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  %Compute the bandwidth
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);



figure;
subplot(6, 1, 1);  %This is the second subplot
plot(t, fi_emitted, t+delay(t), fi_received);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
grid on

subplot(6, 1, 2);  %This is the second subplot
plot(t, phi_emitted, t+delay(t), phi_received);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
grid on

subplot(6, 1, 3);  %Create three subplots vertically, this is the first one
plot(t, transmitted_signal(t), t, received_signal(t));
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(6, 1, 4);  %This is the third subplot
plot(frequencies, abs(fftshift(S)), frequencies, abs(fftshift(R_0)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
grid on

subplot(6, 1, 5);  %Create three subplots vertically, this is the first one
plot(t, video_signal(t), 'g');
title('video signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(6, 1, 6);  %Create three subplots vertically, this is the first one
plot(frequencies, abs(fftshift(X)))
title('video signal spectrum')
xlabel('Frequency (Hz))')
ylabel('|X(f)|')
grid on

