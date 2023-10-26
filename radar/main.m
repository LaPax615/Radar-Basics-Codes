close all; clear ; clc;

%----Constants declaration----
c = 3e8; %speed of light
R_0 = 1e4; %initial range of the object
v = 1e6; %velocity of the object
B = 200e6; %frequency range
F = 512e6; %simulation sampling frequency
fc = 24e9; %carrier frequency
K = 2; %number of chirps
T = 0.4e-3; %chirp duration
N = 512; %fast time FFT size (aka number of samples per chirp in the matrix)
beta = B/T;
alpha = 0.5; %echo signal amplitude
t_prime = 0:1/F:(T - 1/F); %time vector one chirp
t = 0:1/F:(K*T - 1/F); %time vector K chirps
%delay = @(t) t*.100; %delay = @(t) (2/c)*(R_0 + v*t);
delay = @(t) (2*R_0)/c + t.*0;

%----evolution of the carrier frequency with time: f(t)----
fi0 = @(t_prime) beta * t_prime;
fi_emitted = duplicate(fi0(t_prime), K);
fi_received = duplicate(fi0(t_prime-delay(t_prime)), K);

%----φ(t)----
phi_emitted = 2*pi*cumtrapz(t, fi_emitted); %integrale
phi_received = 2*pi*cumtrapz(t, fi_received);
%phi_0 = beta * pi * t_prime.^2;
%phi_01 = beta * pi * (t_prime-delay(t_prime)).^2;
%phi_emitted = generate_phases(phi_0, K, beta, T);
%phi_received = generate_phases(phi_01, K, beta, T-delay(t_prime));

%----radar processing----
transmitted_signal = @(t) cos(2 * pi * fc * t + phi_emitted); %Transmitted signal s(t)
received_signal = @(t) alpha * cos(2 * pi * fc * (t-delay(t)) + phi_received); %echo signal
video_signal =  @(t) (alpha/2) * exp(-1j * 2 * pi * fc * delay(t)) .* exp(1j * (phi_received - phi_emitted)); %video signal x(t)
x = video_signal(t);
x_matrix = reshape(x(1:K*N), K, N);
range_doppler_map = abs(fft2(x_matrix));

%----Fast fourier transforms computing----
N_samples_tot = length(t);
frequencies = F * (-N_samples_tot/2:N_samples_tot/2-1) / N_samples_tot;  %Frequency vector
S = fft(transmitted_signal(t));
R_0 = fft(received_signal(t));
X = fft(video_signal(t));
Spectrum = abs(fftshift(S));
threshold = max(Spectrum) / sqrt(2);  %Threshold at -3dB
indices = find(Spectrum > threshold);  %Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  %Compute the bandwidth
delay_max = T - N/F;
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);
fprintf('T = %f.\n', T);
fprintf('K = %f. \n', K);
fprintf('N = %f.\n', N);
fprintf('delay_max = %f.\n\n', delay_max);
fprintf('N_sample_1_chirp = %f.\n', N_samples_tot/K);


%----Plotting----
figure
colormap(jet); % Or any other colormap of your choice
imagesc(range_doppler_map);
xlabel('Doppler Bins');
ylabel('Range Bins');
colorbar; % Add a colorbar to the plot

figure;
subplot(4, 1, 1);  %Create three subplots vertically, this is the first one
plot(t, transmitted_signal(t), t, received_signal(t));
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Emitted', 'Received')
grid on

subplot(4, 1, 2);  %This is the third subplot
plot(frequencies, abs(fftshift(S)), frequencies, abs(fftshift(R_0)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
legend('Emitted', 'Received')
grid on

subplot(4, 1, 3);  %Create three subplots vertically, this is the first one
plot(t, video_signal(t), 'g');
title('video signal')
xlabel('Time (s)')
ylabel('Amplitude')
grid on

subplot(4, 1, 4);  %Create three subplots vertically, this is the first one
plot(frequencies, abs(fftshift(X)))
title('video signal spectrum')
xlabel('Frequency (Hz))')
ylabel('|X(f)|')
grid on


figure;
subplot(2, 1, 1);  %This is the second subplot
plot(t, fi_emitted, t+delay(t), fi_received);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted', 'Received')
grid on

subplot(2, 1, 2);  %This is the second subplot
plot(t, phi_emitted, t+delay(t), phi_received);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
legend('Emitted', 'Received')
grid on


%{
figure;
subplot(6, 1, 1);  %This is the second subplot
plot(t, fi_emitted, t+delay(t), fi_received);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted', 'Received')
grid on

subplot(6, 1, 2);  %This is the second subplot
plot(t, phi_emitted, t+delay(t), phi_received);
title("φ(t) = π·k·β·T² + π·β·t'")
xlabel('Time (s)')
ylabel('phi')
legend('Emitted', 'Received')
grid on

subplot(6, 1, 3);  %Create three subplots vertically, this is the first one
plot(t, transmitted_signal(t), t, received_signal(t));
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Emitted', 'Received')
grid on

subplot(6, 1, 4);  %This is the third subplot
plot(frequencies, abs(fftshift(S)), frequencies, abs(fftshift(R_0)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
legend('Emitted', 'Received')
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
%}
