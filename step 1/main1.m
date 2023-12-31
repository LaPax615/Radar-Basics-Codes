close all; clear ; clc;

%----Radar parameters----
F = 512e6; %simulation sampling frequency
B = 200e3; %frequency range
T = 0.4e-3; %chirp duration -> produce (204800 ≈ 2^18) samples per chirp
beta = B/T;
K = 3; % number of chirps
fc = 100; %carrier frequency
t_prime = 0:1/F:(T - 1/F); %time vector one chirp
t = 0:1/F:(K*T - 1/F); %time vector K chirps



%----evolution of the carrier frequency with time: f(t)----
fi0 = @(t_prime) beta * t_prime;
fi_emitted = duplicate(fi0(t_prime), K);


%----φ(t)----
phi_emitted = 2*pi*cumtrapz(t, fi_emitted); %integrale


%----radar processing----
transmitted_signal = @(t) cos(2 * pi * fc * t + phi_emitted); %Transmitted signal s(t)


%----Fast fourier transforms computing----
N_samples_tot = length(t);
frequencies = F * (-N_samples_tot/2:N_samples_tot/2-1) / N_samples_tot;  %Frequency vector
S = fft(transmitted_signal(t));
Spectrum = abs(fftshift(S));
threshold = max(Spectrum) / sqrt(2);  %Threshold at -3dB
indices = find(Spectrum > threshold);  %Find indices above the threshold
bandwidth = (frequencies(indices(end)) - frequencies(indices(1)))/2;  %Compute the bandwidth


%----Printing----
fprintf('The bandwidth of the signal is %f MHz.\n', bandwidth/1e6);
fprintf('T = %f.\n', T);
fprintf('N_sample_1_chirp = %f.\n', N_samples_tot/K);
fprintf('N_sample_tot = %f.\n', N_samples_tot);


%----Plotting----
figure;
subplot(3,1,1)
plot(t, fi_emitted);
title('fi(t) = beta*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted')
grid on

subplot(3,1,2)
plot(t, transmitted_signal(t));
title('Time Domain Signal')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Emitted')
grid on

subplot(3,1,3)
plot(frequencies, abs(fftshift(S)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
legend('Emitted')
grid on