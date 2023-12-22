close all; clear ; clc; 
K = 512;
T = 0.4e-3;
Fs = 512e6; Ts = 1/Fs;
t_prime = 0:Ts:(T - Ts); %time vector one chirp
t = 0:Ts:(K*T -  Ts); %time vector K chirps
fc = 24e6;

x = sin(2*pi*fc*t);
N_samples = length(t);

fshift = (-N_samples/2:N_samples/2-1)*(Fs/N_samples);
y = fft(x);
yshift = fftshift(y);


figure
subplot(2,1,1)
plot(t, x);
title('x')
xlabel('Time (s)')
ylabel('s(t)')
legend('Emitted')
grid on

subplot(2,1,2)
plot(fshift/1e6, abs(yshift)); % Convert frequencies to megahertz
title('X')
xlabel('Frequencies (MHz)')
ylabel('|X|')
legend('Emitted')
grid on

%{

%----Fast fourier transforms computing----
N_samples_tot = length(t);
frequencies = F * (-N_samples_tot/2:N_samples_tot/2-1) / N_samples_tot;  %Frequency vector
S = fft(transmitted_signal);
R = fft(received_signal);
X = fft(video_signal);
threshold_emitted = max(abs(fftshift(S)))/sqrt(2); threshold_received = max(abs(fftshift(R)))/sqrt(2); threshold_video =  max(abs(fftshift(X)))/sqrt(2);  %Threshold at -3dB
index_emitted = find(abs(fftshift(S)) > threshold_emitted); index_received = find(abs(fftshift(R)) > threshold_received); index_video = find(abs(fftshift(X)) > threshold_video); %Find indices above the threshold
bandwidth_emitted = (frequencies(index_emitted(end)) - frequencies(index_emitted(1)))/2; bandwidth_received = (frequencies(index_received(end)) - frequencies(index_received(1)))/2; %Compute the bandwidth
bandwidth_video = (frequencies(index_video(end)) - frequencies(index_video(1)))/2;
delay_max = T - N/F;

%----Printing----
fprintf('S: f_max = %f Hz, f_min = %f Hz.\n', max(abs(S)), min(abs(S)));
fprintf('R: f_max = %f Hz, f_min = %f Hz.\n', max(abs(R)), min(abs(R)));
fprintf('X1: f_max = %f Hz, f_min = %f Hz.\n', max(abs(X1)), min(abs(X1)));
fprintf('X2: f_max = %f Hz, f_min = %f Hz.\n', max(abs(X2)), min(abs(X2)));


fprintf('BWD emitted = %f MHz.\n', bandwidth_emitted/1e6);
fprintf('BWD received = %f MHz.\n', bandwidth_received/1e6);
fprintf('BWD video = %f Hz.\n', bandwidth_video);
fprintf('T = %f.\n', T);
fprintf('N = %f.\n', N);
fprintf('K = %f.\n', K);
fprintf('delay_max = %f.\n\n', delay_max);
fprintf('N_sample_1_chirp = %f.\n', N_samples_tot/K);
fprintf('N_sample_tot = %f.\n', N_samples_tot);

%----Plotting----
figure
colormap(jet); % Or any other colormap of your choice
imagesc(range_doppler_map);
xlabel('K');
ylabel('N');
colorbar; % Add a colorbar to the plot


figure
plot(frequencies, abs(fftshift(X)))
title('Frequency Domain Spectrum')
xlabel('Frequency (Hz)')
ylabel('|S(f)|')
legend('Emitted', 'Video')
grid on

figure
subplot(5,1,1)
plot(t, fi_emitted);
title('fi(t) = β*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted')
grid on

subplot(5,1,2)
plot(t, phi_emitted);
title("φ(t) = π·k·β·T² + π·β·t'²")
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted')
grid on

subplot(5,1,3)
plot(t, transmitted_signal, t, received_signal);
title('signals')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Emitted', 'Received')
grid on

subplot(5,1,4)
plot(t, video_signal)
title('Video signal 1')
xlabel('Time (s)')
ylabel('Amplitude')
legend('x1')
grid on

subplot(5,1,5)
plot(t, video_signal_2)
title('Video signal 2')
xlabel('Time (s)')
ylabel('Amplitude')
legend('x2')
grid on


%}
