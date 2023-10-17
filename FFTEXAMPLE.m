


T = 100*(1/100);
fs = 100;
t = 0:1/fs:T-1/fs;
f_sawtooth = 20;


f = @(T, fs, f_sawtooth) 0.5*sawtooth(2*pi*f_sawtooth*t);
g = @(t) pi*0.5*t.^2;

y = f(T, fs, f_sawtooth);
z = g(t);
transmitted_signal = y .* exp(1i*z);


L = length(transmitted_signal);
Y = fft(transmitted_signal);

%Frequency range
Fs = 1 / (t(2) - t(1));  % Sampling frequency
frequencies = Fs * (0:(L/2)) / L;

% FFT plot
plot(frequencies, abs(Y(1:L/2+1)), 'b', 'LineWidth', 1.5);

%plot(t, Y, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Transmitted Signal');
grid on;

% BAndwidth calcus
bandwidth = max(frequencies) - min(frequencies);
disp("Bandwidth:")
disp(bandwidth);
