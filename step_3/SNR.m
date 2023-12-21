% Radar Specifications
Range_max = 109.1;  % Maximum range in meters
Range_resolution = 1;  % Range resolution in meters
Velocity_max = 10;  % Maximum velocity in m/s
c = 3e8;  % Speed of light in m/s

% FMCW Waveform Generation
B = 200e6;  % Bandwidth
Tchirp = 5.5 * (2 * Range_max / c);  % Chirp Time
disp(Tchirp);

slope = B / Tchirp;
fc = 24e9;  % Carrier frequency: 24 GHz

Nr = 256;  % Number of range cells
Nd = 512;  % Number of Doppler cells
SweepTime = Nr * Tchirp;

t = linspace(0, SweepTime, Nr * Nd);  % Total time for Nr sweeps

% Define multiple targets
targets = [
    struct('range', 50, 'velocity', 5),
    struct('range', 100, 'velocity', -3),
    struct('range', 150, 'velocity', 2),
    struct('range', 75, 'velocity', 1),
    struct('range', 125, 'velocity', -2),
    struct('range', 60, 'velocity', 0)
];

% Compute FFT

% Define velocity and range bins based on your data
velocity_bins = linspace(-Velocity_max, Velocity_max, Nd);
range_bins = linspace(0, Range_max, Nr);

% Noise power as used in your code
noise_power = 0.01;
SNRs = [10, 1, 0.01, 0.001, 0.0001, 0.00001];



% Calculate probabilities and plot
for i = 1:length(SNRs)
    RDM = generate_rdm_for_snr(targets, SNRs(i), Nr, Nd, t, fc, c, slope, noise_power);
    subplot(2, 3, i);
    imagesc(velocity_bins, range_bins, RDM);
    colormap('jet');
    title(sprintf('RDM avec SNR = %.2f dB', 10 * log10(SNRs(i))));
    xlabel('V axis');
    ylabel('R axis');
end



function RDM = compute_fft(mix_signal, Nr, Nd)
    fft2d = fftshift(fft2(reshape(mix_signal, Nr, Nd)));
    RDM = abs(fft2d);
end


function RDM = generate_rdm_for_snr(targets, SNR, Nr, Nd, t, fc, c, slope, noise_power)
    Tx = zeros(1, length(t));
    Rx = zeros(1, length(t));
    Mix = zeros(1, length(t));

    for i = 1:length(targets)
        Target_range = targets(i).range;
        Target_velocity = targets(i).velocity;
        f_doppler = 2 * Target_velocity * fc / c;

        for j = 1:length(t)
            range_at_t = Target_range + (Target_velocity * t(j));
            delay = 2 * range_at_t / c;

            Tx_signal = exp(1i * 2 * pi * (fc * t(j) + 0.5 * slope * t(j) ^ 2));
            Rx_signal = exp(1i * 2 * pi * ((fc + f_doppler) * (t(j) - delay) + 0.5 * slope * (t(j) - delay) ^ 2));

            Tx(j) = Tx(j) + Tx_signal;
            % Adjust noise level based on SNR
            noise = (randn(1) + 1i * randn(1)) * sqrt(noise_power / SNR);
            Rx(j) = Rx(j) + Rx_signal + noise;
            Mix(j) = Mix(j) + Tx_signal * conj(Rx(j));
        end
    end

    RDM = compute_fft(Mix, Nr, Nd);
end

