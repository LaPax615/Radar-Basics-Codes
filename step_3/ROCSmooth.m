clear;
close all;
clc;

% Radar Specifications
Range_max = 300; % Maximum Range in meters
Range_resolution = 1; % Range Resolution in meters
Velocity_max = 40; % Maximum Velocity in m/s
c = 3e8; % Speed of light in m/s

% Define Targets
Target_range = [60, 80, 120, 150, 175, 200, 220, 250, 270]; % Target ranges in meters
Target_velocity = [20, -15, 30, 10, -5, 5, -10, 15, -20]; % Target velocities in m/s

% Radar Signal Parameters
B = 200e6; % Bandwidth
Tsweep = 5.5;
Tchirp = 0.4e-6;
slope = B / Tchirp;
fc = 24e9; % Carrier frequency
Nd = 256; % Number of Doppler cells
Nr = 512; % Number of range cells

% Time vector
t = linspace(0, Nd * Tchirp, Nr * Nd); % Total time for samples

% Define axes for the RDM
doppler_axis = linspace(-Velocity_max, Velocity_max, Nd);
range_axis = linspace(0, Range_max, Nr/2);
        
% Different SNR values for ROC curves
SNR_dB_values = [20, 5, -5, -20];

%SNR_dB_values = [9.4, 8, 4, 2, -5, -10, -15];
threshold_levels = linspace(-50, 200, 1000);

% Initialize TPR and FPR values for ROC curves
TPR_values = zeros(length(SNR_dB_values), length(threshold_levels));
FPR_values = zeros(length(SNR_dB_values), length(threshold_levels));

% Manually compute indices for targets in the RDM
range_step = Range_max / (Nr/2);
velocity_step = (2 * Velocity_max) / Nd;

for snr_idx = 1:length(SNR_dB_values)
    SNR_dB = SNR_dB_values(snr_idx);
    
    % Signal Processing with current SNR
    Tx = zeros(1, length(t));
    Rx = zeros(1, length(t));
    Mix = zeros(1, length(t));

    for i = 1:length(t)
        Rx_i = zeros(1, length(Target_range));
        for j = 1:length(Target_range)
            range_at_t = 2 * (Target_range(j) + (Target_velocity(j) * t(i))) / c;
            Rx_i(j) = cos(2 * pi * (fc * (t(i) - range_at_t) + slope * (t(i) - range_at_t)^2 / 2));
        end

        Tx(i) = cos(2 * pi * (fc * t(i) + slope * t(i)^2 / 2));
        Rx(i) = sum(Rx_i); % Sum of received signals before adding noise
    end

    % Adding White Gaussian Noise using awgn function
    Rx = awgn(Rx, SNR_dB, 'measured'); % Rx signal with added noise
    Mix = Tx .* Rx; % Mixing the transmitted and received signals

    % Reshape the Mix signal and perform FFT
    Mix = reshape(Mix, [Nr, Nd]);
    sig_fft2 = fft2(Mix, Nr, Nd);
    sig_fft2 = sig_fft2(1:Nr/2, 1:Nd);
    sig_fft2 = fftshift(sig_fft2);
    RDM = abs(sig_fft2);
    RDM = 10 * log10(RDM);
    
    % Create a logical matrix for true targets in the RDM
    true_targets_RDM = false(size(RDM));
    for i = 1:length(Target_range)
        range_bin = round(Target_range(i) / range_step) + 1;
        velocity_bin = round((Target_velocity(i) + Velocity_max) / velocity_step) + 1;

        range_bin = max(1, min(range_bin, Nr/2));
        velocity_bin = max(1, min(velocity_bin, Nd));

        true_targets_RDM(range_bin, velocity_bin) = true;
    end
    
    % Calculate ROC curve for the current SNR
    for thresh_idx = 1:length(threshold_levels)
        threshold = threshold_levels(thresh_idx);
        detected = RDM > threshold;
        
        TP = sum(sum(detected & true_targets_RDM));
        FP = sum(sum(detected & ~true_targets_RDM));
        FN = sum(sum(~detected & true_targets_RDM));
        TN = sum(sum(~detected & ~true_targets_RDM));
        
        TPR_values(snr_idx, thresh_idx) = TP / (TP + FN);
        FPR_values(snr_idx, thresh_idx) = FP / (FP + TN);
    end
end

% Plot ROC curves for each SNR
% Plot ROC curves for each SNR
figure; hold on;
colors = jet(length(SNR_dB_values));

% Number of points for interpolation
numInterpPoints = 100;

for snr_idx = 1:length(SNR_dB_values)
    % Original FPR and TPR values
    originalFPR = FPR_values(snr_idx, :);
    originalTPR = TPR_values(snr_idx, :);

    % Remove duplicate FPR points and corresponding TPR points
    [uniqueFPR, uniqueIdx] = unique(originalFPR, 'stable');
    uniqueTPR = originalTPR(uniqueIdx);

    % Interpolated FPR and TPR values
    interpFPR = linspace(min(uniqueFPR), max(uniqueFPR), numInterpPoints);
    interpTPR = interp1(uniqueFPR, uniqueTPR, interpFPR, 'spline');

    % Plot the interpolated ROC curve
    plot(interpFPR, interpTPR, 'Color', colors(snr_idx, :), 'LineWidth', 2);
end

xlabel('False Positive Rate - FPR');
ylabel('True Positive Rate - TPR');
title('ROC Curves for Different SNR Values');
legend(strcat('SNR = ', string(SNR_dB_values), ' dB'));
hold off;
