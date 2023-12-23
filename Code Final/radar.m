close all; clear ; clc; 

%----Simulation parameters----
K = 26; % Slow-time FFT size (aka number of chirps in the matrix)
N = 52; %fast time FFT size (aka number of samples per chirp in the matrix)
T = 0.4e-3; %chirp duration -> produce (204800 ≈ 2^18) samples per chirp
Fss = 512e6; Tss = 1/Fss; %simulation sampling frequency
t_prime = 0:Tss:(T - Tss); %time vector one chirp
t = 0:Tss:(K*T -  Tss); %time vector K chirps


%----Radar parameters----
B = 200e6; %frequency range
fc = 24e9; %carrier frequency
beta = B/T;
alpha = 0.5; %echo signal amplitude
kapa = 0.25;


%----Target parameters----
scale = 1e5;
c = 3e8; %speed of light
r_max = 20; %initial range of the target
v_max = 2; %velocity of the target

target_1 = [r_max; v_max]; 
target_2 = [r_max; v_max]*0.9;
target_3 = [r_max; v_max]*0.8;
target_4 = [r_max; v_max]*0.4;
target_5 = [r_max; v_max]*0.5;
target_6 = [r_max; v_max]*0.6;
target_7 = [r_max; v_max]*0.7;
target_8 = [r_max; v_max]*0.8;
target_9 = [r_max; v_max]*0.9;

target__1 = [r_max; v_max]*rand; 
target__2 = [r_max; v_max]*rand;
target__3 = [r_max; v_max]*rand;
target__4 = [r_max; v_max]*rand;
target__5 = [r_max; v_max]*rand;
target__6 = [r_max; v_max]*rand;
target__7 = [r_max; v_max]*rand;
target__8 = [r_max; v_max]*rand;
target__9 = [r_max; v_max]*rand;


targets = {target_1*scale, target_2*scale, target_3*scale};
%targets = {target__1*scale, target__2*scale, target__3*scale, target__4*scale, target__5*scale, target__6*scale, target__7*scale, target__8*scale, target__9*scale};
n_targets = length(targets);


%----Emitted signal----
fi0 = @(t) beta * t;
fi_emitted = duplicate(fi0(t_prime), K);
for k = 0:K-1
    phi_emitted(k+1,:) = pi*k*beta*T^2 + pi*beta*t_prime.^2;
end
phi_emitted = reshape(phi_emitted', 1, []);
transmitted_signal = cos(2*pi*fc*t + phi_emitted); %transmitted signal s(t)


%----Received signal for each target----
fi_received_array = zeros(n_targets, length(t));
phi_received_array = zeros(n_targets, length(t));
received_signal_array = zeros(n_targets, length(t));
video_signal_1_array = zeros(n_targets, length(t));
%video_signal_2_array = zeros(n_targets, length(t));

delay_max = 0;
for i = 1:n_targets
    target_current = targets{i};
    f_B = (2*target_current(1)*beta)/c; f_D = (2*target_current(2)*fc)/c;
    delay_target_current = @(t) (2/c)*(target_current(1) + target_current(2)*t);
    delay_max = max(delay_target_current(t));
    %----Received signal with the target "i"----
    fi_received_current = duplicate(fi0(t_prime-delay_target_current(t_prime)), K);
    for k = 0:K-1
        phi_received_current_current(k+1,:) = pi*k*beta*T^2 + pi*beta*(t_prime-delay_target_current(t_prime)).^2; %#ok<SAGROW>
    end
    phi_received_current = reshape(phi_received_current_current', 1, []);
    received_signal_current = alpha * cos(2*pi*fc*(t-delay_target_current(t)) + phi_received_current); %echo signal
    
    %----Video signal with the target "i"----
    for k = 0:K-1
        video_signal_1_current_current(k+1,:) = kapa * exp(1j*2*pi*f_B*t_prime) .* exp(1j*2*pi*f_D*k*T); %#ok<SAGROW> %video signal x(t)
    end
    video_signal_1_current = reshape(video_signal_1_current_current', 1, []);
    %video_signal_2_current = (alpha/2) * exp(-1j * 2 * pi * fc * delay_target_current(t)) .* exp(1j * (phi_received_current - phi_emitted));

    fi_received_array(i,:) = fi_received_current;
    phi_received_array(i,:) = phi_received_current;
    received_signal_array(i,:) = received_signal_current; %echo signal
    video_signal_1_array(i,:) = video_signal_1_current;
    %video_signal_2_array(i,:) = video_signal_2_current;

end

if(n_targets>1)
    fi_received = sum(fi_received_array);
    phi_received = sum(phi_received_array);
    received_signal = sum(received_signal_array);
    video_signal_1 =  sum(video_signal_1_array);
    %video_signal_2 =  sum(video_signal_2_array);
else
    fi_received = fi_received_array(1,:);
    phi_received = phi_received_array(1,:);
    received_signal = received_signal_array(1,:);
    video_signal_1 =  video_signal_1_array(1,:);
    %video_signal_2 = video_signal_2_array(1,:);
end


%----radar processing----
video_signal_matrix_1 = reshape(video_signal_1(1:K*N), N, K); %video_signal_matrix_2 = reshape(video_signal_2(1:K*N), N, K); %creating a (N x K) matrix
range_doppler_map_1 = abs(fft2(video_signal_matrix_1)); %range_doppler_map_2 = abs(fft2(video_signal_matrix_2));


%----Fast fourier transforms computing----
N_samples = length(t);
frequencies = (-N_samples/2:N_samples/2-1)*(Fss/N_samples);  %Frequency vector
S = fft(transmitted_signal); S_shift = fftshift(S);
R = fft(received_signal); R_shift = fftshift(R);
X1 = fft(video_signal_1); X1_shift = fftshift(X1);
%X2 = fft(video_signal_2); X2_shift = fftshift(X2);

threshold_emitted = max(abs(S_shift))/sqrt(2); threshold_received = max(abs(R_shift))/sqrt(2);
threshold_video_1 =  max(abs(X1_shift))/sqrt(2); %threshold_video_2 =  max(abs(X2_shift))/sqrt(2);

indexes_emitted = find(abs(S_shift) > threshold_emitted); indexes_received = find(abs(R_shift) > threshold_received);
indexes_video_1 = find(abs(X1_shift) > threshold_video_1); %indexes_video_2 = find(abs(X2_shift) > threshold_video_2);

bandwidth_emitted = (frequencies(indexes_emitted(end)) - frequencies(indexes_emitted(1)))/2; bandwidth_received = (frequencies(indexes_received(end)) - frequencies(indexes_received(1)))/2;
bandwidth_video_1 = (frequencies(indexes_video_1(end)) - frequencies(indexes_video_1(1)))/2; %bandwidth_video_2 = (frequencies(indexes_video_2(end)) - frequencies(indexes_video_2(1)))/2;


%----Printing----
fprintf('S: amplitude_min = %f, amplitude_max = %f, BWD = %f MHz.\n', min(abs(S)), max(abs(S)), bandwidth_emitted/1e6);
fprintf('R: amplitude_min = %f, amplitude_max = %f, BWD = %f MHz.\n\n', min(abs(R)), max(abs(R)), bandwidth_received/1e6);
fprintf('X1: amplitude_min = %f, amplitude_max = %f, BWD = %f Hz.\n', min(abs(X1)), max(abs(X1)), bandwidth_video_1);
%fprintf('X2: amplitude_min = %f, amplitude_max = %f, BWD = %f Hz.\n\n', min(abs(X2)), max(abs(X2)), bandwidth_video_2)


%----Ploting----
plotting_time_domain(t, fi_emitted, fi_received, phi_emitted, phi_received, transmitted_signal, received_signal, video_signal_1);
plotting_frequency_domain(frequencies, S, R, 'Emitted', 'Received');
%plotting_frequency_domain(frequencies, X1, X2, 'X1', 'X2');

%plotting_RDM(range_doppler_map_1, 'video 1');
%plotting_RDM(range_doppler_map_2, 'video 2');

