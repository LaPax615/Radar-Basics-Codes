function output = plotting_time_domain(t, fi_emitted, fi_received, phi_emitted, phi_received, transmitted_signal, received_signal, video_signal_1);

output = 0;

figure
subplot(4,1,1)
plot(t, fi_emitted, t, fi_received);
title('fi(t) = β*t')
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted', 'Received')
grid on

subplot(4,1,2)
plot(t, phi_emitted, t, phi_received);
title("φ(t) = π·k·β·T² + π·β·t'²")
xlabel('Time (s)')
ylabel('fi (Hz)')
legend('Emitted', 'Received')
grid on

subplot(4,1,3)
plot(t, transmitted_signal, t, received_signal);
title('signals')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Emitted', 'Received')
grid on

subplot(4,1,4)
plot(t, video_signal_1)
title('Video signal 1')
xlabel('Time (s)')
ylabel('Amplitude')
legend('x1')
grid on

end