function output = plotting_time_domain(t, fi_emitted, fi_received, phi_emitted, phi_received, transmitted_signal, received_signal, video_signal_1);

output = 0;

figure
subplot(4,1,1)
plot(t/1e-3, fi_emitted/1e6, t/1e-3, fi_received/1e6);
title('fi(t) = β*t')
xlabel('Time (ms)')
ylabel('Frequency (MHz)')
legend('Emitted', 'Received')
grid on

subplot(4,1,2)
plot(t/1e-3, phi_emitted, t/1e-3, phi_received);
title("φ(t) = π·k·β·T² + π·β·t'²")
xlabel('Time (ms)')
ylabel('Phase')
legend('Emitted', 'Received')
grid on

subplot(4,1,3)
plot(t/1e-3, transmitted_signal, t/1e-3, received_signal);
title('signals')
xlabel('Time (ms)')
ylabel('Amplitude')
legend('Emitted', 'Received')
grid on

subplot(4,1,4)
plot(t/1e-3, video_signal_1)
title('Video signal 1')
xlabel('Time (ms)')
ylabel('Amplitude')
legend('x1')
grid on

end