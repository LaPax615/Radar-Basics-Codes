function output = plotting_frequency_domain(frequencies, spectrum1, spectrum2, legend1, legend2)
output = 0;
figure
plot(frequencies/1e6, abs(fftshift(spectrum1)), frequencies/1e6, abs(fftshift(spectrum2)))
title('Frequency Domain Spectrum')
xlabel('Frequency (MHz)')
ylabel('Amplitude')
legend(legend1, legend2)
grid on
end