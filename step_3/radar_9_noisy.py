import numpy as np
import matplotlib.pyplot as plt

# Radar Specifications
Range_max = 109.1  # Maximum range in meters
Range_resolution = 1  # Range resolution in meters
Velocity_max = 10  # Maximum velocity in m/s
c = 3e8  # Speed of light in m/s

# FMCW Waveform Generation
B = 200e6  # Bandwidth
Tchirp = 5.5 * (2 * Range_max / c)  # Chirp Time
print(Tchirp)
#Tchirp = 0.4e-3

slope = B / Tchirp
fc = 24e9  # Carrier frequency: 24 GHz

Nr = 256  # Number of range cells
Nd = 512  # Number of Doppler cells
SweepTime = Nr * Tchirp

t = np.linspace(0, SweepTime, Nr * Nd)  # Total time for Nr sweeps

r = 40
v = 5


# Define multiple targets
targets = [
    {"range": 50, "velocity": 5},
    {"range": 100, "velocity": -3},
    {"range": 150, "velocity": 2},
    {"range": 75, "velocity": 1},
    {"range": 125, "velocity": -2},
    {"range": 60, "velocity": 0}
]


'''
targets = [
    {"range": r, "velocity": v},
    {"range": 2*r, "velocity": -1*v},
    {"range": 2*r, "velocity": 2*v},
    {"range": 4*r, "velocity": 3*v},
    {"range": 3*r, "velocity": -1*v},
    {"range": r, "velocity": v}
]'''


# Generate Tx and Rx signals (as in your code)

# Compute FFT
def compute_fft(mix_signal, Nr, Nd):
    fft2d = np.fft.fft2(mix_signal.reshape(Nr, Nd))
    fft2d = np.fft.fftshift(fft2d)
    RDM = np.abs(fft2d)
    return RDM

# Define velocity and range bins based on your data
velocity_bins = np.linspace(-Velocity_max, Velocity_max, Nd)
range_bins = np.linspace(0, Range_max, Nr)

# Noise power as used in your code
noise_power = 0.01
noises = [0.1, 0.2, 0.4, 1.2, 2]
# Assumed range of signal powers (you may adjust these based on your system)
signal_p = 1
signal_powers = [0.1, 0.5, 1, 5, 10]

# Calculating SNR values in dB
#SNRs = [10 * np.log10(signal_power / noise_power) for signal_power in signal_powers]

#SNRs = [10 * np.log10(signal_p / noise) for noise in noises]
SNRs = [10, 1, 0.01, 0.001, 0.0001, 0.00001]

def generate_rdm_for_snr(targets, SNR, Nr, Nd, t, fc, c, slope, noise_power):
    Tx = np.zeros(len(t), dtype=complex)
    Rx = np.zeros(len(t), dtype=complex)
    Mix = np.zeros(len(t), dtype=complex)

    for target in targets:
        Target_range = target["range"]
        Target_velocity = target["velocity"]
        f_doppler = 2 * Target_velocity * fc / c

        for i in range(len(t)):
            range_at_t = Target_range + (Target_velocity * t[i])
            delay = 2 * range_at_t / c

            Tx_signal = np.exp(1j * 2 * np.pi * (fc * t[i] + 0.5 * slope * t[i] ** 2))
            Rx_signal = np.exp(1j * 2 * np.pi * ((fc + f_doppler) * (t[i] - delay) + 0.5 * slope * (t[i] - delay) ** 2))

            Tx[i] += Tx_signal
            # Adjust noise level based on SNR
            Rx[i] += Rx_signal + (np.random.randn() + 1j * np.random.randn()) * np.sqrt(noise_power / SNR)
            Mix[i] += Tx_signal * np.conj(Rx[i])

    return compute_fft(Mix, Nr, Nd)






# Calculate probabilities (using the previously defined calculate_probabilities function)
plt.figure(figsize=(15, 10))
for i, snr in enumerate(SNRs, 1):
    RDM = generate_rdm_for_snr(targets, snr, Nr, Nd, t, fc, c, slope, noise_power)
    plt.subplot(2, 3, i)
    plt.imshow(RDM, cmap='jet', aspect='auto', extent=[velocity_bins[0], velocity_bins[-1], range_bins[0], range_bins[-1]])
    plt.title(f'RDM avec SNR = {10 * np.log10(snr):.2f} dB')
    plt.xlabel('V axis')
    plt.ylabel('R axis')

plt.tight_layout()
plt.show()
