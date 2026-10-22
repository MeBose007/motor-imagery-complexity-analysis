"""Reference signals for Fig. S2 (periodic and 1/f noise), generated as in the authors'
code but with a fixed random seed."""
import numpy as np


def sine(fs, duration_s, freq=1.0, amplitude=0.002):
    t = np.arange(0, duration_s, 1.0 / fs)
    return amplitude * np.sin(2 * np.pi * freq * t)


def pink_noise(fs, duration_s, rng):
    """1/f noise: Gaussian Fourier amplitudes with variance 1/f, zero DC, inverse FFT.
    (The authors drew the phase values from a standard normal distribution; kept here.)"""
    n = len(np.arange(0, duration_s, 1.0 / fs))
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    f[0] = 1.0
    power = 1.0 / f
    power[0] = 0.0
    amp = np.sqrt(power) * rng.standard_normal(len(f))
    phase = rng.standard_normal(len(f))
    return np.fft.irfft(amp * np.exp(1j * phase), n=n)[:n]


def white_noise(fs, duration_s, rng):
    return rng.standard_normal(len(np.arange(0, duration_s, 1.0 / fs)))
