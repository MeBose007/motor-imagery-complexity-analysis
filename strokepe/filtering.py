"""Zero-phase Butterworth band-pass filter (used for the 8-13 Hz analysis, Fig. S3)."""
from scipy.signal import butter, filtfilt


def bandpass(data, low, high, fs, order=4, axis=-1):
    nyq = 0.5 * fs
    b, a = butter(order, [low / nyq, high / nyq], btype="band")
    return filtfilt(b, a, data, axis=axis)
