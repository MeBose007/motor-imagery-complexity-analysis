"""Readers for the three EEG datasets. All return arrays shaped (channels, samples)
with channels in the canonical order of config/channels.yaml."""
import warnings
import numpy as np
import pandas as pd


def _key(name):
    return "".join(ch for ch in str(name) if ch.isalnum()).lower()


# ---------------------------------------------------------------- stroke (CSV, 500 Hz)
def load_stroke_csv(path, drop_columns, canonical):
    """Preprocessed stroke recording: one column per channel, header = channel names."""
    df = pd.read_csv(path)
    drop = set(drop_columns)
    df = df.iloc[:, [i for i in range(df.shape[1]) if i not in drop]]
    lookup = {_key(c): c for c in df.columns}
    missing = [ch for ch in canonical if _key(ch) not in lookup]
    extra = [c for c in df.columns if _key(c) not in {_key(ch) for ch in canonical}]
    if missing or extra:
        raise ValueError(f"{path}: after removing columns {sorted(drop)} expected the 25 "
                         f"canonical channels; missing {missing}, unexpected {extra}")
    return df[[lookup[_key(ch)] for ch in canonical]].to_numpy(dtype=float).T


# ---------------------------------------------------------------- EDF (PhysioNet, Dortmund)
def read_edf(path):
    """Return (data, labels, sampling_rates) for every signal in an EDF file."""
    import pyedflib
    f = pyedflib.EdfReader(str(path))
    try:
        data = np.array([f.readSignal(i) for i in range(f.signals_in_file)])
        labels = f.getSignalLabels()
        rates = list(f.getSampleFrequencies())
    finally:
        f.close()
    return data, labels, rates


def load_edf_by_index(path, index_map, canonical, fs_expected=None):
    data, _, rates = read_edf(path)
    idx = [index_map[ch] for ch in canonical]
    if fs_expected is not None and any(abs(rates[i] - fs_expected) > 1e-6 for i in idx):
        warnings.warn(f"{path}: sampling rate {rates[idx[0]]} Hz, expected {fs_expected} Hz")
    return data[idx]


def load_edf_by_name(path, names):
    data, labels, rates = read_edf(path)
    lookup = {_key(lab): i for i, lab in enumerate(labels)}
    missing = [n for n in names if _key(n) not in lookup]
    if missing:
        raise ValueError(f"{path}: channels {missing} not found among {labels}")
    idx = [lookup[_key(n)] for n in names]
    return data[idx], rates[idx[0]]


# ---------------------------------------------------------------- EEGLAB .set (OpenNeuro SRM)
def read_set_epochs(path):
    """Return the EEGLAB data array as (channels, samples_per_epoch, epochs)."""
    import scipy.io
    data = None
    try:
        mat = scipy.io.loadmat(str(path))
        if "data" in mat and np.issubdtype(np.asarray(mat["data"]).dtype, np.number):
            data = np.asarray(mat["data"], dtype=float)
    except (NotImplementedError, ValueError):
        pass  # MATLAB v7.3 file or data stored in a separate .fdt file
    if data is None:
        import mne
        try:
            ep = mne.read_epochs_eeglab(str(path), verbose="error")
            data = ep.get_data().transpose(1, 2, 0) * 1e6
        except Exception:
            raw = mne.io.read_raw_eeglab(str(path), preload=True, verbose="error")
            data = raw.get_data()[:, :, None] * 1e6
    if data.ndim == 2:
        data = data[:, :, None]
    return data


def concat_epochs(epochs):
    """(channels, samples, epochs) -> (channels, samples*epochs), epoch after epoch."""
    ch, n, ne = epochs.shape
    return epochs.transpose(0, 2, 1).reshape(ch, n * ne)


def select_by_index(data, index_map, canonical):
    return data[[index_map[ch] for ch in canonical]]
