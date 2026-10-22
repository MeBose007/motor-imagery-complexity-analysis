"""Subject-level (PE_sub) and sensor-level (PE_sensor) averages, Sec. II B."""
import pandas as pd

KEYS = ["dataset", "subject", "severity", "lesion", "condition"]


def subject_channel(df):
    """Mean over segments (a no-op in 'paper' mode, where there is one segment)."""
    return df.groupby(KEYS + ["channel"], sort=False)["pe"].mean().reset_index()


def pe_sub(df):
    """PE_sub: one value per subject and condition = mean over the 25 channels."""
    return subject_channel(df).groupby(KEYS, sort=False)["pe"].mean().reset_index()


def pe_sensor(df, condition, canonical, **filters):
    """PE_sensor: mean over subjects for each channel, in canonical channel order."""
    sc = subject_channel(df)
    sc = sc[sc["condition"] == condition]
    for col, val in filters.items():
        sc = sc[sc[col] == val]
    return sc.groupby("channel")["pe"].mean().reindex(canonical)


def select(sub, condition, **filters):
    out = sub[sub["condition"] == condition]
    for col, val in filters.items():
        out = out[out[col] == val]
    return out
