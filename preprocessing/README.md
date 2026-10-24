# Preprocessing (EEGLAB)

`eeg_preprocessing_using_interpolation_V1.m` preprocesses the stroke (Figshare) and
PhysioNet (run R08) recordings. Run it once per dataset, setting at the top:

| Dataset | `DATA_DIR` | `CHAN_LOCS` |
|---|---|---|
| Stroke | folder with the 36 stroke EDF files | <stroke location file> |
| PhysioNet | folder with the 18 R08 EDF files | Standard-10-10-Cap64.locs |
