# Data

No EEG data are stored in this repository. All datasets are public; download them and
set the folders in `config/paths.yaml` (copy `config/paths.example.yaml`).

| Folder (`paths.yaml` key) | Dataset | Used for | Expected files |
|---|---|---|---|
| `stroke_dir` | Acute stroke motor-imagery EEG, Xuanwu Hospital (Liu et al., *Sci. Data* 11, 131, 2024, doi:10.1038/s41597-023-02787-8), Figshare | all stroke results | `sub-XX_task-motor-imagery_eeg_cleaned_EOG_avgref_interp_preprocessed_final_channels_as_columns.csv` |
| `eo_mi_dir` | EEG Motor Movement/Imagery Database, PhysioNet (Schalk et al. 2004), https://physionet.org/content/eegmmidb/1.0.0/ | EO-MI | `SXXXR08_avgref_interp_preprocessed_final.edf` |
| `ec_rs_dir` | SRM resting-state EEG, OpenNeuro ds003775 (doi:10.18112/openneuro.ds003775.v1.2.1), cleaned epochs (`derivatives`) | EC-RS | `sub-XXX_ses-t1_task-resteyesc_desc-epochs_eeg.set` |
| `dortmund_dir` | Dortmund Vital Study resting-state EEG, OpenNeuro ds005385 (doi:10.18112/openneuro.ds005385.v1.0.2) | Figs. S7-S8 | `sub-XXX_ses-1_task-EyesClosed_acq-pre_eeg.edf`, `..._task-EyesOpen_...` |
| `electrodes_tsv` | electrode positions from the stroke dataset (BIDS `*_electrodes.tsv`) | Table S3 | columns `name` (or `channel`), `x`, `y`, `z` |

Subject lists are in `config/subjects.yaml`; file-name patterns can be changed there.

## Preprocessing

The stroke and PhysioNet recordings were preprocessed in EEGLAB before this analysis
(Appendix A2 of the paper); the file names above are the outputs of that step. The
OpenNeuro SRM epochs and the Dortmund recordings are used as distributed. See
`preprocessing/README.md`.

The stroke CSVs hold one column per channel with channel names in the header. The loader
removes the columns listed in `config/channels.yaml` (`stroke_drop_columns`) and checks
that exactly the 25 analysed channels remain.
