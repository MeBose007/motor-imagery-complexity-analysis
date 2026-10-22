# Preprocessing (EEGLAB)

The preprocessing script will be added here. The steps (Appendix A2 of the paper):

1. Import EDF; assign electrode locations from a standard cap layout.
2. Band-pass 1-50 Hz, zero-phase FIR.
3. `clean_rawdata` bad-channel detection: flat-line > 5 s, channel correlation 0.8,
   line-noise criterion 4, burst criterion 20; high-pass and window rejection disabled.
4. Extended Infomax ICA; ICLabel; reject components with probability >= 0.9 of eye or
   muscle activity.
5. Re-reference to the common average.
6. Interpolate removed channels (spherical spline).
7. Save as EDF / SET (stroke recordings were then exported to CSV, one column per channel).
