%% [Folders] overview

[Jeff_scripts] - All scripts shared by Jeff Wang.

[Plot.R] - R scripts used for plotting the R01 prelim figures

[Preprocess.m] - Matlab scripts used for R01 prelim iEEG data preprocessing. Most adapted from Jeff's scripts.

[Preprocess.xlsx] - Spreadsheets used as queue file in the preprocessing pipeline (see notes for .m scripts below).

[Statistics.R] - R scripts for prelim statistical analyses.


%% Scripts in [Preprocess.m]

LL_Process_NLX_Files.m/LL_Process_NLX_Files_esTT.m - Uses LL_Process_Queue.xlsx in [Preprocess.xlsx] to pull nlx files and do filtering, artifact rejection, etc. and convert everything to MATLAB formats.

LL_Fix_Trigs.m - splits triggers into sham vs TMS trials, does channel rejection, etc. Uses LL_FixTrigs_Queue.xlsx in [Preprocess.xlsx] to queue up files.

LL_epoch_trials.m - does detrending and epochs the time trace based on trigger times. Uses LL_Epoch_Queue.xlsx in [Preprocess.xlsx].

LL_epoch_stats.m - nonparametric clustering to determine significant iTEPs. Uses xl_Epoch_Queue.xlsx in [Preprocess.xlsx].

xl_get_iTEPtable.m/xl_get_CCEPtable.m - Export preprocessed iEEG data into .csv format for statistical analyses.