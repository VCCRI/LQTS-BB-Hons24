Files used for 2024 Honours Project - Exploring the effect of beta-blockers on repolarisation parameters in Long QT Syndrome

There are two folders currently:
  Batch_Ishne_Conversion:
    This folder contains all critical functions to convert ishne .ecg and .ann file formats to Arash formats for processing with the LQTS Holter Analysis folder. Developed by AS, refined by GP. BATCH_CONVERSION is the main script to run. Can handle multiple files if the folder is the selected path.
  New_LQTS_Analysis_Program-Hons-2024
    This folder is a bit of a misnomer. It can be used to analyse any Arash format ecg file, not just LQTS. BATCH_ARASH_FORMAT_ANALYSIS is the main script to run. analyse_ecg extracts parameters and LQTS_Program_Function perform tasks like heatmaps (LQTS_Program_Function).