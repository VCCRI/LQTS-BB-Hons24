Files used for 2024 Honours Project - Exploring the effect of beta-blockers on repolarisation parameters in Long QT Syndrome

There are two folders currently:
  Batch_Ishne_Conversion:
    This folder contains all critical functions to convert ishne .ecg and .ann file formats to Arash formats for processing with the LQTS Holter Analysis folder. Developed by AS, refined by GP. BATCH_CONVERSION is the main script to run. Can handle multiple files if the folder is the selected path.
  New_LQTS_Analysis_Program-Hons-24 folder is a bit of a misnomer. It can be used to analyse any Arash format ecg file, not just LQTS. BATCH_ARASH_FORMAT_ANALYSIS is the main script to run.


Magic numbers in SSIM_recs using the repol savers. Pre-assigned the repol_savers size which is problematic if it changes. TODO: Fix such that summary_data and repol_savers are not magic numbers and then modify SSIM_recs to accept new modifications.

Magic numbers also used in the assingment of nbins for pixel size of heatmaps - which are then used in difference_maps. Easy to fix, by assigning them as parameters in LQTS_Porgram_function instead (save it into ecg_analysis to pull it out of the function for difference_maps later).