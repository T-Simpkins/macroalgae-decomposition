# Macroalgae detritus decomposition and cross-shelf carbon export from shallow and deep reefs
This repository contains data and scripts relevant to *Limnology & Oceanography* article doi: XXXX. Description of data columns are in the following sections below.

# Decomposition #
- '2021 Marmion Kelp Litter Data - Raw.csv'# - decomposition measurements (biomass and chemical composition) in 2 study species
   - Sample_ID = naming convention for each sample: 'site_timepoint_replicate.first letter of Genus'
   - Bag_no = unique litterbag identifier; each litterbag contained 7 replicates of each species
   - species = categorical variable with 2 levels: *Ecklonia radiata*, and *Scytothalia dorycarpa*
   - site = categorical variable with levels SM5 (collection site), SHALLOW 1, SHALLOW 2, SHALLOW 3, MID 1, MID 2, MID 3, DEEP 1, DEEP 2, and DEEP 3
   - depth = categorical variable with levels 10, 20, and 50 representing depth (m) below sea level (lowest astronomical tide)
   - time = detrital age (days)
   - time_int = sampling interval as a categorical variable with levels T0, T1, T2, T3, and T4
   - ww = wet weight (g)
   - dw = dry weight (g)
   - rb = remaining biomass as a percentage of the initial biomass
   - d15N
   - d13C
   - N
   - C
   - C.N
- '2021 Marmion Kelp Decomp.R'# - Code to analyse and visualize biomass and carbon and nitrogen decomposition in 2 study species
  - Input = '2021 Marmion Kelp Litter Data - Raw.csv' from ##Decomposition##
  - Output = Figures 2-4 
  


