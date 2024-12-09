# Macroalgae detritus decomposition and cross-shelf carbon export from shallow and deep reefs
This repository contains data and scripts relevant to *Limnology & Oceanography* article doi: XXXX. Ocean model and particle tracking simulation data are freely available and can be requested from MvdM (mirjam.vandermheen@uwa.edu.au). Description of data columns are in the sections below.

# Decomposition #
#### '2021 Marmion Kelp Litter Data - Raw.csv' - decomposition measurements (biomass and chemical composition) in 2 study species
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
#### '2021 Marmion Kelp Decomp.R' - Code to analyse and visualize biomass and carbon and nitrogen decomposition in 2 study species
- Input = '2021 Marmion Kelp Litter Data - Raw.csv' from Decomposition
- Output = Figures 2-4 

# Particle tracking #
#### 'export estimates 0-60 days.csv' - macroalgae detritus age and deep ocean export estimates (200 m depth +)
- Particle age (days) = detrital age (days since detachment)
- source reef depth = categorical variable with levels 10 - 20 m and 50 m representing depth (m) of source reef below sea level (lowest astronomical tide)
- percentage of particles past shelf break (200m) = continuous variable; predicted percentage of detritus exported to the deep ocean (see data availability in description above)
#### '2021 Marmion export model.R' - Code to analyse and visualize biomass and carbon and nitrogen decomposition in 2 study species
- Input = 'export estimates 0-60 days.csv' from Particle tracking
- Output = Figure 5


