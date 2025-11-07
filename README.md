# Hg_SOM

This repository contains Octave/MATLAB scripts for training and visualizing Self-Organizing Maps (SOMs) applied to Rosati's mercury (Hg) and methylmercury (MeHg) observation dataset.  
The workflow loads, normalizes, and models environmental and biogeochemical data to explore patterns between MeHg concentration, temperature, salinity, and oxygen.

# Usage

To fit, evaluate, and project the species distribution models, open and run **Hg_SOM.m**.

## Data Source 

The data and functions required to run **Hg_SOM.m** are organized as follows:

Cossarini_Model_Output/

│ 20200101_y-CMCC--PSAL-MFSe3r1-MED-b20220901_re-sv01.00.nc

│ 20200101_y-CMCC--TEMP-MFSe3r1-MED-b20220901_re-sv01.00.nc

│ 20200101_y-OGS--BIOL-MedBFM3-MED-b20221120_re-sv05.00.nc

Rosati_Data_Model_Output/

│ ave.20140616-00_00_00.DMHg.nc

│ ave.20140616-00_00_00.MMHg.nc

│ data_ancillary_DGM.csv

│ data_ancillary_HgT.csv

│ data_ancillary_MeHg.csv

somtoolbox/

## Software used
Octave-10.3.0


