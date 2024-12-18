#Up-scaling: Code for the Publication
The importance of inland water CO2, CH4, and N2O for summertime greenhouse gas exchange with the atmosphere in Arctic tundra lowlands

##Overview
This repository contains the code used in the study "The importance of inland water CO2, CH4 and N2O for summertime greenhouse gas exchange with the atmosphere in Arctic tundra lowlands", currently under review at AGU Journal of Geophysical Research - Biogeosciences. The scripts and notebooks provided here replicate the data processing, analysis, and modeling used to produce the results in the manuscript.

##Data Availability
The input data used for this analysis is publicly available at:
https://doi.pangaea.de/10.1594/PANGAEA.939906

Please ensure that you have downloaded this dataset before running the scripts in this repository.

##Repository Structure
The repository includes the following files:

data_loading.py
Script to load and preprocess the input data. Use this script as the starting point to prepare the data for analysis.

LMM_siberia_concs.ipynb
Jupyter notebook for performing the linear mixed model (LMM) analysis described in the manuscript.

RF_classification.ipynb
Jupyter notebook for random forest classification of satellite images to estimate the surface areas of inland water systems.

flux_calculation_functions.py
Python script containing functions to calculate greenhouse gas fluxes based on the methods described in the manuscript.

up_scaling_2016.py and up_scaling_2017.py
Scripts to perform up-scaling of emissions by combining the calculated fluxes with surface area estimates for the years 2016 and 2017.

##Requirements
pip install -r requirements.txt

#Usage
Data Loading:
Run data_loading.py to load and preprocess the data. Ensure that the PANGAEA dataset is downloaded and available in the specified directory.

Analysis Steps:

Use LMM_siberia_concs.ipynb to perform linear mixed model analysis.
Use RF_classification.ipynb for random forest classification of satellite images.
Flux Calculation:
Call the functions in flux_calculation_functions.py to calculate greenhouse gas fluxes as described in the paper.

Emission Up-Scaling:
Execute up_scaling_2016.py and up_scaling_2017.py to calculate the up-scaled emissions for 2016 and 2017.

#License
This repository is licensed under the GNU General Public License v3.0. See the LICENSE file for details.

#Citation
If you use this code in your work, please cite the associated publication (citation to be updated once the paper is published).

