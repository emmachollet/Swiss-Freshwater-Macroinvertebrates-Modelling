Archive for PhD Project

# "Bridging the gap between data science and mechanistic modelling to gain knowledge about macroinvertebrate community assembly"

Author: Emma Chollet Ramampiandra (ecr)
Supervisor: Nele Ina Schuwirth (nis)
Workflow first author: Peter Vermeiren (pv)
GIS data support: Rosi Siber (rs)
Data science support: Andreas Scheidegger (as)

-----------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------

## Part I : Data processing

-----------------------------------------------------------------------------------------------------------------------------------

We need to process  A. invertebrate data and B. environmental data, and C. explore final datasets (to find missing values, problems)

Directories:
R scripts: *Swiss Freshwater Macroinvertebrates Modelling/Data/R scripts/*
Original data/input files: *Swiss Freshwater Macroinvertebrates Modelling/Data/Original data/*
Invertebrate input/output files: *Swiss Freshwater Macroinvertebrates Modelling/Data/Processed data/Invertebrate data/*
Environmental input/output files: *Swiss Freshwater Macroinvertebrates Modelling/Data/Processed data/Environmental data/*
Output plots: *Swiss Freshwater Macroinvertebrates Modelling/Analysis/Plots/Explorative plots/*

-----------------------------------------------------------------------------------------------------------------------------------

### A. Process invertebrate data

#### Step 1
Add SiteId and SampleId, correct taxonomic spelling, separate abundance and abundance class, construct taxon name, 
aggregate duplicated entries (due to "species complexes") ,
add Longitude and Latitude, produce wide format, remove useless columns

Input:
	- *MIDAT-data-export-20200625_ALL_2_transm_2.csv*   # provided by the CSCF (contact: Maxime Chèvre)

Script:
	- **convert_invertebrates_data_2020-06-25.r**

Utilities:
	- *utilities.r*

Outputs:
	- *summary_processed_WideInvData _ 2020-06-25 .txt*
	- *invertebrates_2020-06-25.dat*
	- *invertebrates_wide_2020-06-25.dat*
	- *invertebrates_taxonomy_2020-06-25.dat*  
	- *invertebrates_wide_occ_homogenizedtaxonomy_2020-06-25.dat*
	- *invertebrates_metadata_2020-06-25.txt*

#### Step 2
Add Sampling.window, remove samples out of official sampling window, harmonize taxonomic spelling (with CH and EU database),
identify mixed taxonomic levels (ALL and BDM), chose taxa level with more than 5% of prevalence, produce ALL and BDM datasets,
produce information file for rs, plot  prevalence for ALL and BDM datasetts

Inputs:	
	- *invertebrates_wide_2020-06-25.dat*    # output of I.A.1
	- *invertebrates_taxonomy_2020-06-25.dat*    # output of I.A.1
	- *taxa_CH_FreshEco.csv*    #  invertebrate taxa list based on freshwater ecology database
	- *taxa_EU_FreshEco.csv*    #  invertebrate taxa list based on freshwater ecology database
	
Script:	
	- **get_tax_dataset_2020-06-25.r**

Utilities:	
	- *utilities.r*
	- *get_tax_level.r*
	
Outputs:
	- *SitesData_for_RS_2020-06-25.dat*    # only sites-samples information table to get environmental data from rs
	- *ALL_occ_data_2020-06-25.dat*
	- *BDM_occ_data_2020-06-25.dat*
	- *BDM.fam_occ_data_2020-06-25.dat*
	- *mixed_taxonomy_BDM_2020-06-25 .txt*
	- *mixed_taxonomy_ALL_2020-06-25 .txt*
	- *mixed_taxonomy_BDM_fam_2020-06-25 .txt*
	- *All_prev_perc_2020-06-25.dat*
	- *BDM_prev_perc_2020-06-25.dat*
	- two plots of prevalence percentage (ALL and BDM)

**Scripts checked 29 June 2021**

-----------------------------------------------------------------------------------------------------------------------------------

### B. Process environmental data

#### Step 1
Define substrate types and covering classes, convert BDM substrate data in substrate and substrate fraction coverage 

Input:	
	- *BDM_EPT_Kopfdaten_vonChristianStickelberger_20200811.csv*

Script:	
	- **convert_habitatdata_BDM_2020-08-11.r**

Outputs:
	- *BDM_Habitat_samplesubstratefraction_2020-08-11.dat*
	- *BDM_Habitat_substrate_coverage_2020-08-11.dat*

#### Step 2
Convert BDM ecomorpholgy data column names (e.g. EINDOL, GSBREITE, ...) and contents in numerical values

Inputs:	
	- *BDM_EPT_Kopfdaten_vonChristianStickelberger_20200811.csv*
	- *invertebrates_wide_2020-06-25.dat*
	- *SitesData_for_RS_2020-06-25_result.csv* # provided by rs based on *SitesData_for_RS_2020-06-25.dat* output of A. Step 2 
				(contains FRI and bFRI factors previously provided by Bogdan Caradima)

Script:	
	- **convert_ecomorph_BDM_2020-08-11.r**

Output:	
	- *BDM_EcoMorph_data_HintermannWeber_2020-08-11.dat*

Note: ecomorph and substrate data for NAWA SPEZ and TREND available for 2018 and 2019 respectively (see folders below in Original data), 
but is not worth to be investigated (won't be compared to BDM)
- SPEZ18_für Nele.zip    # seperated excel files for each site
- TREND19_Nele.zip    # seperated excel files for ecomorph for each site 
- TREND19_Nele.zip/für Nele/GesamtDatentabelle_MZB_20210318.xlsx    # contains invertebrate and substrate data

#### Step 3
Add substrate and ecomorphology data, calculate environmental factor (e.g. temperature, discharge, ...), 
combine data from CH and InS, calculate values of morphology msk between 0 and 1, produce env dataset (while removing useless factors)

Inputs:	
	- *SitesData_for_RS_2020-06-25_result.dat*    # provided by rs based on *SitesData_for_RS_2020-06-25.dat* (output of A. Step 2)
	- *BDM_Habitat_samplesubstratefraction_2020-08-11.dat*    # output of I.B.1
	- *BDM_Habitat_substrate_coverage_2020-08-11.dat*    # output of I.B.1
	- *BDM_EcoMorph_data_HintermannWeber_2020-08-11.dat*    # output of I.B.2
	- *invertebrates_wide_2020-06-25.dat*    # output of I.A.1
	- *ranking_env_data.csv*    # produced manually with colnames of old environmental dataset and ordered by nis and ecr in levels 0, 1, 2 and 3, 
				if colnames don't match need to be updated with the env dataset produced by the new workflow

Script:	
	- **convert_environmental_data_2021-05-26.r**

Outputs:
	- *All_environmental_data_2021-05-26.dat*
	- *BDM_environmental_data_2020-06-25.dat*

**Scripts checked 30 July 2021**

-----------------------------------------------------------------------------------------------------------------------------------

### C. Explore data
Classify environmental factors by their importance from prior knowledge and select the ones to be used in priority, compute and plot dataframe with main summary statistics, plot correlation matrix, plot factors on Swiss map, plot boxplots for distribution of factors per taxa presence/absence

Inputs:	
	- *All_occ_data_2020-06-25.dat*
	- *BDM_occ_data_2020-06-25.dat*
	- *All_environmental_data_2020-06-25.dat*
	- *BDM_environmental_data_2020-06-25.dat*

	- *ranking_env_data.csv*
	- *environmentaldata_documentation_20210728.csv*
	- 

Script:
	- **explorative_plots.r**

Outputs:



**Scripts checked 01 July 2021**

-----------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------

## Part II : Create explorative plots

-----------------------------------------------------------------------------------------------------------------------------------

We need to produce A. explorative plots of environmental factors in Switzerland 
and outputs of different models.

Directories:
R script: *Swiss Freshwater Macroinvertebrates Modelling/Analysis/*
Invertebrate intput files: *Swiss Freshwater Macroinvertebrates Modelling/Data/Processed data/Invertebrate data/*
Environmental and Swiss map intput files: *Swiss Freshwater Macroinvertebrates Modelling/Data/Processed data/Environmental data/*
Output files: *Swiss Freshwater Macroinvertebrates Modelling/Analysis/Plots/*

-----------------------------------------------------------------------------------------------------------------------------------

### A. Producing plots

#### Step 1
To be written ...