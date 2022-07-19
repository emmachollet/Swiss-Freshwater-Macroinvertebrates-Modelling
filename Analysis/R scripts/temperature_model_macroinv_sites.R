# Investigate differences of two temperature models (linear and mixed linear effects) temperature for sites 
# of the "Bridging the gap between data science and mechanistic modelling to gain knowledge about community assembly" 
# project of Emma Chollet Ramampiandra and Nele Schuwirth.

# Author: Jonas Wydler, 17.07.2022

# libraries
if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots

# Set if we want to compute for the All or the BDM dataset
BDM <- F
file.prefix <- ifelse(BDM, "BDM_", "All_")

# Load data ####

# Define directory and files
dir.env.data      <- "../../Data/Processed data/Environmental data/"
dir.inv.data      <- "../../Data/Processed data/Invertebrate data/"
dir.orig.data     <- "../../Data/Original data/"
dir.workspace     <- "../Intermediate results/"
dir.expl.plots.output <- "../Plots/Explorative plots/"
dir.models.output <- "../Intermediate results/Trained models/"

file.env.data     <- "environmental_data_2020-06-25.dat"
file.env.data.lme <- "environmental_data_lme_2020-06-25.dat"
file.inv.data     <- "occ_data_2020-06-25.dat"
file.ibch         <- "ListTaxaIBCH2019.csv"
file.prev         <- "prevalence_2020-06-25.dat"
file.stations     <- "temperature_stations.dat"
file.stations.used     <- "temperature_dataset_11_07.rds"
# Load datasets
data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data),header=T,sep="\t", stringsAsFactors=T)

# Load functions ####
source("plot_functions.r")
source("utilities.r")


  
# load dataset containing information on the temperature measuring stations
data.stations.used     <- readRDS(paste0(dir.orig.data,file.stations.used))
data.stations.used.ID  <- unique(data.stations.used$ID)
  
# Extract X, Y values from coordinates
data.stations.used$X   <- as.numeric(sub("/.*", "", data.stations.used$coordinates))  
data.stations.used$Y   <- as.numeric(sub(".*/", "", data.stations.used$coordinates))
  
# Some stations are at the wrong place (ID 2481 too far west, ID 2613 too far south)
data.stations.used[data.stations.used$ID == 2481,]$coordinates <- c("673551/202870")
data.stations.used[data.stations.used$ID == 2481,]$X <- 673551
 
data.stations.used[data.stations.used$ID == 2613,]$coordinates <- c("611737, 272317")
data.stations.used[data.stations.used$ID == 2613,]$Y <- 272317 
  
# Prepare inputs for geographic plots
inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)
data.env$temperature.models.diff <- data.env$temperature.lme - data.env$temperature.lm

  
pdf(file = paste0(dir.expl.plots.output,"Comparision_temp_models.pdf"))
 
  temperature.plots(inputs, data.env, info = "based on lme temperature model", variable = "temperature.lme")
  
  temperature.plots(inputs, data.env, info = "based on lm temperature model", variable =  "temperature.lm")
  
  temperature.plots(inputs, data.env, info = "difference between lme and lm", variable = "temperature.models.diff")
  
dev.off()
  
sink(file = paste0(dir.expl.plots.output, "summary_statistics_temperature_models.txt"))
 
  print(paste0("Temperature based on linear model:"))
  print(summary(data.env$temperature.lm))
  print(paste0("Temperature based on linear mixed effects model:"))
  print(summary(data.env$temperature.lme))
  print(paste0("Difference between the two temperature models:"))
  print(summary(data.env$temperature_models.diff))
  
sink()
  



