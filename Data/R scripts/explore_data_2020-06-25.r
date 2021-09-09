## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Explore Env data (if missing values, distribution, ...) ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 30, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Import libraries and data ----

# if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } 
# !! Should maybe be written like above and not like below ####

# Load packages
library(dplyr) # to sort, join, merge data
library(ggplot2) # to do nice plots
library(skimr) # to show key descriptive stats

# Check and set working directory
getwd() # show working directory
setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Data/R scripts") # set the working directory to this folder

# Free workspace
rm(list=ls())
graphics.off()

# Define directory and files
dir.env.data      <- "../Processed data/Environmental data/"
dir.inv.data      <- "../Processed data/Invertebrate data/"
# dir.output        <- "../Processed data/"
dir.output        <- "../../Analysis/Plots/"

file.inv.data     <- "All_occ_data_2020-06-25.dat"
file.inv.BDM.data <- "BDM_occ_data_2020-06-25.dat"
file.env.data     <- "All_environmental_data_2020-06-25.dat"
file.env.BDM.data <- "BDM_environmental_data_2020-06-25.dat"

file.prev         <- "All_prev_perc_2020-06-25.dat"
file.prev.BDM     <- "BDM_prev_perc_2020-06-25.dat"

# Read data
# BDM data
data.inv.BDM      <- read.delim(paste(dir.inv.data,file.inv.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
prev.inv.BDM      <- read.delim(paste(dir.inv.data,file.prev.BDM,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
data.env.BDM      <- read.delim(paste(dir.env.data,file.env.BDM.data,sep=""), na = c("<Null>", "NA"), # directly replace <Null> by Na
                                header=T, sep="\t", stringsAsFactors=FALSE)
file.prefix.BDM   <- "BDM_"

# All data
data.inv          <- read.delim(paste(dir.inv.data,file.inv.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
prev.inv          <- read.delim(paste(dir.inv.data,file.prev,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
data.env          <- read.delim(paste(dir.env.data,file.env.data,sep=""), na = c("<Null>", "NA"),  # directly replace <Null> by Na
                                header=T, sep="\t", stringsAsFactors=FALSE)
file.prefix       <- "All_"


## ---- Format data sets ----

## ---- Replace Null by NA ----

# # for null in number columns (that are character vectors for now because of the Null)
# data.env$mqn.Jahr
# as.numeric(data.env$mqn.Jahr)
# 
# # for null in factor columns
# skimmed <- skim(data.env) # produce dataframe containing descriptive statistics for each column
# data.env$mqn.Jahr[data.env$mqn.Jahr=="<Null>"] <- NA
# 
# # EASIER, do this everytime we read it
# read.delim(blabla, na.strings= c("<Null>", "NA", blabla) )

# Convert character columns to factors
data.env[sapply(data.env, is.character)] <- lapply(data.env[sapply(data.env, is.character)], 
                                                   as.factor)
# Construct BDM data set
data.BDM <- data.env %>%
    filter(MonitoringProgram == "BDM") # filter the rows, for columns use select
    # data[which(data.env$MonitoringProgram =="BDM") , ]
dim(data.BDM)

# BDM SampId and row indices for additional checks
BDM.samples <- data.env[which(data.env$MonitoringProgram =="BDM") , "SampId"]
length(BDM.samples) # should have 886 samples
BDMind <- which(data.env$SampId %in% BDM.samples)


## ---- Explore data ----


# Produce table with main infos
skimmed <- skim(data.env) # produce dataframe containing descriptive statistics for each column
skimmed.BDM <- skim(data.BDM)

# Based on observation in 'skimmed' explore data

# Missing River Typology info
head(data.env[which(is.na(data.env$HOEHE)),])
data.env[which(is.na(data.env$HOEHE)), "Canton"]

data.BDM[which(is.na(data.BDM$HOEHE)),]
# CSCF_698283_2014-03-07 BDM sample doesn't have HOEHE ABFLUSS GEFAELLE  GEO

summary(data.env$NAME.FGTy)

# Lots of Null values, that are missing values (should be NA)
# as.numeric(as.character(data.env$mqn.Jahr))

# Missing Slope
# data.env[which(is.na(data.env$Slope)), "SampId"]
data.env[which(is.na(data.env$Slope)), "Canton"] # lots of Hors frontières

# catchment area (and mean elevation)
data.env[which(is.na(data.env$A.EZG)), "Canton"] # because are cantons on the frontières ?
data.env[which(is.na(data.env$DEM.EZG)), "Canton"]

# Problematic because then missing temperature, cow density
# Same with all .ANT
# because then missing IAR
# But no missing values for BDM

# Discharge
hist(data.env$Discharge)
summary(data.env$Discharge)

data.env[which(is.na(data.env$Discharge)), "Canton"]

data.BDM[which(is.na(data.BDM$Discharge)), "SampId"]
# Problematic for these sites because then Discharge used for other things (ara.fract -> sapr cond, velocity)
# Espcially for 3 BDM sites:
# [1] CSCF_764139_2011-06-26 CSCF_638123_2014-07-02 CSCF_764139_2016-06-22

data.BDM[which(is.na(data.BDM$Discharge)), "Canton"]
data.BDM[which(data.BDM$SampId == "CSCF_764139_2011-06-26"), "Discharge"]

# EINDOL
data.BDM[which(is.na(data.BDM$EINDOL)), "SampId"]
# [1] CSCF_758147_2013-06-16
# problematic for morphology etc
