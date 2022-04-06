## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- PhD Project ---
## 
##                          --- February 04, 2022 ---
##
## --- Emma Chollet, Jonas Wydler, Andreas Scheidegger and Nele Schuwirth ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# next two lines are for jonas, dont run these
# setwd("Q:/Abteilungsprojekte/siam/Jonas Wydler/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/R scripts")
#.libPaths("C:/Program Files/R/R-4.1.1/library")

## ---- PACKAGES, DATA & FCTS ----

# Check and set working directory
getwd() # show working directory
# setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts") # set the working directory to this folder

# Free workspace
rm(list=ls())
graphics.off()

# Setup options ####

# Set if we want to compute for the All or the BDM dataset
BDM <- F
file.prefix <- ifelse(BDM, "BDM_", "All_")

# Set date for file names
d <- Sys.Date()    # e.g. 2021-12-17

# Fit models to entire dataset or perform cross-validation (CV)
CV <- F # Cross-Validation
extrapol <- ifelse(CV, FALSE, # Extrapolation
                   F
)
extrapol.info <- c(training.ratio = 0.8, variable = "IAR")
dl <- F # Data Leakage
if(!CV){ dl <- F } # if it's only fitting, we don't need with or without dataleakage

# Set number of cores for Stat and ML models
n.cores.splits <-  3 # a core for each split, 3 in our case
n.cores.stat.models <- 1 # a core for each stat model 2 in our case (UF0, and CF0)
# Settings Stat models 
# Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models,
# also make sure the cross-validation (CV) flag is set correctly
sampsize <- 4000 #10000 # This needs to be an even number for some reason (stan error)
n.chain  <- 2 #2

# Select taxa
# TRUE: taxa with 5% < prev < 95%
# FALSE: taxa with intermediate prev (default 25% < prev < 75%, but can be changed)
all.taxa <- T

# Set analysis 
server <- F # Run the script on the server (and then use 3 cores for running in parallel)
run.ann <- F # Run ANN models or not (needs administrative rights)
analysis.dl <- F
analysis.ml <- F
analysis.ann <- F
analysis.training <- F

lme.temp <- F # Select if you want to use the linear mixed effect temperature model or not

# Load libraries ####

# Set a checkpoint to use same library versions, which makes the code repeatable over time
if ( !require("checkpoint") ) { install.packages("checkpoint"); library("checkpoint") }
checkpoint("2022-01-01") # replace with desired date
# checkpoint("2020-01-01", r_version="3.6.2") # replace with desired date and R version

if ( !require("parallel") ) { install.packages("parallel"); library("parallel") } # need to run things in parallel

# Data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data
if ( !require("splitTools") ) { install.packages("splitTools"); library("splitTools") } # to split the data

# Plots
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page
if ( !require("cowplot") ) { install.packages("cowplot"); library("cowplot") } # to arrange multiple plots on a page
if ( !require("pdp") ) { install.packages("pdp"); library("pdp") } # to plot partial dependance plots
if ( !require("gt") ) { install.packages("gt"); library("gt") } # to plot nice tables
if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") } # to plot nice tables
if ( !require("viridis")) {install.packages("viridis", repos="http://cloud.r-project.org"); library("viridis")} # to do even nicer plots
if ( !require("scales") ) { install.packages("scales"); library("scales") } # to look at colors
if ( !require("reshape2") ) { install.packages("reshape2"); library("reshape2") } # to reshape dataframes
if ( !require("gt") ) { install.packages("gt"); library("gt") } # to make tables

if(!server){ # packages having problems on the server
  if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map
  if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") } # to arrange multiple plots on a page
}

# Stat model
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") } # to read layers for map

# ANN model
if(run.ann){  # packages having problems with administrative rights
  library("reticulate")
  # install_miniconda() # run this the very first time reticulate is installed
  # 
  # install.packages("tensorflow")
  library("tensorflow")
  # install_tensorflow() # run this line only when opening R
  # 
  # install.packages("keras")
  library("keras")
  # install_keras() # run this line only when opening R
  # use_condaenv()
}

# ML models
if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") } # to run generalized additive model (GAM) algorithm
if ( !require("gam") ) { install.packages("gam"); library("gam") } # to run generalized additive model (GAM) algorithm
# if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") } # to run support vector machine (svm) algorithm
# if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
if ( !require("randomForest") ) { install.packages("randomForest"); library("randomForest") } # to run random forest (RF)
if ( !require("RRF") ) { install.packages("RRF"); library("RRF") } # to run RF and additional features
if ( !require("caret") ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models

# Load functions ####

source("ml_model_functions.r")
source("stat_model_functions.r")
source("plot_functions.r")
source("utilities.r")
if(run.ann){ source("ann_model_functions.r")}

# Load data ####

# Define directory and files
dir.env.data      <- "../../Data/Processed data/Environmental data/"
dir.inv.data      <- "../../Data/Processed data/Invertebrate data/"
dir.workspace     <- "../Intermediate results/"
dir.plots.output  <- "../Plots/Models analysis plots/"
dir.models.output <- "../Intermediate results/Trained models/"

file.env.data     <- "environmental_data_2020-06-25.dat"
file.env.data.lme <- "environmental_data_lme_2020-06-25.dat"
file.inv.data     <- "occ_data_2020-06-25.dat"
file.prev         <- "prevalence_2020-06-25.dat"

# Load datasets

data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data.lme),header=T,sep="\t", stringsAsFactors=T)
data.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.inv.data),header=T,sep="\t", stringsAsFactors=F)
prev.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.prev),header=T,sep="\t", stringsAsFactors=F)

## ---- DATA WRANGLING  ----

# Select env. factors ####
if (lme.temp)
{ env.fact <- c("temperature.lme",     # Temp
                "velocity",          # FV
                "A10m",              # A10m
                "cow.density",       # LUD
                "IAR",               # IAR
                "urban.area",        # Urban
                "FRI",               # FRI
                "bFRI",              # bFRI
                "width.variability") # WV

env.fact.full <- c(env.fact,
                   "temperature2.lme",
                   "velocity2")


}else
{ env.fact <- c("temperature",     # Temp
                "velocity",          # FV
                "A10m",              # A10m
                "cow.density",       # LUD
                "IAR",               # IAR
                "urban.area",        # Urban
                "FRI",               # FRI
                "bFRI",              # bFRI
                "width.variability") # WV

env.fact.full <- c(env.fact,
                   "temperature2",
                   "velocity2")
}

# env.fact <- env.fact.full
no.env.fact <- length(env.fact)

# Preprocess data ####

prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, prev.inv = prev.inv,
                               env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                               BDM = BDM, dl = dl, CV = CV, extrapol = extrapol, extrapol.info = extrapol.info)
data <- prepro.data$data
if(CV | extrapol){ splits <- prepro.data$splits
# rem.taxa <- prepro.data$rem.taxa # if we want to check if same taxa in each split
list.splits <- names(splits)
no.splits <- length(list.splits)
}

list.taxa <- prepro.data$list.taxa # list taxa after data pre-processing
list.taxa <- list.taxa[order(match(list.taxa, prev.inv$Occurrence.taxa))] # reorder taxa by prevalence
centered.data <- prepro.data$centered.data
centered.data.factors <- prepro.data$centered.data.factors
normalization.data <- prepro.data$normalization.data
prev.inv <- prepro.data$prev.inv
remove(prepro.data)

# Select taxa ####

cind.taxa <- which(grepl("Occurrence.",colnames(data)))
list.taxa.full <- colnames(data)[cind.taxa] # list of all taxa, before removing these too unbalanced or with too many NA
list.taxa.full <- list.taxa.full[order(match(list.taxa.full, prev.inv$Occurrence.taxa))] # reorder taxa by prevalence
remove(cind.taxa)

# 2 taxa
# list.taxa.int       <- list.taxa[list.taxa %in% c("Occurrence.Gammaridae", "Occurrence.Heptageniidae")]

# 5 taxa
# list.taxa.int       <- list.taxa[list.taxa %in% prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55),
#                            "Occurrence.taxa"]] # Select only few taxa

# Intermediate prevalence taxa
list.taxa.int     <- list.taxa[list.taxa %in% prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25),
                                                       "Occurrence.taxa"]] # Select with prevalence percentage between 25 and 75%

# Select taxa for prediction and analyses
if (all.taxa == F){
  list.taxa <- list.taxa.int
}

no.taxa.full <- length(list.taxa.full)
no.taxa.int <- length(list.taxa.int)
no.taxa <- length(list.taxa)

# Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c( 
  "deepskyblue4" = 'glm'#, # Generalized Linear Model
  #"deepskyblue" = 'gamLoess',
  #"#7B1359" = 'svmRadial', # Support Vector Machine
  # "darkmagenta" = 'RRF'#, # Regularized Random Forest
  #"hotpink3" = 'rf' # Random Forest
)
no.algo <- length(list.algo)


# Null model ####

null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]

# Machine Learning models ####

ptm <- proc.time() # to calculate time of simulation

info.file.ml.name <-  paste0("ML_model_",
                             file.prefix, 
                             no.algo, "algo_",
                             # ifelse(analysis.ml, "RFanalysis_", ""),
                             no.taxa, "taxa_", 
                             ifelse(CV, "CV_", 
                                    ifelse(extrapol, paste(c("extrapol", extrapol.info, "_"), collapse = ""), 
                                           "FIT_")),
                             ifelse(dl, "DL_", "no_DL_"),
                             ifelse(lme.temp, "lme.temp", "lm.tmp")
)

file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
# file.name <- paste0(dir.models.output, "ML_model_All_8tunedRRF_59taxa_CV_no_DL_.rds")
cat(file.name)

if( file.exists(file.name) == T ){
  
  cat("The file already exists. Reading it", file.name, "and uploading it in the environment.")
  if(CV){ ml.outputs.cv <- readRDS(file = file.name)
  } else { ml.outputs <- readRDS(file = file.name) }
  
} else {
  
  cat("No ML outputs exist yet, we produce it and save it in", file.name)
  if(CV == T | extrapol == T){
    
    if(server == T){
      # Compute three splits in paralel (should be run on the server)
      ml.outputs.cv <- mclapply(centered.data.factors, mc.cores = n.cores.splits,
                                FUN = apply.ml.model, list.algo, list.taxa, env.fact, env.fact.full, prev.inv = prev.inv)
    } else {
      # Compute one split after the other
      ml.outputs.cv <- lapply(centered.data.factors, FUN = apply.ml.model, list.algo, list.taxa, env.fact, env.fact.full, CV, prev.inv = prev.inv)
    }
    
    cat("Saving outputs of algorithms in", file.name)
    saveRDS(ml.outputs.cv, file = file.name, version = 2)
    
  } else {
    
    splitted.data <- list("Training data" =  centered.data.factors[[1]], "Testing data" = data.frame())
    ml.outputs <- apply.ml.model(splitted.data = splitted.data, list.algo = list.algo, list.taxa = list.taxa,
                                 env.fact = env.fact, env.fact.full, CV = F, prev.inv = prev.inv)
    
    cat("Saving outputs of algorithms in", file.name)
    saveRDS(ml.outputs, file = file.name, version = 2)
    
  }
}

names(ml.outputs)[[1]] <- paste("GLM_lm")
list.models <- list()
list.models <- append(list.models, ml.outputs)

