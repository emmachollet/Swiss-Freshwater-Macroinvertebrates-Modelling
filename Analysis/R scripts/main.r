## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- PhD Project ---
## 
##                          --- Endless moment, 2022 ---
##
## --- Emma Chollet, Jonas Wydler, Andreas Scheidegger and Nele Schuwirth ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ---- PACKAGES, DATA & FCTS ----

# Check and set working directory
getwd() # show working directory

# Free workspace
rm(list=ls())
graphics.off()

## Setup options ####

# Set if we want to compute for the All or the BDM dataset
BDM <- F
file.prefix <- ifelse(BDM, "BDM_", "All_")

# Set date for file names
d <- Sys.Date()    # e.g. 2021-12-17

# Fit models to entire dataset or perform cross-validation (CV) or out-of-domain generalization (ODG)
CV <- T # Cross-Validation
ODG <- ifelse(CV, FALSE, # If CV is TRUE, then no extrapolation
                  F # Set to TRUE for out-of-domain generalization
                  )
ODG.info <- c(training.ratio = 0.8, 
                   variable = "temperature")
                   # variable = "IAR")

# Set if we allow data leakage (changes the way data is normalized)
dl <- F
if(!CV){ dl <- F } # if it's only fitting, we don't need with or without dataleakage

# Set number of cores for Stat and ML models
n.cores.splits <-  1 # a core for each split, 3 in our case
n.cores.stat.models <- 1 # a core for each stat model 2 in our case (UF0, and CF0)

# Settings Stat models 
# Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models,
# also make sure the cross-validation (CV) flag is set correctly
sampsize <- 4000 #10000 # This needs to be an even number for some reason (stan error)
n.chain  <- 2 #2

# Select taxa
# TRUE: taxa with 5% < prev < 95%
# FALSE:  with intermediate prev (default 25% < prev < 75%, but can be changed)
all.taxa <- T

# Set analysis 
server <- F # Run the script on the server (and then use 3 cores for running in parallel)
run.ann <- T # Run ANN models or not (needs administrative rights)
analysis.dl <- F
analysis.ml <- F # Hyperparameter tuning or check randomness (mainly for RF) 
analysis.ann <- F # Hyperparameter tuning
analysis.training <- F
comparison.cv.odg <- F
comparison.lm.lme <- F
analysis.temp <- T

lme.temp <- T # Select if you want to use the linear mixed effect temperature model or not

## Load libraries ####

# Set a checkpoint to use same library versions, which makes the code repeatable over time
if ( !require("checkpoint") ) { install.packages("checkpoint"); library("checkpoint") }
checkpoint("2022-01-01") # replace with desired date
# checkpoint("2020-01-01", r_version="3.6.2") # replace with desired date and R version

if ( !require("parallel") ) { install.packages("parallel"); library("parallel") } # need to run things in parallel

# Data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data
if ( !require("splitTools") ) { install.packages("splitTools"); library("splitTools") } # to split the data
if ( !require("vtable") ) { install.packages("vtable"); library("vtable") } # to make table with summary statistics

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
if ( !require("DiagrammeR")) { install.packages("DiagrammeR"); library("DiagrammeR") } # to plot trees of BCT
if ( !require("skimr") ) { install.packages("skimr"); library("skimr") } # to show key descriptive stats

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
  # install_tensorflow() # run this line only when opening new R session
  # 
  # install.packages("keras")
  library("keras")
  # install_keras() # run this line only when opening new R session
  # use_condaenv()
}

# ML models
if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") } # to run generalized additive model (GAM) algorithm
if ( !require("gam") ) { install.packages("gam"); library("gam") } # to run generalized additive model (GAM) algorithm
## if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") } # to run support vector machine (svm) algorithm
## if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
if ( !require("randomForest") ) { install.packages("randomForest"); library("randomForest") } # to run random forest (RF)
## if ( !require("RRF") ) { install.packages("RRF"); library("RRF") } # to run RF and additional features

if(!server){ # packages having problems on the server
 if( !require("xgboost") ) { install.packages("xgboost"); library("xgboost") } # to run Boosted Classification Trees
 if( !require("ada") ) { install.packages("ada"); library("ada") } # to run Boosted Classification Trees
}
# have to be loaded at the end to not cache function 'train'
if ( !require("caret") ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models


## Load functions ####

source("ml_model_functions.r")
source("stat_model_functions.r")
source("plot_functions.r")
source("utilities.r")
if(run.ann){ source("ann_model_functions.r")}


## Load data ####

# Define directory and files
dir.env.data      <- "../../Data/Processed data/Environmental data/"
dir.inv.data      <- "../../Data/Processed data/Invertebrate data/"
dir.orig.data     <- "../../Data/Original data/"
dir.workspace     <- "../Intermediate results/"
dir.plots.output  <- "../Plots/Models analysis plots/"
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

if(lme.temp){
  data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data.lme),header=T,sep="\t", stringsAsFactors=T)
}else{
  data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data),header=T,sep="\t", stringsAsFactors=T)
}

data.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.inv.data),header=T,sep="\t", stringsAsFactors=F)
taxa.ibch         <- read.csv(paste0(dir.inv.data, file.ibch), header = T, sep = ";", stringsAsFactors = F)
prev.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.prev),header=T,sep="\t", stringsAsFactors=F)

# Prepare inputs for geographic plots
inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA WRANGLING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Select env. factors ####
env.fact <- c("temperature",     # Temp
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

no.env.fact <- length(env.fact)

## Preprocess data ####

# Remove NAs, normalize data, split data for CV or extrapolation
prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, prev.inv = prev.inv, 
                               remove.na.env.fact = T, remove.na.taxa = T, prev.restrict = 0.05,# 0.05,
                               env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                               BDM = BDM, dl = dl, CV = CV, ODG = ODG, ODG.info = ODG.info, lme.temp = lme.temp)
data <- prepro.data$data
if(CV | ODG){ splits <- prepro.data$splits
  # rem.taxa <- prepro.data$rem.taxa # if we want to check if same taxa in each split
  list.splits <- names(splits)
  no.splits <- length(list.splits)
}
list.taxa <- prepro.data$list.taxa # list taxa after data pre-processing
prev.inv <- prepro.data$prev.inv # dataframe prevalence after data pre-processing
standardized.data <- prepro.data$standardized.data
standardized.data.factors <- prepro.data$standardized.data.factors
normalization.data <- prepro.data$normalization.data
remove(prepro.data)

if(CV|ODG){
  
  # Plot analysis splits
  list.plots <- explore.splits(inputs = inputs, splits = splits, env.fact = env.fact)
  file.name <- "SplitsAnalysis"
  info.file.name <- paste0(file.prefix,
                           dim(data)[1], "samples_",
                           no.env.fact, "envfact_",
                           ifelse(CV, "CV_",
                                  ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                         "FIT_")),
                           ifelse(dl, "DL_", "no_DL_"),
                           ifelse(lme.temp, "lme_temp_", ""),
                           "")
  print.pdf.plots(list.plots = list.plots, width = 15, height = 8, dir.output = dir.expl.plots.output, info.file.name = info.file.name, file.name = file.name, png = TRUE)
}

## Select taxa ####

cind.taxa <- which(grepl("Occurrence.",colnames(data)))
list.taxa.full <- colnames(data)[cind.taxa] # list of all taxa
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

# Analyze which taxa needed to calculate the IBCH index
temp.list.taxa <- gsub("Occurrence.", "", list.taxa)
temp.list.taxa <- gsub("group.", "", temp.list.taxa)
temp.list.taxa <- sort(temp.list.taxa)

list.taxa.ibch <- taxa.ibch[,"Taxon_IBCH_2019..actuel."]
list.taxa.ibch <- gsub("Tachet", "", list.taxa.ibch)
list.taxa.ibch <- gsub(" ", "", list.taxa.ibch)
list.taxa.ibch <- gsub("[^[:alnum:]]", "", list.taxa.ibch)
list.taxa.ibch <- sort(list.taxa.ibch)

ind <- which(list.taxa.ibch %in% temp.list.taxa)
ind2 <- which(temp.list.taxa %in% list.taxa.ibch)
cat("Following", length(ind), "taxa needed for IBCH is present in our dataset:\n",
    list.taxa.ibch[ind],
    "\n\nFollowing", length(list.taxa.ibch) - length(ind),"taxa are missing in our dataset to calculate IBCH:\n",
    list.taxa.ibch[-ind],
    # "\n\nFollowing", length(ind2), "are in our dataset and needed for IBCH:\n",
    # temp.list.taxa[ind2],
    "\n\nFollowing", length(temp.list.taxa) - length(ind2), "are in our dataset and seem to not be needed for IBCH:\n",
    temp.list.taxa[-ind2]
)

## Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c( 
  "deepskyblue4" = 'glm', # Generalized Linear Model
  "deepskyblue" = 'gamLoess', # Generalized Additive Model
  "#7B1359" = 'svmRadial', # Support Vector Machine
  "hotpink1" = "ada" , # Boosted Classification Tree
  "hotpink3" = 'rf'  # Random Forest
                )
# # "darkmagenta" = 'RRF'#, # Regularized Random Forest

no.algo <- length(list.algo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
# APPLY MODELS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Null model ####

# "Predicts" taxon prevalence everywhere
null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]


## Statistical models ####

ptm <- proc.time() # to calculate time of simulation
# file.name <- paste0(dir.models.output, "Stat_model_100iterations_corr_22taxa_CV_no_DL.rds") #to test

comm.corr.options <- c(T,F)
# comm.corr.options <- c(F)

names(comm.corr.options) <- c("CF0", "UF0")

stat.outputs <- mclapply(comm.corr.options, mc.cores = n.cores.stat.models, function(comm.corr){

  # comm.corr <- comm.corr.options[[1]]
  info.file.stat.name <- paste0("Stat_model_",
                                file.prefix,
                                no.taxa.full, "taxa_",
                                sampsize,"iterations_",
                                ifelse(comm.corr,"CF0_","UF0_"),
                                ifelse(CV, "CV_",
                                       ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                              "FIT_")),
                                ifelse(dl, "DL_", "no_DL_"),
                                ifelse(lme.temp, "lme_temp",""))

  file.name <- paste0(dir.models.output, info.file.stat.name, ".rds")
  cat(file.name)

  # If the file with the output already exist, just read it
  if (file.exists(file.name)){

      if(!exists("stat.output")){
          cat("\nFile with statistical model output already exists, we read it from:\n", file.name, "\nand save it in object 'stat.output'.\n")
          stat.output <- readRDS(file = file.name)
          }
      else{
          cat("\nList with statistical model output already exists as object 'stat.output' in this environment.\n")
      }
  } else {

      cat("\nNo statistical model output exist yet, we produce it and save it in", file.name)

      if(CV|ODG){

        stat.output <- mclapply(standardized.data, mc.cores = n.cores.splits, FUN = stat_mod_cv, CV, ODG, comm.corr, sampsize, n.chain, list.taxa, env.fact.full)

        cat("\nSaving output of statistical models in:\n", file.name)
        saveRDS(stat.output, file = file.name, version = 2) #version two here is to ensure compatibility across R versions

      } else {
        # apply temporary on training data of Split1, later take whole dataset centered etc

        stat.output <- stat_mod_cv(standardized.data, CV, ODG, comm.corr, sampsize, n.chain, list.taxa, env.fact.full)
        cat("\nSaving output of statistical models in\n", file.name)
        saveRDS(stat.output, file = file.name, version = 2)
        }

  }
  return(stat.output)
})

print(paste("Simulation time of statistical model "))
print(proc.time()-ptm)


# if(!server){

# Process output from stat models to fit structure of ml models (makes plotting easier)
stat.outputs.transformed <- transfrom.stat.outputs(stat.outputs = stat.outputs, list.taxa = list.taxa, CV = CV, ODG = ODG)
stat.outputs.transformed <- stat.outputs.transformed[2:1]

# Make vector statistical models
names(stat.outputs.transformed) <- c("hGLM", "chGLM")
list.stat.mod <- names(stat.outputs.transformed)
names(list.stat.mod) <- c("#256801", "#59AB2D")



## Machine Learning models ####

ptm <- proc.time() # to calculate time of simulation

# Split computation of ML algorithm, because it' too heavy to read in
if(CV|ODG){ temp.list.algo <- list(list.algo[1:3], list.algo[4:5])
} else { temp.list.algo <- list(list.algo) }

temp.outputs <- list()

for (n in 1:length(temp.list.algo)) {
  # n = 1
  info.file.ml.name <-  paste0("ML_model_",
                               file.prefix,
                               paste0(paste(temp.list.algo[[n]], collapse = "_"), "_"),
                               # ifelse(analysis.ml, "RFanalysis_", ""),
                               no.taxa, "taxa_",
                               ifelse(CV, "CV_",
                                      ifelse(ODG, paste(c("ODG", ODG.info, "_"), collapse = ""),
                                             "FIT_")),
                               ifelse(dl, "DL_", "no_DL_"),
                               ifelse(lme.temp, "lme_temp", "")
                               )
  
  file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
  # file.name <- paste0(dir.models.output, "ML_model_All_8tunedRRF_59taxa_CV_no_DL_.rds")
  cat(file.name)
  
  if(file.exists(file.name)){
  
    cat("\nThe file already exists:\n", file.name, "\nReading it and uploading it in the environment.\n")
    
    if(CV | ODG){ temp.outputs[[n]] <- readRDS(file = file.name)
    } else { temp.outputs[[n]] <- readRDS(file = file.name) }
  
  } else {
  
    cat("\nNo ML outputs exist yet, we run the models and save the results in:\n", file.name)
    
    if(CV|ODG){
  
        if(server){
          # Compute three splits in parallel (should be run on the server)
          temp.outputs[[n]] <- mclapply(standardized.data.factors, mc.cores = n.cores.splits,
                                   FUN = apply.ml.model, temp.list.algo[[n]], list.taxa, env.fact, env.fact.full, CV, ODG, prev.inv = prev.inv)
        } else {
          # Compute one split after the other
          temp.outputs[[n]] <- lapply(standardized.data.factors, FUN = apply.ml.model, temp.list.algo[[n]], list.taxa, env.fact, env.fact.full, CV, ODG, prev.inv = prev.inv)
        }
  
        cat("\nSaving outputs of algorithms in:\n", file.name)
        saveRDS(temp.outputs[[n]], file = file.name, version = 2)
  
      } else {
  
        splitted.data <- list("Training data" =  standardized.data.factors[[1]], "Testing data" = data.frame())
        temp.outputs[[n]] <- apply.ml.model(splitted.data = splitted.data, list.algo = temp.list.algo[[n]], list.taxa = list.taxa,
                                      env.fact = env.fact, env.fact.full, CV = F, ODG, prev.inv = prev.inv)
  
        cat("\nSaving outputs of algorithms in:\n", file.name)
        saveRDS(temp.outputs[[n]], file = file.name, version = 2)
  
      }
  }
}

# Check if we have problem with some ml performance (>1.5)
# if(CV|ODG){check.outputs.perf(outputs.cv = ml.outputs.cv, list.taxa = list.taxa, CV = CV, ODG = ODG)}

print(paste("Simulation time of different models ", info.file.ml.name))
print(proc.time()-ptm)

source("ml_model_functions.r")

if(analysis.ml){
  if(CV|ODG){
    # # Try to regularize RF
    # 
    # # mtry.vect <- c(1,2,3,4,8)
    # # mtry.vect <- c(1)
    # # names(mtry.vect) <- paste("rf_", mtry.vect, "mtry", sep = "")
    # 
    # list.tuned.grid <- list()
    # tuned.grid <- expand.grid(mtry = c(1,2), coefReg = c(1), coefImp = c(0, 0.5))
    # for (n in 1:nrow(tuned.grid)) {
    #   list.tuned.grid[[n]] <- tuned.grid[n,]
    #   names(list.tuned.grid)[[n]] <- paste(c("RRF_", paste(tuned.grid[n,], colnames(tuned.grid), sep="")), collapse = "")
    # }
    # list.tuned.algo <- names(list.tuned.grid)
    # no.tuned.algo <- length(list.tuned.algo)
    # names(list.tuned.algo) <- rainbow(no.tuned.algo)
    # print(list.tuned.algo)
    # 
    # tuned.ml.outputs.cv <- lapply(standardized.data.factors, function(split){
    #   # split <- standardized.data.factors[[1]]
    #   lapply(list.tuned.grid, FUN = apply.tuned.ml.model, splitted.data = split, algorithm = "RRF", list.taxa = list.taxa,
    #          env.fact = env.fact, CV = CV, prev.inv = prev.inv)})

    # Check randomness of BCT and RF
    temp.list.algo <- c(
                        "ada" , # Boosted Classification Tree
                        "rf"  # Random Forest
                        )
    seeds <- c(
               2020,
               2021,
               2022,
               2023, 
               2024)
    
    grid.algo.seeds <- expand.grid("algorithm" = temp.list.algo, "seed" = seeds)
    grid.algo.seeds$algorithm <- as.character(grid.algo.seeds$algorithm)
    grid.algo.seeds$seed <- as.integer(grid.algo.seeds$seed)
    no.algo.seeds <- nrow(grid.algo.seeds)
    list.algo.seeds <- vector("list", no.algo.seeds)
    for (n in 1:no.algo.seeds) {
      list.algo.seeds[[n]] <- grid.algo.seeds[n,]
      names(list.algo.seeds)[n] <- paste(grid.algo.seeds[n,], collapse = "")
    }
    
    seeds.ml.outputs.cv <- lapply(standardized.data.factors, function(split){
      # split <- standardized.data.factors[[1]]
      lapply(list.algo.seeds, function(algo.seed){apply.tuned.ml.model(splitted.data = split, algorithm = algo.seed[,"algorithm"], list.taxa = list.taxa,
                                                                 env.fact = env.fact, CV = CV, ODG = ODG, prev.inv = prev.inv, seed = algo.seed[,"seed"])
      })
    })
    
    
    info.file.ml.name <-  paste0("ML_model_",
                                 file.prefix,
                                 no.algo.seeds, "AlgSeed",
                                 # no.tuned.algo, "tunedRRF_",
                                 no.taxa, "taxa_",
                                 ifelse(CV, "CV_", "FIT_"),
                                 ifelse(dl, "DL_", "no_DL_"),
                                 ifelse(lme.temp, "lme_temp", ""))

    # list.models <- names(list.algo.seeds)
    # no.models <- length(list.models)
    # names(list.models) <- rainbow(no.algo.seeds)
    # outputs.cv <- seeds.ml.outputs.cv
    list.algo <- names(list.algo.seeds)
    no.algo <- length(list.algo)
    names(list.algo)[which(grepl("ada", list.algo))] <- colorRampPalette(c("blue", "cadetblue1"))(5)
    names(list.algo)[which(grepl("rf", list.algo))] <- colorRampPalette(c("blueviolet", "hotpink3"))(5)
    ml.outputs.cv <- seeds.ml.outputs.cv
    
    file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
    cat(file.name)
    saveRDS(seeds.ml.outputs.cv, file = file.name, version = 2)
  }
}


# if(!server){

## Neural Networks ####

source("ann_model_functions.r")

# Set hyperparameters
learning.rate <- 0.01
batch.size <-  64
no.epo <- 50
act.fct <- c(#"tanh",
  # "swish",
  "leakyrelu") #,
  # "relu")

# ECR: For hyperparam grid search
# grid.hyperparam <- expand.grid(layers = c(3,5), units = c(32, 64), act.fct = act.fct, no.epo = c(150,200))
# grid.hyperparam <- expand.grid(layers = c(3,5), units = c(32, 64), act.fct = act.fct, no.epo = c(150,200))

# ECR: For specific hyperparam selection
grid.hyperparam <- expand.grid(layers = c(3), units = c(32), act.fct = act.fct, no.epo = no.epo)

grid.hyperparam$act.fct <- as.character(grid.hyperparam$act.fct)
no.hyperparam <- nrow(grid.hyperparam)
list.hyper.param <- vector("list", no.hyperparam)

if(no.hyperparam == 1){ # only one ANN selected

  if(analysis.ann){
    # Number of runs to check randomness in results
    no.ann.runs = 10
    for (n in 1:no.ann.runs) {
      list.hyper.param[[n]] <- grid.hyperparam[1,]
      names(list.hyper.param)[n] <- paste0("ANN_rand", n)
    }
    list.ann <- names(list.hyper.param)
    no.ann <- length(list.ann)
    names(list.ann) <- hcl.colors(no.ann, palette = "Oranges")
      # rainbow(no.ann)
  } else {
    list.hyper.param[[1]] <- grid.hyperparam[1,]
    names(list.hyper.param) <- c("ANN")
    list.ann <- names(list.hyper.param)
    no.ann <- length(list.ann)
    names(list.ann) <- "#FFB791" # if only one selected ANN
  }

} else { # try different hyperparameters

  for (n in 1:no.hyperparam) {
    list.hyper.param[[n]] <- grid.hyperparam[n,]
    names(list.hyper.param)[n] <- paste(paste0(grid.hyperparam[n,], c("L", "U", "FCT", "epo")), collapse = "")
  }

  names(list.hyper.param) <- paste("ANN_", names(list.hyper.param), sep = "")
  list.ann <- names(list.hyper.param)
  no.ann <- length(list.ann)
  names(list.ann) <- rainbow(no.ann) # assign colors
}

info.file.ann.name <-  paste0("ANN_model_",
                              file.prefix,
                              no.ann, "ann_",
                              ifelse(analysis.ann, "RandAnalysis_", ""),
                              no.taxa, "taxa_",
                              ifelse(CV, "CV_",
                                     ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                            "FIT_")),
                              ifelse(dl, "DL_", "no_DL_"),
                              ifelse(lme.temp, "lme_temp", ""))

file.name <- paste0(dir.models.output, info.file.ann.name, ".rds")
cat(file.name)

if(file.exists(file.name)){
  cat("The file already exists. Reading it", file.name)
  if(CV | ODG){
    ann.outputs.cv <- readRDS(file = file.name)
  } else {
    ann.outputs <- readRDS(file = file.name)
  }
} else {
    if(run.ann){

    cat("This ANN output doesn't exist yet, we produce it and save it in", file.name)
    if(CV|ODG){
      ann.outputs.cv <- lapply(standardized.data, function(split){
        lapply(list.hyper.param, FUN = build_and_train_model, split = split,
               env.fact = env.fact, list.taxa = list.taxa,
               learning.rate = learning.rate, batch.size = batch.size,
               CV = CV, ODG = ODG)
      })
      saveRDS(ann.outputs.cv, file = file.name)

    } else {
      ann.outputs <- lapply(list.hyper.param, FUN = build_and_train_model, split = standardized.data,
                            env.fact = env.fact, list.taxa = list.taxa,
                            learning.rate = learning.rate, batch.size = batch.size,
                            CV = CV, ODG = ODG)
      saveRDS(ann.outputs, file = file.name)
    }
    } else {
    cat("This ANN output doesn't exist yet and the run.ann condition is set to false.")
  }
}

# ECR: For ANN analysis, to compare more ann outputs
# file.name <- file.name <- paste0(dir.models.output, "ANN_model_All_1ann50epo_59taxa_CV_no_DL_.rds")
# ann.outputs.cv2 <- readRDS(file = file.name)
# ann.outputs.cv <- ann.outputs.cv

## Merge outputs ####

# Change names of algorithms but keeping colors selected at the beginning
print(list.algo)
list.algo.orig <- list.algo
list.algo <- c("iGLM", # Generalized Linear Model
               "GAM", # Generalized Additive Model
               "SVM", # Support Vector Machine
               "BCT", # Boosted Classification Tree
               "RF"   # Random Forest
               )

names(list.algo) <- names(list.algo.orig) # keep the same colors
print(list.algo)
if(length(list.algo) != no.algo){ cat("The new list of ML algo dosen't match with the original one.") }
if(length(temp.outputs) == 1){
  list.algo <- list.algo[1:length(temp.outputs[[1]])]
}
no.algo <- length(list.algo)

if(CV | ODG){
  # Merge all CV outputs in one
  outputs.cv <- vector(mode = "list", length = no.splits)
  names(outputs.cv) <- list.splits

  for (s in list.splits) {
    #s = "Split3"
    
    if(exists("temp.outputs")){
      outputs.cv[[s]][[list.algo[1]]] <- temp.outputs[[1]][[s]][[list.algo.orig[1]]]
    }
    
    if(exists("stat.outputs.transformed")){
        outputs.cv[[s]][[list.stat.mod[1]]] <- stat.outputs.transformed[[1]][[s]]
        outputs.cv[[s]][[list.stat.mod[2]]] <- stat.outputs.transformed[[2]][[s]]
    }

    if(exists("temp.outputs")){
      for (n in 1:length(temp.outputs)) {
        if(n == 1){
          for (l in 2:length(temp.outputs[[n]][[s]]) ) {
            outputs.cv[[s]][[list.algo[l]]] <- temp.outputs[[n]][[s]][[list.algo.orig[l]]]
          }
        } else {
          for (l in (length(temp.outputs[[n-1]][[s]])+1):(length(temp.outputs[[n-1]][[s]])+length(temp.outputs[[n]][[s]])) ) {
            outputs.cv[[s]][[list.algo[l]]] <- temp.outputs[[n]][[s]][[list.algo.orig[l]]]
          }
        }
      }
    }

    if(exists("ann.outputs.cv")){
      names(ann.outputs.cv[[s]]) <- list.ann
      outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv[[s]])

        # outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv2[[s]])
        # outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv3[[s]])
        # outputs.cv[[s]] <- append(outputs.cv[[s]], tuned.ml.outputs.cv[[s]])
    }

  }
  
} else {
  # Make final outputs as list
  outputs <- list()
  
  if(exists("temp.outputs")){
    outputs[[list.algo[1]]] <- temp.outputs[[1]][[1]]
  }
  if(exists("stat.outputs.transformed")){
    outputs[[list.stat.mod[1]]] <- stat.outputs.transformed[[1]]
    outputs[[list.stat.mod[2]]] <- stat.outputs.transformed[[2]]
  } 
  if(exists("temp.outputs")){
    for (l in list.algo[2:no.algo]) {
      outputs[[l]] <- temp.outputs[[1]][[l]]
    } 
  }
  if (exists("ann.outputs")){
    outputs <- append(outputs, ann.outputs)
  }
}

remove(temp.outputs)

# Make final list of models
if(CV|ODG){
  vec1 <- names(outputs.cv$Split1)
  vec2 <- c(list.algo, list.stat.mod, list.ann)
  list.models <- vec2[match(vec1, vec2)]
} else {
  vec1 <- names(outputs)
  vec2 <- c(list.algo, list.stat.mod, list.ann)
  list.models <- vec2[match(vec1, vec2)]
}
print(list.models)
no.models <- length(list.models)

# names(list.models) <- rainbow(no.models)

# See the final colors for each model
show_col(names(list.models))

info.file.name <- paste0(file.prefix,
                         no.models, "models_",
                         # no.models, "tunedRRF_",
                         ifelse(analysis.ann, "AnalysisANN_", ""),
                         ifelse(analysis.ml, "AnalysisRandomness_", ""),
                         no.taxa, "taxa_",
                         # no.env.fact, "envfact_",
                         ifelse(CV, "CV_",
                                ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                       "FIT_")),
                         ifelse(dl, "DL_", "no_DL_"),
                         ifelse(lme.temp, "lme_temp_", ""),
                         "")
cat(info.file.name)

source("utilities.r")

# Produce final outputs with mean performance across splits
# ECR: Have to fix CF0
if(CV | ODG){
  # Make final outputs as list
  # outputs <- make.final.outputs.cv(outputs.cv = outputs.cv, list.models = list.models, list.taxa = list.taxa)

  # Make final outputs as tables
  df.cv <- make.df.outputs(outputs = outputs.cv, 
                           list.models = list.models,
                           list.taxa = list.taxa, list.splits = list.splits,
                           null.model = null.model, prev.inv = prev.inv, CV = CV, ODG = ODG)
  df.pred.perf.cv <- df.cv$`Table predictive performance CV`
  df.pred.perf <- df.cv$`Table predictive performance`
  df.fit.perf.cv <- df.cv$`Table fit performance CV`
  df.fit.perf <- df.cv$`Table fit performance`
  df.perf <- df.cv$`Table merged`
  remove(df.cv)
    
} else {
  # Make final outputs as tables
  df.perf <- make.df.outputs(outputs = outputs, list.models = list.models,
                           list.taxa = list.taxa, list.splits = list.splits,
                           null.model = null.model, prev.inv = prev.inv, CV = CV, ODG = ODG)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("utilities.r")
source("ml_model_functions.r")
source("plot_functions.r")
# rm(list=ls())
# graphics.off()

## Cases comparison ####

# CV vs ODG comparison
if(comparison.cv.odg){
  
  names.appcase <- c("Cross-validation", "Generalization")
  
  if(CV){
    df.perf1 <- df.perf
  } else if(ODG){
    df.perf2 <- df.perf
  }
  
  CV <- !CV
  ODG <- !ODG
  
} else if(comparison.lm.lme){
  
  names.appcase <- c("lm", "lme")
  df.perf1 <- df.perf
  lme.temp <- !lme.temp
  
}

temp.info.file.name <- paste0(file.prefix,
                         no.models, "models_",
                         # no.models, "tunedRRF_",
                         ifelse(analysis.ann, "AnalysisANN_", ""),
                         ifelse(analysis.ml, "AnalysisRandomness_", ""),
                         no.taxa, "taxa_",
                         ifelse(CV, "CV_",
                                ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                       "FIT_")),
                         ifelse(dl, "DL_", "no_DL_"),
                         ifelse(lme.temp, "lme_temp_", ""),
                         "")
cat(temp.info.file.name)

if(comparison.cv.odg){
  
  if(CV){
    df.perf1 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableResults", ".csv", sep=""), header=TRUE, sep=";", stringsAsFactors=FALSE)
  } else if(ODG){
    df.perf2 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableResults", ".csv", sep=""), header=TRUE, sep=";", stringsAsFactors=FALSE)
  }
  
  CV <- !CV
  ODG <- !ODG
  
} else if(comparison.lm.lme) {
  
  df.perf2 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableResults", ".csv", sep=""), header=TRUE, sep=";", stringsAsFactors=FALSE)
  lme.temp <- !lme.temp
  
}

if(exists("df.perf1")){
  list.df.perf <- list(df.perf1, df.perf2)
  names(list.df.perf) <- names.appcase
}

# Select taxa for further analyses
if(exists("df.perf") & ("RF" %in% list.models)){
  # Reorder performance dataframe according to expl.pow of RF and select taxa with intermediate prevalence
  
  if(CV|ODG){
    temp.df <- df.perf[,which(grepl("expl.pow_", colnames(df.perf)) & grepl("pred", colnames(df.perf)))]
    colnames(temp.df) <- list.models
    temp.df$Taxa <- df.perf$Taxa
    temp.df <- arrange(temp.df, desc(RF))
  } else {
    temp.df <- df.perf
    temp.df <- temp.df[order(temp.df[,"RF"]),]
  }
  temp.df <- temp.df[which(temp.df$Taxa %in% list.taxa.int),]
  select.taxa <- c(temp.df$Taxa[1:3], "Occurrence.Psychodidae")
  cat("List of selected taxa:", sub("Occurrence.", "", select.taxa))
  # select.taxa <- "Occurrence.Gammaridae"
  # select.taxa <- "Occurrence.Psychodidae"
} else {
  select.taxa <- c("Occurrence.Gammaridae", "Occurrence.Psychodidae")
}

source("plot_functions.r")

if(exists("list.df.perf")){
  
  # Boxplots standardized deviance
  
  list.plots <- plot.boxplots.compar.appcase(list.df.perf = list.df.perf, list.models = list.models)
  file.name <- "ModelCompar_Boxplots"
  print.pdf.plots(list.plots = list.plots, width = 15, height = 8, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name, png = TRUE)
  
  # # Print a jpeg file
  # file.name <- paste0(name, ".png")
  # png(paste0(dir.plots.output,info.file.name,file.name), width = 1280, height = 720) # size in pixels (common 16:9 --> 1920???1080 or 1280???720)
  # print(list.plots[[1]])
  # dev.off()
  
  # Standardized deviance vs prevalence
  
  list.plots <- plot.perfvsprev.compar.appcase(list.df.perf = list.df.perf, list.models = list.models,
                                               list.taxa = list.taxa, select.taxa = select.taxa)
  file.name <- "ModelCompar_PerfVSPrev"
  print.pdf.plots(list.plots = list.plots, width = 15, height = 12, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name, png = TRUE)
  
  # # Print a jpeg file
  # file.name <- paste0(name, ".png")
  # png(paste0(dir.plots.output,info.file.name,file.name), width = 1080, height = 864) # size in pixels (common 5:4 --> 1080???864 )
  # print(list.plots[[1]])
  # dev.off()
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Models comparison ####

## GLM parameters comparison ####

list.glm <- list.models[1:3]

if(CV){
  list.plots <- plot.glm.param(outputs.cv, df.perf, list.glm, list.taxa, env.fact.full, CV)
} else {
  list.plots <- plot.glm.param(outputs, df.perf, list.glm, list.taxa, env.fact.full, CV)
}

file.name <- "GLMParametersComparison"
cat(file.name)
print.pdf.plots(list.plots = list.plots, width = 8.4, height = 11.88, # A4 proportions  
                dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)


## Performance table ####

# Make csv file with performance table
file.name <- paste(dir.workspace, info.file.name, "TableResults", ".csv", sep="")
cat(file.name)
write.table(df.perf, file.name, sep=";", row.names=F, col.names=TRUE)

# Make csv file with summary statistics
file.name <- paste(dir.workspace, info.file.name, "TableSummaryStat", ".csv", sep="")
cat(file.name)
df.summary <- sumtable(df.perf, out='return')
write.table(df.summary, file.name, sep=";", row.names=F, col.names=TRUE)

# Boxplots and Performance against prevalence

list.plots <- model.comparison(df.perf = df.perf, list.models = list.models, CV = CV, ODG = ODG, select.taxa = select.taxa)
name <- "PredFitModelsCompar_Boxplots"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 15, height = 6, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Performance training vs prediction

list.plots <- plot.perf.fitvspred(df.fit.perf = df.fit.perf, df.pred.perf = df.pred.perf, list.models = list.models, select.taxa = select.taxa)
file.name <- "PredFitModelsCompar_Scatterplot"
print.pdf.plots(list.plots = list.plots, width = 12, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# # Print a jpeg file
# file.name <- paste0(name, ".jpg")
# jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
# print(list.plots[[1]])
# dev.off()

# Plots specifically related to trained models (and not to CV)

if(CV | ODG){
  # ECR: For ML analysis, if CV = T, take just the first split for trained models analysis (instead of running everything for CV=F, i.e. FIT, again)
  outputs <- outputs.cv[[1]]
  normalization.data.cv <- normalization.data
  # normalization.data.ODG <- normalization.data
  normalization.data <- normalization.data.cv[[1]]
  
}


# produce table with main infos
df.summaries <- skim(data[,env.fact]) # produce dataframe containing descriptive statistics for each column
View(df.summaries)

mosaicplot(table(data$temperature, data$velocity),
           color = TRUE,
           xlab = "Species", # label for x-axis
           ylab = "Size" # label for y-axis
)

## Individual Cond. Exp. ####

source("plot_functions.r")
source("utilities.r")

no.samples <- 100
no.steps <- 200
subselect <- 1

data.orig <- data

subselect.taxa <- select.taxa
select.env.fact <- env.fact[c(1,2,5)]

list.list.plots <- lapply(subselect.taxa, FUN= plot.ice.per.taxa, outputs = outputs, ref.data = standardized.data[[1]], 
                          list.models = list.models, list.taxa = list.taxa,
                          env.fact = env.fact, select.env.fact = select.env.fact,
                          normalization.data = normalization.data, ODG = ODG, 
                          no.samples = no.samples, no.steps = no.steps, subselect = subselect)

for (j in 1:length(subselect.taxa)) {
  # j = 1
  for (k in select.env.fact) {
    # k = 1
    taxon <- sub("Occurrence.", "", subselect.taxa[j])
    file.name <- paste0("ICE_", no.samples, "samp_", taxon)
    print.pdf.plots(list.plots = list.list.plots[[j]][[k]], width = 8, height = 12, 
                    dir.output = paste0(dir.plots.output, "ICE/"), 
                    info.file.name = info.file.name, file.name = file.name)
  }
}

## Response shape ####

# list.list.plots <- lapply(select.taxa[1], FUN = plot.rs.taxa.per.factor, outputs, list.models, env.fact, CV, ODG)
list.list.plots <- lapply(select.taxa[1], FUN = plot.rs.taxa.per.model, outputs, list.models, env.fact, CV, ODG)


for (j in 1:length(select.taxa[1])) {
  taxon <- sub("Occurrence.", "", select.taxa[j])
  file.name <- paste0("RS_", taxon)
  print.pdf.plots(list.plots = list.list.plots[[j]], width = 8, height = 12,
                  dir.output = paste0(dir.plots.output, "ICE/"),
                  info.file.name = info.file.name, file.name = file.name)
}

# 10 samples - 200 steps - 5 taxa - 1 env fact --> take 2 minutes
# 100 samples - 200 steps - 5 taxa - 1 env fact --> take 5 minutes

# Recover/amalysis other metrics and information from ml outputs

# Variable importance
varImp(outputs[[6]][[select.taxa[1]]][["Trained model"]]) # use model specific method
varImp(outputs[[7]][[select.taxa[1]]][["Trained model"]], useModel = FALSE) # use model agnostic method, with AUC metric

# Other metrics (need to build dataframe with observation, predictions as levels and as probabilities)
test_set <- data.frame(pred = as.factor(outputs[[7]][[select.taxa[1]]][["Prediction factors training set"]]), # prediction as factors (present and absent)
                       obs = outputs[[7]][[select.taxa[1]]][["Observation training set"]][,select.taxa[1]], # observation as factors
                       present = outputs[[7]][[select.taxa[1]]][["Prediction probabilities training set"]]["present"], # predicted probability class 1 (present)
                       absent = outputs[[7]][[select.taxa[1]]][["Prediction probabilities training set"]]["absent"]) # predicted probability class 2 (absent)
confusionMatrix(data = test_set$pred, reference = as.factor(test_set$obs)) # calculate confusion matrix, accuracy, sensitivity, specificity
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall") # calculate precision, recall, F1
twoClassSummary(test_set, lev = levels(test_set$obs)) # calculate AUC, sensitivity, specificity
# bizarre

# file.name <- "testICE.pdf"
# pdf(paste0(dir.plots.output, info.file.name, file.name), paper = 'special', width = 20, # height = height,
#     onefile = TRUE)
# print(q)
# dev.off()

# PDP ####

# source("plot_functions.r")

ptm <- proc.time() # to calculate time of simulation

# PDP of one model
list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)


# PDP of one model
# list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
#                       list.taxa = list.taxa, env.fact = env.fact)
#
file.name <- "rf_PDP.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# PDP of all models
list.plots <- plot.pdp(outputs = outputs[[1]], list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "allPDP.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# PDP of all models
# We sub-select taxa and env.fact because it takes a lot of time
list.plots <- plot.pdp(outputs = outputs, list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "allPDP.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# ICE ####

ptm <- proc.time() # to calculate time of simulation

# ICE of one model
list.plots <- plot.ice(outputs = outputs, algo = list.algo[2], list.algo = list.algo,
                       list.taxa = list.taxa[1:2], env.fact = env.fact[1:2])

file.name <- "ICE.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# multiple PDP ####

# We just make one example because it's computationally heavy

ptm <- proc.time() # to calculate time of simulation

select.algo <- list.algo[4]

# Multiple (2) predictors PDP (of one model) for now just for 1 algo and 1 taxa
list.plots <- plot.mult.pred.pdp(outputs = outputs, list.algo = select.algo,
                                 list.taxa = select.taxa, env.fact = env.fact[c(1,5)])

file.name <- paste0("multpredPDP_", select.algo[1], "_", sub("Occurrence.", "", select.taxa[1]), ".pdf")
print.pdf.plots(list.plots = list.plots, width = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)



# Map predictions ####

ptm <- proc.time() # to calculate time of pdf production


# make a list with all plots and plot them in a pdf
list.plots <- lapply(list.taxa, FUN = map.ml.pred.taxa, inputs, outputs, list.algo, CV)

name <- "ObsvsPred_map"

# Print a pdf file
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()

print(paste(file.name, "printing:"))
print(proc.time()-ptm)



# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()
print("Producing PDF time:")
print(proc.time()-ptm)

source("plot_functions.r")

# Comparison temp models ####
if(analysis.temp & lme.temp){
  data.env.lme <- data
  
  # load dataset based on lm temperature model
  
  data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data),header=T,sep="\t", stringsAsFactors=T)
  
  data.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.inv.data),header=T,sep="\t", stringsAsFactors=F)
  prev.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.prev),header=T,sep="\t", stringsAsFactors=F)
  
  # load dataset containing information on the temperature measuring stations
  data.stations    <- read.table(paste(dir.orig.data,file.stations,sep=""),header=T,sep="\t", quote = "" )
  data.stations.used     <- readRDS(paste0(dir.orig.data,file.stations.used))
  data.stations.used.ID  <- unique(data.stations.used$ID)
  
  
# 
#   data.stations <- data.stations[data.stations$ID %in% data.stations.used.ID, ]
# 
#   # Extract X, Y values from coordinates
#   data.stations$X   <- as.numeric(sub("/.*", "", data.stations$coordinates))  
#   data.stations$Y   <- as.numeric(sub(".*/", "", data.stations$coordinates))
#   
#   # Some stations are at the wrong place (ID 2481 too far west, ID 2613 too far south)
#   data.stations[data.stations$ID == 2481,]$coordinates <- c("673551/202870")
#   data.stations[data.stations$ID == 2481,]$X <- 673551
#   
#   data.stations[data.stations$ID == 2613,]$coordinates <- c("611737, 272317")
#   data.stations[data.stations$ID == 2613,]$Y <- 272317 
  
  # load dataset containing information on the temperature measuring stations
  data.stations.used     <- readRDS(paste0(dir.orig.data,file.stations.used))
  data.stations.used.ID  <- unique(data.stations.used$ID)
  
  # Extract X, Y values from coordinates
  data.stations.used$X   <- as.numeric(sub("/.*", "", data.stations.used$coordinates))  
  data.stations.used$Y   <- as.numeric(sub(".*/", "", data.stations.used$coordinates))
  
  # Some stations are at the wrong place (ID 2481 too far west, ID 2613 too far south)
  data.stations.used[data.stations.used$ID == 2481,]$coordinates <- c("673551/202870")
  data.stations.used[data.stations.used$ID == 2481,]$X <- 673551
  #   
  data.stations.used[data.stations.used$ID == 2613,]$coordinates <- c("611737, 272317")
  data.stations.used[data.stations.used$ID == 2613,]$Y <- 272317 
  
  # Select env. factors 
  env.fact <- c("temperature",     # Temp
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
  
  
  # env.fact <- env.fact.full
  no.env.fact <- length(env.fact)
  
  
  # Preprocess data
  
  # Remove NAs, normalize data, split data for CV or extrapolation
  
  
  prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, prev.inv = prev.inv, 
                                 remove.na.env.fact = T, remove.na.taxa = T, prev.restrict = 0.05,# 0.05,
                                 env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                                 BDM = BDM, dl = dl, CV = CV, ODG = ODG, ODG.info = ODG.info, lme.temp = lme.temp)
  
  data.env.lm <- prepro.data$data
  
  data.diff <- data.env.lm[1:15]
  data.diff[5:15] <- data.env.lm[5:15] - data.env.lme[5:15]
  
  inputs1 <- map.inputs(dir.env.data = dir.env.data, data.env = data.env.lm)
  inputs2 <- map.inputs(dir.env.data = dir.env.data, data.env = data.env.lme)
  inputs3 <- map.inputs(dir.env.data = dir.env.data, data.env = data.diff)
  
  
  # Summary statistics
  summary(data.env.lm$temperature)
  summary(data.env.lme$temperature)
  summary(data.diff$temperature)
  
  
  pdf(file = paste0(dir.expl.plots.output,"Comparision_temp_models.pdf"))
  #
  map.env.fact.2(inputs1, env.fact, data.env.lm, info = "linear temperature model", data.stations = data.stations)
  map.env.fact.2(inputs2, env.fact, data.env.lme, info = "linear mixed effect temperature model", data.stations = data.stations)
  map.env.fact.2(inputs3, env.fact, data.diff, info = "Difference between temperature models", data.stations = data.stations)
  dev.off()
 
  sink(file = paste0(dir.expl.plots.output, "summary_statistics_temperature_models.txt"))
  print(paste0("Temperature based on linear model:"))
  print(summary(data.env.lm$temperature))
  print(paste0("Temperature based on linear mixed effects model:"))
  print(summary(data.env.lme$temperature))
  print(paste0("Difference between the two temperature models:"))
  print(summary(data.diff$temperature))
  sink()
  
  
}


# } # closing server braket

