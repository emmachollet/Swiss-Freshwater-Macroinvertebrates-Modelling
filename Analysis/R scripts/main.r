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
CV <- T # Cross-Validation
extrapol <- ifelse(CV, FALSE, # Extrapolation
                  F
                  )
extrapol.info <- c(training.ratio = 0.8, variable = "IAR")

# Set if we allow data leakage (changes the way data is normalized)
dl <- F
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
run.ann <- T # Run ANN models or not (needs administrative rights)
analysis.dl <- F
analysis.ml <- F # Hyperparameter tuning (mainly for RF)
analysis.ann <- F # Hyperparameter tuning
analysis.training <- F

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
# if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") } # to run support vector machine (svm) algorithm
# if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
if ( !require("randomForest") ) { install.packages("randomForest"); library("randomForest") } # to run random forest (RF)
# if ( !require("RRF") ) { install.packages("RRF"); library("RRF") } # to run RF and additional features
if ( !require("xgboost") ) { install.packages("xgboost"); library("xgboost") } # to run Boosted Classification Trees
if ( !require("ada") ) { install.packages("ada"); library("ada") } # to run Boosted Classification Trees
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
file.inv.data     <- "occ_data_2020-06-25.dat"
file.prev         <- "prevalence_2020-06-25.dat"

# Load datasets
data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data),header=T,sep="\t", stringsAsFactors=T)
data.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.inv.data),header=T,sep="\t", stringsAsFactors=F)
prev.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.prev),header=T,sep="\t", stringsAsFactors=F)

## ---- DATA WRANGLING  ----

# Select env. factors ####

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
prev.inv <- prepro.data$prev.inv # dataframe prevalence after data pre-processing
list.taxa <- list.taxa[order(match(list.taxa, prev.inv$Occurrence.taxa))] # reorder taxa by prevalence
centered.data <- prepro.data$centered.data
centered.data.factors <- prepro.data$centered.data.factors
normalization.data <- prepro.data$normalization.data
remove(prepro.data)

# Select taxa ####

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

# Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c( 
  "deepskyblue4" = 'glm', # Generalized Linear Model
  "deepskyblue" = 'gamLoess',
  "#7B1359" = 'svmRadial', # Support Vector Machine
  # "darkmagenta" = 'RRF'#, # Regularized Random Forest
  "hotpink3" = 'rf' # Random Forest
  # "hotpink1" = "xgbTree" # Boosted Classification Tree
                )
no.algo <- length(list.algo)
                                       
## ---- APPLY MODELS ----

# Null model ####

# "Predicts" taxon prevalence everywhere
null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]

# Statistical models ####

ptm <- proc.time() # to calculate time of simulation
#file.name <- paste0(dir.models.output, "Stat_model_100iterations_corr_22taxa_CV_no_DL.rds") #to test

comm.corr.options <- c(T,F)
#comm.corr.options <- c(F)

names(comm.corr.options) <- c("CF0", "UF0")

stat.outputs <- mclapply(comm.corr.options, mc.cores = n.cores.stat.models, function(comm.corr){

  # comm.corr <- comm.corr.options[[1]]
  info.file.stat.name <- paste0("Stat_model_",
                                file.prefix,
                                no.taxa.full, "taxa_",
                                #126, "taxa_",
                                sampsize,"iterations_",
                                ifelse(comm.corr,"CF0_","UF0_"),
                                ifelse(CV, "CV_", "FIT_"),
                                ifelse(dl, "DL_", "no_DL_"))

  file.name <- paste0(dir.models.output, info.file.stat.name, ".rds")
  cat(file.name)

  # If the file with the output already exist, just read it
  if (file.exists(file.name) == T){

      if(exists("stat.output") == F){
          cat("\nFile with statistical model output already exists, we read it from", file.name, "and save it in object 'stat.output'\n")
          stat.output <- readRDS(file = file.name)

          }
      else{
          cat("\nList with statistical model output already exists as object 'stat.output' in this environment.")
      }
  } else {

      cat("\nNo statistical model output exist yet, we produce it and save it in", file.name)

      if(CV == T){

        stat.output <- mclapply(centered.data, mc.cores = n.cores.splits, FUN = stat_mod_cv, CV, comm.corr, sampsize, n.chain)

        cat("\nSaving output of statistical models in", file.name)
        saveRDS(stat.output, file = file.name, version = 2) #version two here is to ensure compatibility across R versions

      } else {
        # apply temporary on training data of Split1, later take whole dataset centered etc

        stat.output <- stat_mod_cv(data.splits = centered.data, CV, comm.corr, sampsize, n.chain = n.chain)
        cat("\nSaving output of statistical models in", file.name)
        saveRDS(stat.output, file = file.name, version = 2)
        }

  }
  return(stat.output)
})

print(paste("Simulation time of statistical model "))
print(proc.time()-ptm)

# Process output from stat models to fit structure of ml models (makes plotting easier)
stat.outputs.transformed <- transfrom.stat.outputs(CV, stat.outputs)

# Make vector statistical models
names(stat.outputs.transformed) <- c("chGLM", "hGLM")
list.stat.mod <- names(stat.outputs.transformed)
names(list.stat.mod) <- c("#59AB2D","#256801")

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
                             ifelse(dl, "DL_", "no_DL_"))

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
                                   FUN = apply.ml.model, list.algo, list.taxa, env.fact, prev.inv = prev.inv)
        } else {
            # Compute one split after the other
            ml.outputs.cv <- lapply(centered.data.factors, FUN = apply.ml.model, list.algo, list.taxa, env.fact, CV, prev.inv = prev.inv)
        }
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(ml.outputs.cv, file = file.name, version = 2)
    
        } else {
        
        splitted.data <- list("Training data" =  centered.data.factors[[1]], "Testing data" = data.frame())
        ml.outputs <- apply.ml.model(splitted.data = splitted.data, list.algo = list.algo, list.taxa = list.taxa,
                                      env.fact = env.fact, CV = F, prev.inv = prev.inv)
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(ml.outputs, file = file.name, version = 2)
        
        }
}

# Check if we have problem with some ml performance (>1.5)
check.outputs.perf(outputs.cv = ml.outputs.cv, list.taxa = list.taxa, CV = CV, extrapol = extrapol)

print(paste("Simulation time of different models ", info.file.ml.name))
print(proc.time()-ptm)

if(analysis.ml){
  if(CV){
    # Try to regularize RF ####
    
    # mtry.vect <- c(1,2,3,4,8)
    # mtry.vect <- c(1)
    # names(mtry.vect) <- paste("rf_", mtry.vect, "mtry", sep = "")
    
    list.tuned.grid <- list()
    tuned.grid <- expand.grid(mtry = c(1,2), coefReg = c(1), coefImp = c(0, 0.5))
    for (n in 1:nrow(tuned.grid)) {
      list.tuned.grid[[n]] <- tuned.grid[n,]
      names(list.tuned.grid)[[n]] <- paste(c("RRF_", paste(tuned.grid[n,], colnames(tuned.grid), sep="")), collapse = "")
    }
    list.tuned.algo <- names(list.tuned.grid)
    no.tuned.algo <- length(list.tuned.algo)
    names(list.tuned.algo) <- rainbow(no.tuned.algo)
    print(list.tuned.algo)
    
    tuned.ml.outputs.cv <- lapply(centered.data.factors, function(split){
      #split <- centered.data.factors[[1]]
      lapply(list.tuned.grid, FUN = apply.tuned.ml.model, splitted.data = split, algorithm = "RRF", list.taxa = list.taxa,
             env.fact = env.fact, CV = CV, prev.inv = prev.inv)})
  
    info.file.ml.name <-  paste0("ML_model_",
                                 file.prefix, 
                                 no.tuned.algo, "tunedRRF_",
                                 no.taxa, "taxa_", 
                                 ifelse(CV, "CV_", "FIT_"),
                                 ifelse(dl, "DL_", "no_DL_"))
    
    file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
    cat(file.name)
    saveRDS(tuned.ml.outputs.cv, file = file.name, version = 2)
  }
}


if(!server){

  
# Neural Networks ####

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
                              ifelse(CV, "CV_", "FIT_"),
                              ifelse(dl, "DL_", "no_DL_"))

file.name <- paste0(dir.models.output, info.file.ann.name, ".rds")
cat(file.name)

if( file.exists(file.name) == T ){
  cat("The file already exists. Reading it", file.name)
  if(CV){ 
    ann.outputs.cv <- readRDS(file = file.name)
  } else {
    ann.outputs <- readRDS(file = file.name) 
  }
} else {
    if(run.ann){
    
    cat("This ANN output doesn't exist yet, we produce it and save it in", file.name)
    if(CV){
      ann.outputs.cv <- lapply(centered.data, function(split){
        lapply(list.hyper.param, FUN = build_and_train_model, split = split,
               env.fact = env.fact,
               list.taxa = list.taxa,
               learning.rate = learning.rate,
               batch.size = batch.size,
               CV = CV)
      })
      saveRDS(ann.outputs.cv, file = file.name)
      
    } else {
      ann.outputs <- lapply(list.hyper.param, FUN = build_and_train_model, split = centered.data,
                            env.fact = env.fact,
                            list.taxa = list.taxa,
                            learning.rate = learning.rate,
                            batch.size = batch.size,
                            CV = CV)
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

# Merge outputs ####

# Change names of algorithms
print(list.algo)
list.algo.temp <- c( "iGLM", # Generalized Linear Model
                    "GAM", # Generalized Additive Model
                    "SVM", # Support Vector Machine
                    "RF" # Random Forest
)
names(list.algo.temp) <- names(list.algo) # keep the same colors
list.algo <- list.algo.temp
remove(list.algo.temp)
print(list.algo)
if(length(list.algo) != no.algo){ cat("The new list of ML algo dosen't match with the original one.") }
no.algo <- length(list.algo)

if(CV | extrapol){
  # Merge all CV outputs in one
  outputs.cv <- ml.outputs.cv
  # outputs.cv <- ann.outputs.cv1
  for (s in list.splits) {
    #s = "Split2"
    names(outputs.cv[[s]]) <- list.algo
    # names(ml.outputs.cv.bct[[s]]) <- "BCT"
    
    # outputs.cv[[s]] <- append(outputs.cv[[s]], tuned.ml.outputs.cv[[s]])
    # outputs.cv[[s]] <- append(outputs.cv[[s]], ml.outputs.cv.bct[[s]])
    
    
    if(exists("stat.outputs.transformed")){
        outputs.cv[[s]][[list.stat.mod[1]]] <- stat.outputs.transformed[[1]][[s]]
        outputs.cv[[s]][[list.stat.mod[2]]] <- stat.outputs.transformed[[2]][[s]]
    }
    
    if(exists("ann.outputs.cv")){
        # ECR: For ANN analysis
        # names(ann.outputs.cv2[[s]]) <- gsub("1FCT", "tanhFCT", names(ann.outputs.cv2[[s]]))
        # names(ann.outputs.cv2[[s]]) <- gsub("2FCT", "leakyreluFCT", names(ann.outputs.cv2[[s]]))
        names(ann.outputs.cv[[s]]) <- list.ann
    
        outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv[[s]])
    
        # outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv2[[s]])
        # outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv3[[s]])
        # outputs.cv[[s]] <- append(outputs.cv[[s]], tuned.ml.outputs.cv[[s]])
    }
  }
} else {
  # Make final outputs as list
  outputs <- append(append(ml.outputs, stat.outputs.transformed), ann.outputs)
  #outputs <- append(appendml.outputs, ann.outputs)
}
  
# ECR: For analysis
# To analyze ANN
# list.ann <- c(names(ann.outputs.cv[[1]]), names(ann.outputs.cv2[[1]]))
# no.ann <- length(list.ann)
# names(list.ann) <- rainbow(no.ann)
# To analyze ML or RRF
# outputs.cv <- tuned.ml.outputs.cv
# list.models <- c(list.algo, list.tuned.algo)
 
# Make final list of models

if( all(names(outputs.cv$Split1) %in% c(list.algo, list.stat.mod, list.ann))){
    list.models <- c(list.algo, list.stat.mod, list.ann)
} else {
    # list.models <- c(list.algo, list.stat.mod)
    # list.models <- names(outputs.cv$Split1)
    # list.models <- c(list.algo, list.ann)
    list.models <- c(list.algo,"hotpink" = "BRT", list.stat.mod, list.ann)
    list.models <- list.algo
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
                         ifelse(analysis.ml, "AnalysisML_", ""),
                         no.taxa, "taxa_", 
                         # no.env.fact, "envfact_",
                         ifelse(CV, "CV_", 
                                ifelse(extrapol, paste(c("extrapol", extrapol.info, "_"), collapse = ""), 
                                       "FIT_")),
                         ifelse(dl, "DL_", "no_DL_"),
                         "")
cat(info.file.name)

source("utilities.r")

# Produce final outputs with mean performance across splits
# ECR: I will fix CF0 soon ####
if(CV | extrapol){ 
    # Make final outputs as list
    # outputs <- make.final.outputs.cv(outputs.cv = outputs.cv, list.models = list.models, list.taxa = list.taxa)

    # Make final outputs as tables
    df.cv <- make.df.outputs(outputs = outputs.cv, list.models = list.models, 
                                 list.taxa = list.taxa, list.splits = list.splits,
                                 null.model = null.model, prev.inv = prev.inv, CV = CV, extrapol = extrapol)
    df.pred.perf.cv <- df.cv$`Table predictive performance CV`
    df.pred.perf <- df.cv$`Table predictive performance`
    df.fit.perf.cv <- df.cv$`Table fit performance CV`
    df.fit.perf <- df.cv$`Table fit performance`
    df.merged.perf <- df.cv$`Table merged`
    remove(df.cv)
} else {
    # Make final outputs as tables
    df.fit.perf <- make.df.outputs(outputs = outputs, list.models = list.models, 
                             list.taxa = list.taxa, list.splits = list.splits,
                             null.model = null.model, prev.inv = prev.inv, CV = CV)
}

## ---- PLOTS ----
source("utilities.r")
source("ml_model_functions.r")
source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# Models comparison ####

source("utilities.r")

# Table with performance

# Make csv file
file.name <- paste(dir.workspace, info.file.name,"TableResults",".csv", sep="")
write.table(df.merged.perf, file.name, sep=";", row.names=F, col.names=TRUE)

# HTML file with numbers
file.name <- paste0(info.file.name, "ModelsCompar_")

# tab.model.comp <- make.table(df.pred.perf = df.pred.perf, df.fit.perf = df.fit.perf, list.models = list.models)
# gtsave(data = tab.model.comp, filename = paste0(file.name, "Table_forAll.html"), path =  dir.plots.output)

tab.model.comp.species <- make.table.species(df.merged.perf = df.merged.perf, list.models = list.models)
gtsave(data = tab.model.comp.species, filename = paste0(file.name, "Table_perTaxonforAll.html"), path =  dir.plots.output)

tab.model.comp.species <- make.table.species.rearranged(df.merged.perf = df.merged.perf, list.models = list.models)
gtsave(data = tab.model.comp.species, filename = paste0(file.name, "Table_perTaxonperModel.html"), path =  dir.plots.output)

tab.model.comp.species <- make.table.species.rearranged.order(df.merged.perf = df.merged.perf, list.models = list.models)
gtsave(data = tab.model.comp.species, filename = paste0(file.name, "Table_perTaxonperModel_OrderedbyDiff.html"), path =  dir.plots.output)


# PDF file with colors
if(CV){
        # Tables prediction
        list.plots1 <- plot.df.perf(df.perf = df.pred.perf.cv, list.models = list.models, list.taxa = list.taxa, CV)
        list.plots2 <- plot.df.perf(df.perf = df.pred.perf, list.models = list.models, list.taxa = list.taxa, CV)
        
        temp.df.merged <- df.merged.perf[,which(grepl("likelihood.ratio", colnames(df.merged.perf)))]
        colnames(temp.df.merged) <- list.models
        temp.df.merged$Taxa <- df.merged.perf$Taxa
        list.plots3 <- plot.df.perf(df.perf = temp.df.merged, list.models = list.models, list.taxa = list.taxa, CV,
                                   title = "Comparison likelihood ratio")
        
        list.plots <- append(append(list.plots1, list.plots2), list.plots3)
        name <- "VisualTablesResults"
        file.name <- paste0(name, ".pdf")
        print.pdf.plots(list.plots = list.plots, width = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
        
        # if(dl){
        # list.plots.dl <- plot.dl.perf(df.pred.perf.dl.comb, list.models = list.models)
        # }
    } else {
        list.plots <- plot.df.perf(df.perf = df.fit.perf, list.models = list.models, list.taxa = list.taxa, CV)
}

temp.df.merged <- arrange(df.merged.perf, desc(CF0.pred.expl.pow.diff))
select.taxa <- temp.df.merged$Taxa[1:5]

# Performance against prevalence and boxplots

list.plots <- model.comparison(df.merged.perf = df.merged.perf, list.models = list.models, CV = CV, extrapol = extrapol, select.taxa = select.taxa)
name <- "PredFitModelsCompar_Boxplots"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 25, height = 14, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Performance training vs prediction

list.plots <- plot.perf.fitvspred(df.fit.perf = df.fit.perf, df.pred.perf = df.pred.perf, list.models = list.models, select.taxa = select.taxa)
name <- "PredFitModelsCompar_Scatterplot"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 12, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# # Print a jpeg file
# file.name <- paste0(name, ".jpg")
# jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
# print(list.plots[[1]])
# dev.off()

# Plots specifically related to trained models (and not to CV)

if(CV){
    
  # ECR: For ML analysis, if CV = T, take just the first split for trained models analysis (instead of running everything for CV=F, i.e. FIT, again)
  outputs <- outputs.cv[[1]]
  normalization.data.cv <- normalization.data
  normalization.data <- normalization.data.cv[[1]]

}  else {
  
    # ICE Manual ####
    
    source("plot_functions.r")
    
    list.list.plots <- lapply(select.taxa, FUN= plot.ice.per.taxa, outputs, list.algo = list.models[c(2:4,7)], env.fact, normalization.data, BDM)
    
    for (j in 1:length(select.taxa)) {
        taxon <- sub("Occurrence.", "", j)
        file.name <- paste0("ICE_", taxon, ".pdf")
        print.pdf.plots(list.plots = list.list.plots[[j]], width = 15, height = 10, 
                        dir.output = paste0(dir.plots.output, "ICE/"), 
                        info.file.name = info.file.name, file.name = file.name)
    }
    
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
    
    # Multiple (2) predictors PDP (of one model) for now just for 1 algo and 1 taxa
    list.plots <- plot.mult.pred.pdp(outputs = outputs, list.algo = list.algo,
                                     list.taxa = list.taxa, env.fact = env.fact)
    
    file.name <- "multpredPDP.pdf"
    print.pdf.plots(list.plots = list.plots, width = 17, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
    
    print(paste(file.name, "printing:"))
    print(proc.time()-ptm)

} # closing bracket of !CV clause

# Map predictions ####

ptm <- proc.time() # to calculate time of pdf production

inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

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

# Response shape ####

ptm <- proc.time() # to calculate time of pdf production

# select an algorithm to plot
algo <- list.algo[4]

# make a list with all plots and plot them in a pdf
list.plots <- lapply(list.taxa, FUN = response.ml.pred.taxa, outputs, list.algo, env.fact, algo, CV)

name <- paste0(algo, "_Resp_EnvFactvsTax")

# Print a pdf file
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()
print("Producing PDF time:")
print(proc.time()-ptm)

} # closing server braket

