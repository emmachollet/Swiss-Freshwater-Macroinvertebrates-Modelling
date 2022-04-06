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
 # "deepskyblue4" = 'glm',#, # Generalized Linear Model
  #"deepskyblue" = 'gamLoess',
  #"#7B1359" = 'svmRadial', # Support Vector Machine
  # "darkmagenta" = 'RRF'#, # Regularized Random Forest
  "hotpink3" = 'rf' # Random Forest
                )
no.algo <- length(list.algo)
                                       
## ---- APPLY MODELS ----

# Null model ####

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
list.stat.mod <- names(stat.outputs.transformed)
names(list.stat.mod) <- c("#59AB2D","#256801")

#res2 <- stat.outputs[[1]][[1]]
#res.extracted   <- rstan::extract(res2,permuted=TRUE,inc_warmup=FALSE)

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
                             ifelse(lme.temp, "lme_temp_", "lm_temp_")
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

# Check if we have problem with some ml performance, which would suggest that the model didn't converge
if(CV | extrapol){
  cat("Check if performance > 1.5 in\n")
    for (s in list.splits) {
        for (l in list.algo) {
            cat( s, "for model", l, "\n")
            list.taxa.temp <- names(ml.outputs.cv[[s]][[l]])
            if(length(list.taxa.temp) != no.taxa){
              cat(length(list.taxa.temp), "taxa for this algorithm\n")
            }
            for (j in list.taxa.temp) {
                    
                    perf.test <- ml.outputs.cv[[s]][[l]][[j]][["Performance testing set"]]
                    perf.train <- ml.outputs.cv[[s]][[l]][[j]][["Performance training set"]]
                    
                    if(perf.train > 1.5){
                      cat(j, "has training performance bigger than 1.5", perf.train, "\n")
                    }
                    
                    if(perf.test > 1.5){
                      cat(j, "has testing performance bigger than 1.5", perf.test, "\n")
                    }
                    
                    # if(length(ml.outputs.cv[[s]][[l]][[j]][["Trained model"]]) == 1){
                    #   ml.outputs.cv[[s]][[l]][[j]][["Trained model"]]
                    #   cat("This model has NULL MODEL instead of trained algorithm\n")
                    # }
                    # to replace performance having a problem
                    # ml.outputs.cv[[s]][[l]][[j]][["Performance testing set"]] <- ifelse(perf.test > 1.5, Inf, perf.test)
                    # ml.outputs.cv[[s]][[l]][[j]][["Performance training set"]] <- ifelse(perf.train > 1.5, Inf, perf.train)
                }
        }
    }
} else {
    for (l in list.algo) {
        print(l)
        list.taxa.temp <- names(ml.outputs[[l]])
        for (j in list.taxa.temp) {
            perf <- ml.outputs[[l]][[j]][["Performance training set"]]
            if(perf > 1.5){
              cat(j, "has training performance bigger than 1.5", perf, "\n")
            }
            # ml.outputs[[l]][[j]][["Performance training set"]] <- ifelse(perf > 1.5, Inf, perf)
        }
    }
}

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
list.algo.temp <- c( "GLM", # Generalized Linear Model
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
    # names(outputs.cv[[s]]) <- list.algo
    
    outputs.cv[[s]] <- append(outputs.cv[[s]], tuned.ml.outputs.cv[[s]])
    
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

# JW analyis for different temp models, need to move somewhere else

outputs <- list()
outputs <- append(outputs, ml.outputs)
names(outputs)[[1]] <- "GlM_LM"
names(outputs)[[2]] <- "RF_LM"
#saveRDS(outputs, file = paste0(dir.models.output,"temp_comp_glm"), version = 2)
temp.lm <- readRDS(file = paste0(dir.models.output, info.file.ml.name,".rds"))
names(temp.lm) <- c("glm_lm_temp", "rf_lm_temp")

temp.lme <- readRDS(file = paste0(dir.models.output, info.file.ml.name,".rds"))
names(temp.lme) <- c("glm_lme_temp", "rf_lme_temp")

names(outputs)
outputs <- append(outputs, temp.lme[1:2])
outputs <- append(outputs, temp.lm[1:2])

outputs <- outputs[c("glm_lm_temp", "glm_lme_temp", "rf_lm_temp", "rf_lme_temp")]
#saveRDS(outputs, file = paste0(dir.models.output,"temp_comp_2models.rds"), version = 2)


list.models <- names(outputs)
no.models <- length(list.models)
# ECR: For analysis
# To analyze ANN
# list.ann <- c(names(ann.outputs.cv[[1]]), names(ann.outputs.cv2[[1]]))
# no.ann <- length(list.ann)
# names(list.ann) <- rainbow(no.ann)
# To analyze ML or RRF
# outputs.cv <- tuned.ml.outputs.cv
# list.models <- c(list.algo, list.tuned.algo)
 
# Make final list of models

if(exists("ann.outputs.cv")){
 list.models <- c(list.algo, list.stat.mod, list.ann)
} else {
    # list.models <- c(list.algo, list.stat.mod)
    # list.models <- names(outputs.cv$Split1)
    # list.models <- c(list.algo, list.ann)
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
                             null.model = null.model, prev.inv = prev.inv, CV = CV, extrapol = extrapol)
}
#saveRDS(df.fit.perf, file = paste0(dir.models.output,"rf_glm_lm.rds"), version = 2)

# Look into the effect of data leakage ####

# save intermediate output to compare dl no dl
# saveRDS(df.pred.perf, file = paste0(dir.models.output, info.file.name, "df_perf_.rds"), version = 2)
#df.pred.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_7models_59taxa_CV_no_DL_df_perf_.rds"))
#df.pred.perf.dl <- readRDS(file = paste0(dir.models.output, "All_7models_59taxa_CV_DL_df_perf_.rds"))

# df.pred.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_no_DL_df_perf_.rds"))
# df.pred.perf.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_DL_df_perf_.rds"))
#
# df.pred.perf.no.dl$DL <- F
# df.pred.perf.dl$DL <- T
# df.pred.perf.dl.comb <- rbind(df.pred.perf.no.dl,df.pred.perf.dl)
#
# #remove infite vlaues to take mean later
# df.pred.perf.dl.comb$rf[!is.finite(df.pred.perf.dl.comb$rf)] <- NA
# df.pred.perf.dl.comb$gamSpline[!is.finite(df.pred.perf.dl.comb$gamSpline)] <- NA
# df.pred.perf.dl.comb$glm[!is.finite(df.pred.perf.dl.comb$glm)] <- NA
#
# colnames(df.pred.perf.dl.comb)
#
# df.pred.perf.dl.comb.table <- as.data.frame(df.pred.perf.dl.comb %>% group_by(DL) %>%
#                                       summarise(glm.mean = mean(glm, na.rm = T),
#                                       glm.sd = sd(glm, na.rm = T),
# 
#                                       gamLoess.mean = mean(gamLoess,na.rm = T),
#                                       gamLoess.sd = sd(gamLoess,na.rm = T),
# 
#                                       svmRadial.mean = mean(svmRadial, na.rm = T),
#                                       svmRadial.sd = sd(svmRadial, na.rm = T),
# 
#                                       rf.mean = mean(rf, na.rm = T),
#                                       rf.sd = sd(rf, na.rm = T),
# 
#                                       CF0.mean = mean(CF0, na.rm = T),
#                                       CF0.sd = sd(CF0, na.rm = T),
# 
#                                       UF0.mean = mean(UF0, na.rm = T),
#                                       UF0.sd = sd(UF0, na.rm = T),
# 
# 
#                                       ANN_3L32UleakyreluFCT50epo.mean = mean(ANN_3L32UleakyreluFCT50epo, na.rm = T),
#                                       ANN_3L32UleakyreluFCT50epo.sd = sd(ANN_3L32UleakyreluFCT50epo, na.rm = T)
# 
#                                       #ANN3L64U.mean = mean(ANN3L64U, na.rm = T),
#                                       #ANN3L64U.sd = sd(ANN3L64U, na.rm = T),
# 
#                                       #ANN5L32U.mean = mean(ANN5L32U, na.rm = T),
#                                       #ANN5L32U.sd = sd(ANN5L32U, na.rm = T)
# 
#                                        ))

#saveRDS(df.pred.perf.dl.comb.table, file = paste0(dir.plots.output, "Table_means_dl.rds"), version = 2)
#df.pred.perf.dl.comb.table <- readRDS(file = paste0(dir.plots.output, "Table_means_dl.rds"))
# diff.means.dl <- df.pred.perf.dl.comb.table[1,c("glm.mean", "gamLoess.mean","svmRadial.mean",
#                          "rf.mean", "CF0.mean", "UF0.mean", "ANN_3L32UleakyreluFCT50epo.mean")] -
#   df.pred.perf.dl.comb.table[2,c("glm.mean", "gamLoess.mean","svmRadial.mean",
#                             "rf.mean", "CF0.mean", "UF0.mean", "ANN_3L32UleakyreluFCT50epo.mean")]
# 
# diff.means.dl <- df.pred.perf.dl.comb.table[1,c("glm.mean", "CF0.mean", "UF0.mean")] -
#   df.pred.perf.dl.comb.table[2,c("glm.mean", "CF0.mean", "UF0.mean")]
# 
# 
# 
# dloui <- subset(df.pred.perf.dl.comb, DL == T)
# dloui <- dloui[,1:7]
# dlno <- subset(df.pred.perf.dl.comb, DL == F)
# dlno <- dlno[,1:7]
# comb <- cbind(dlno,dloui)
# comb2 <- dloui[,1:7] - dlno[,1:7]
# comb2 <- cbind(comb2,dlno[,8])
# colnames(comb2)[8] <- "Taxa"
# 
# comb3 <- comb2
# mean_difference <- colMeans(comb2[1:7])
# sd_difference <- colSds(as.matrix(comb2[1:7]))
# names(sd_difference) <- c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN")
# 
# table2 <- comb3 %>% gt()
# 
# tmp.table <- comb2[c("Taxa", "glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN_3L32UleakyreluFCT50epo")]
# colnames(tmp.table) <- c("Taxa", "glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN")
# #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
# table1 <- tmp.table %>% gt() %>%
#   tab_header(
#     title = md("**Difference in predictive performance between dl and no dl**") # make bold title
#   ) %>%
#   fmt_number(
#     columns = c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN"), # round numbers
#     decimals = 3
#   ) %>% # remove uneccessary black lines
#   tab_options(
#     table.border.top.color = "white",
#     heading.border.bottom.color = "black",
#     row_group.border.top.color = "black",
#     row_group.border.bottom.color = "white",
#     #stub.border.color = "transparent",
#     table.border.bottom.color = "white",
#     column_labels.border.top.color = "black",
#     column_labels.border.bottom.color = "black",
#     table_body.border.bottom.color = "black",
#     table_body.hlines.color = "white")
# 
# table2 <- comb3 %>% gt() %>%
#   tab_header(
#     title = md("**Difference in predictive performance between dl and no dl**") # make bold title
#   ) %>%
#   fmt_number(
#     columns = c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN_3L32UleakyreluFCT50epo"), # round numbers
#     decimals = 3
#   ) %>% # remove uneccessary black lines
#   tab_options(
#     table.border.top.color = "white",
#     heading.border.bottom.color = "black",
#     row_group.border.top.color = "black",
#     row_group.border.bottom.color = "white",
#     #stub.border.color = "transparent",
#     table.border.bottom.color = "white",
#     column_labels.border.top.color = "black",
#     column_labels.border.bottom.color = "black",
#     table_body.border.bottom.color = "black",
#     table_body.hlines.color = "white")
# 


## ---- PLOTS ----
source("utilities.r")
source("ml_model_functions.r")
source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# # Stat model traceplots ####
# res <- stat.outputs[[1]][[1]][[1]]
# res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
# 
# print(traceplot(res,pars=c(names(res)[134:162],"lp__")))
# print(traceplot(res))

# Models comparison ####

source("utilities.r") # JW: Why do we source the utilities all the time

# Table with performance
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
        list.plots.cv <- plot.df.perf(df.perf = df.pred.perf.cv, list.models = list.models, list.taxa = list.taxa, CV)
        list.plots <- plot.df.perf(df.perf = df.pred.perf, list.models = list.models, list.taxa = list.taxa, CV)
        list.plots1 <- append(list.plots.cv, list.plots)

        name <- "TablesPredPerf"
        file.name <- paste0(name, ".pdf")
        print.pdf.plots(list.plots = list.plots1, width = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
        
        temp.df.merged <- df.merged.perf[,which(grepl("likelihood.ratio", colnames(df.merged.perf)))]
        colnames(temp.df.merged) <- list.models
        temp.df.merged$Taxa <- df.merged.perf$Taxa
        list.plots <- plot.df.perf(df.perf = temp.df.merged, list.models = list.models, list.taxa = list.taxa, CV,
                                   title = "Comparison likelihood ratio")
        name <- "TableLikelihoodRatiof"
        file.name <- paste0(name, ".pdf")
        print.pdf.plots(list.plots = list.plots, width = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
        
        # Tables fit
        # list.plots.cv <- plot.df.perf(df.perf = df.fit.perf.cv, list.models = list.models, list.taxa = list.taxa, CV = F)
        # list.plots <- plot.df.perf(df.perf = df.fit.perf, list.models = list.models, list.taxa = list.taxa, CV = F)
        # list.plots <- append(list.plots.cv, list.plots)

        # if(dl){
        # list.plots.dl <- plot.dl.perf(df.pred.perf.dl.comb, list.models = list.models)
        # }
    } else {
        list.plots <- plot.df.perf(df.perf = df.fit.perf, list.models = list.models, list.taxa = list.taxa, CV)
}
# 
# name <- "TablesFitPerf"
# file.name <- paste0(name, ".pdf")
# # list.plots.dl[[1]]
# print.pdf.plots(list.plots = list.plots, width = 12, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
#print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = paste0(info.file.name,"_temp_models"), file.name = file.name)


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

if(!CV){
    
  # ECR: For ML analysis, if CV = T, take just the first split for trained models analysis (instead of running everything for CV=F, i.e. FIT, again)
  outputs <- outputs.cv[[1]]
  normalization.data.cv <- normalization.data
  normalization.data <- normalization.data[[1]]
  
    # Performance vs hyperparameters ####
    
    # Compute plots
    list.plots <- plot.perf.hyperparam(outputs = outputs, 
                                       list.algo = list.algo[4], # GLM algo doesn't have hyperparameters
                                       list.taxa = list.taxa)
    # Print a pdf file
    name <- "PerfvsHyperparam"
    file.name <- paste0(name, ".pdf")
    if( file.exists(paste0(dir.plots.output, info.file.name, file.name)) == F ){ # produce pdf only if it doesn't exist yet (takes too much time)
        print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
    }
    # # Print a jpeg file
    # file.name <- paste0(name, ".jpg")
    # jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
    # print(list.plots[[1]])
    # dev.off()
    
    # Variables importance ####

    list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)
    
    name <- "VarImp"
    file.name <- paste0(name, ".pdf")
    if( file.exists(paste0(dir.plots.output, info.file.name, file.name)) == F ){ # produce pdf only if it doesn't exist yet (takes too much time)
        print.pdf.plots(list.plots = list.plots, width = 10, height = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
    }
    # # Print a jpeg file
    # file.name <- paste0(name, ".jpg")
    # jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
    # print(list.plots[[1]])
    # dev.off()

    # ICE Manual ####
    
    source("plot_functions.r")
    
    list.list.plots <- lapply(select.taxa, FUN= plot.ice.per.taxa, outputs, list.algo = list.algo[2:4], env.fact, normalization.data, BDM)
    
    for (j in 1:length(select.taxa)) {
        taxon <- sub("Occurrence.", "", j)
        file.name <- paste0("ICE_", taxon, ".pdf")
        print.pdf.plots(list.plots = list.list.plots[[j]], width = 15, height = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
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

