## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- TRAINING ML MODELS ON INVERTEBRATE AND ENVIRONMENTAL DATA ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- December 1, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- PRELIMINARIES ----

# Libraries ####

# Load packages
# data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data

# plots
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") } # to arrange multiple plots on a page
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page
if ( !require("cowplot") ) { install.packages("cowplot"); library("cowplot") } # to arrange multiple plots on a page
if ( !require("pdp") ) { install.packages("pdp"); library("pdp") } # to plot partial dependance plots
if ( !require("gt") ) { install.packages("gt"); library("gt") } # to plot nice tables
if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") } # to plot nice tables
if(!require(viridis)) {install.packages("viridis", repos="http://cloud.r-project.org"); library(viridis)} # to do even nicer plots
if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map

# ml
if ( !require("caret") ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models
if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") } # to run generalized additive model (GAM) algorithm
# if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") } # to run support vector machine (svm) algorithm
# if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
if ( !require("randomForest") ) { install.packages("randomForest"); library("randomForest") } # to run random forest (RF)
# if ( !require("randomForestSRC") ) { install.packages("randomForestSRC"); library("randomForestSRC") } # to run RF and additional features
# if ( !require("keras") ) { devtools::install_github("rstudio/keras"); library("keras") } # to run Neural Networks
# if ( !require("tensorflow") ) { devtools::install_github("rstudio/tensorflow"); library("tensorflow") } # to run Neural Networks
# use_condaenv("r-tensorflow")

# # needed for ANN
if( !require("reticulate")){
  # remotes::install_github("rstudio/reticulate"); # this tried to install reticulate in a inaccessible folder
  install.packages("reticulate");
  library("reticulate")}
use_condaenv()
# 
# if( !require("tensorflow")){
#   # devtools::install_github("rstudio/tensorflow");
#   install.packages("tensorflow");
#   library(tensorflow)}
# # tf_version()
# 
# install_tensorflow()
# 
if( !require("keras")){
  devtools::install_github("rstudio/keras");
  # install.packages("keras");
  library(keras)}
# 
# # library("tensorflow")
# # library("keras")
# 
# install_keras()

## Run this to prompt miniconda installation request!
## If no request is prompted, you're ready to start coding! Else, press "n".
# set_random_seed(0)

# Check and set working directory
getwd() # show working directory
# setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts") # set the working directory to this folder

# Free workspace
rm(list=ls())
graphics.off()

# Load data ####

# Define directory and files
dir.env.data      <- "../../Data/Processed data/Environmental data/"
dir.inv.data      <- "../../Data/Processed data/Invertebrate data/"
dir.output        <- "../Plots/Models analysis plots/"
dir.workspace     <- "../Intermediate results/"

file.inv.data     <- "All_occ_data_2020-06-25.dat"
file.inv.BDM.data <- "BDM_occ_data_2020-06-25.dat"
file.env.data     <- "All_environmental_data_2020-06-25.dat"
file.env.BDM.data <- "BDM_environmental_data_2020-06-25.dat"

file.prev         <- "All_prevalence_2020-06-25.dat"
file.prev.BDM     <- "BDM_prevalence_2020-06-25.dat"

file.rank.env     <- "ranking_env_data.csv"

# Read data
# needed to sort environmental factors by their importance
rank.env          <- read.csv(paste(dir.env.data,file.rank.env,sep=""),header=TRUE, sep=";", stringsAsFactors=FALSE)

# set if we want to compute for the All or the BDM dataset
BDM <- F

if( BDM == TRUE){

    data.inv      <- read.delim(paste(dir.inv.data,file.inv.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=F)
    prev.inv      <- read.delim(paste(dir.inv.data,file.prev.BDM,sep=""),header=T,sep="\t", stringsAsFactors=F)
    data.env      <- read.delim(paste(dir.env.data,file.env.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=T)
    file.prefix   <- "BDM_"
    
} else {
    
    data.inv          <- read.delim(paste(dir.inv.data,file.inv.data,sep=""),header=T,sep="\t", stringsAsFactors=F)
    prev.inv          <- read.delim(paste(dir.inv.data,file.prev,sep=""),header=T,sep="\t", stringsAsFactors=F)
    data.env          <- read.delim(paste(dir.env.data,file.env.data,sep=""),header=T,sep="\t", stringsAsFactors=T)
    file.prefix       <- "All_"
    
}

d <- paste0(Sys.Date(), "_")    # date for file names

# Load functions ####

source("ml_model_functions.r")
source("plot_functions.r")
source("utilities.r")

## ---- DATA WRANGLING  ----

# Select env. factors ####

# Replace "_" and " " by "." in colnames to be consistent
colnames(rank.env) <- gsub("_", ".", colnames(rank.env))
colnames(rank.env) <- gsub(" ", ".", colnames(rank.env))

# Sort the columns in 4 categories: the sample/site information = 3,
# the environmental factors to keep in priority = 2, keep to explore = 1, exclude = 0
info    <- colnames(rank.env)[which(rank.env[1,] == 3)]
prio    <- colnames(rank.env)[which(rank.env[1,] == 2)]
explo   <- colnames(rank.env)[which(rank.env[1,] == 1)]
excl    <- colnames(rank.env)[which(rank.env[1,] == 0)]

# decide to explore the selected (Bogdan's) factors or the list "priority"
select = T    

if (select == F){
  
  env.fact <- prio
  
} else if (select == T){
  
  env.fact <- c("temperature",       # Temp
                "velocity",          # FV
                "A10m",              # A10m
                "cow.density",       # LUD
                "IAR",               # IAR
                "urban.area",        # Urban
                "FRI",               # FRI
                "bFRI",              # bFRI
                "width.variability", # WV
                "temperature2",
                "velocity2")
}

if(BDM != TRUE) {
  # remove InS env. fact.
  cind <- c(grep("InS.",env.fact),
            grep("covclass.",env.fact),
            grep("normcov.",env.fact),
            grep("sfract.",env.fact),
            grep("v.",env.fact, fixed = TRUE)) # has to match exactly, to avoid to remove .variability or .vegetation
  if(length(cind) != 0){
    env.fact <- env.fact[-cind]
  }
}

# Remove env. fact. missing in data.env
env.fact <- env.fact[which(env.fact %in% colnames(data.env))]
# env.fact[-which(env.fact %in% colnames(data.env))] # which selected factors are not in data.env ?
# [1] "A.EDO" "F.EDO" are missing

# Standardize env. factors ####

preproc.data.env <- preProcess(data.env[,c("SiteId", "SampId", env.fact)], method = c("center", "scale")) 

# Select taxa ####

cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
 
# Replace "0" and "1" by "absent" and "present" and convert them to factors
for (i in cind.taxa ) {
    data.inv[which(data.inv[,i] == 0),i] <- "absent"
    data.inv[which(data.inv[,i] == 1),i] <- "present"
    data.inv[,i] = as.factor(data.inv[,i])
}

# Construct main dataset (with inv and env)
data <- data.env[, c("SiteId", "SampId", env.fact)] %>%
  left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data)


# Select taxa
all.taxa = F 
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa

# Select taxa for prediction
if (all.taxa == T){
  list.taxa <- colnames(data.inv)[cind.taxa]
} else if (all.taxa == F){
  # list.taxa       <- c("Occurrence.Gammaridae", "Occurrence.Heptageniidae") # two taxa
  list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55), # few taxa
                              "Occurrence.taxa"] # Select only few taxa
  # list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25), # all taxa with intermediate prevalence
  #                             "Occurrence.taxa"] # Select with prevalence percentage between 25 and 75%
}

# - - - - 

# Select ml algorithms ####

# Select models to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c("#030AE8" = 'glm', # Random Forest
               # "#048504" = 'bam', # Generalized Additive Model using splines
               # "#948B8B" = 'adaboost',
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # "#A84E05" = 'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
               "#DB1111" = 'svmRadial', # Support Vector Machine
               "#790FBF" = 'rf') # Random Forest

# Assemble information to insert in file names
no.env.fact <- length(env.fact)
no.taxa <- length(list.taxa)
no.algo <- length(list.algo)

# Write information for file names
# percentage.train.set <- ratio * 100
info.file.name <- paste0(file.prefix, 
                         # d, # don't need to include the date
                         no.taxa, "taxa_", 
                         no.env.fact, "envfact_",
                         no.algo, "algo_",
                         # "trainset", percentage.train.set, 
                         # if( ratio != 1) {split.var}, 
                         "_")

# Summary of prevalence of chosen taxa
for ( i in 1:no.taxa){
    cat("Summary of absence, presence and NA for", list.taxa[i], ":", summary(data[, list.taxa[i]]), "\n")
}


# Split (and save) data ####

# Split for ml purpose (e.g. 80% training, 20% testing)
# ratio <- 0.8 # set a ratio for training dataset (can be 1)
# split.var <- "random" # data splitted according to this variable (default is "random")
# splitted.data <- split.data(data = data, training.ratio = ratio, variable = split.var)

# Split for CV
file.name <- paste0(dir.workspace,"SplitsForCV.rds")

# If the file with the three different splits already exist, just read it
if (file.exists(file.name) == T ){
  
  if(exists("splits") == F){ splits <- readRDS(file = file.name)
    cat("File with data splits already exists, we read it from", file.name, "and save it in object 'splits'")}
    else{
    cat("List with data splits already exists as object 'splits' in this environment.")
    }
  } else {
  
  cat("No data splits exist yet, we produce it and save it in", file.name)
  splits <- split.data(data, 1)
  saveRDS(splits, file = file.name)

}

# Normalize data ####

centered.splits <- lapply(splits, FUN = center.splits, cv = T)

#pp <- preProcess(splits[[]])


## ---- Apply stat model ----


# read in results produced by Jonas
file.name <- paste0(dir.workspace, "Output25112021.rds")

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T ){
    
    if(exists("outputs") == F){
        cat("File with statistical model outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
        stat.outputs <- readRDS(file = file.name)}
    else{
        cat("List with statistical model outputs already exists as object 'outputs' in this environment.")
    }
} else {
    
    cat("No statistical model outputs exist yet, we produce it and save it in", file.name)
    
    # #1) No cross validation ####
    # centered.occ.data <- center.splits(list(inv.occ), cv = F)
    # #a) No cross validation, no comm corr
    # res.1 <- stat_mod_cv(data.splits = centered.occ.data, cv = F, comm.corr = F)
    # 
    # #b) No cross validation, comm corr
    # res.1 <- stat_mod_cv(data.splits = centered.occ.data, cv = F, comm.corr = T)
    # 
    # #2) Cross validation ####
    # centered.splits <- lapply(splits, FUN = center.splits, cv = T)
    # 
    # #b) Cross validation, no comm corr
    # res.3 <- mclapply(centered.splits, mc.cores = 3, FUN = stat_mod_cv, cv = T, comm.corr = F)
    # 
    # #b) Cross validation, comm corr
    # res.4 <- mclapply(centered.splits, mc.cores = 3, FUN = stat_mod_cv, cv = T, comm.corr = T)
    
    # cat("Saving outputs of statistical model in", file.name)
    # saveRDS(stat.outputs, file = file.name)
    
}

## ---- Apply ML models ----

ptm <- proc.time() # to calculate time of simulation

## TEST ANN WITH KERAS ####

# # Prepare training and testing set
# target <- list.taxa[1]
# 
# temp.train <- splits[[1]][[1]][, c(target,env.fact)]
# temp.train <- na.omit(temp.train)
# 
# temp.test <- splits[[1]][[2]][, c(target,env.fact)]
# temp.test <- na.omit(temp.test)
# 
# # Xtrain <- as.matrix(temp.train[, env.fact])
# # Ytrain <- as.matrix(temp.train[, target])
# # 
# # Xtest <- as.matrix(temp.test[, env.fact])
# # Ytest <- as.matrix(temp.test[, target])
# 
# Xtrain <- temp.train[, env.fact]
# Ytrain <- temp.train[, target]
# 
# Xtest <- temp.test[, env.fact]
# Ytest <- temp.test[, target]
# 
# # One Hot Encoding
# # Ytraincat <- to_categorical(as.numeric(Ytrain[,1]) -1)
# 
# # use_session_with_seed(42)
# 
# # # Initialize a sequential model
# model <- keras_model_sequential()
# 
# model %>%
#   layer_dense(units = 12, activation = 'relu', input_shape = c(20)) %>%
#   layer_dense(units = 12, activation = 'relu') %>%
#   layer_dense(units = 10, activation = 'sigmoid')
# 
# model %>% compile(
#   loss = 'binary_crossentropy',
#   optimizer = optimizer_rmsprop()
# )
# 
# model %>% fit(
#   x = Xtrain, y = Ytrain, epochs = 20, batch_size = 32
# )

# intermediately, we need to store one of the splits in splitted data to make 
# the ML algorithms running
splitted.data <- splits[["Split1"]]


# "Apply" null model
null.model <- apply.null.model(data = data, list.taxa = list.taxa, prev.inv = prev.inv)

file.name <- paste0(dir.workspace, no.algo, "MLAlgoTrained.rds")

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T ){
    
    if(exists("outputs") == F){
        cat("File with ML outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
        outputs <- readRDS(file = file.name)}
    else{
        cat("List with ML outputs already exists as object 'outputs' in this environment.")
    }
} else {
    
    cat("No ML outputs exist yet, we produce it and save it in", file.name)
    # Make a list to store the outputs of each model
    outputs <- vector(mode = 'list', length = no.algo)
    names(outputs) <- list.algo
    
    # Apply the ML models to list.taxa and env.fact
    for(k in 1:no.algo){
        outputs[[k]] <- apply.ml.model(splitted.data = splitted.data, list.taxa = list.taxa,
                                       env.fact = env.fact, algorithm = list.algo[k])
    }
    cat("Saving outputs of algorithms in", file.name)
    saveRDS(outputs, file = file.name)
    
}

print(paste("Simulation time of different models ", info.file.name))
print(proc.time()-ptm)

# Training GAM bam for 6 taxa: 5 hours

# source("ml_model_functions.r")
# source("plot_functions.r")
# rm(list=ls())
# graphics.off()

## ---- Plot models comparison ----

ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- model.comparison(outputs = outputs, null.model = null.model, list.algo = list.algo, list.taxa = list.taxa, prev.inv = prev.inv)

# Print the plots in a pdf file
file.name <- "ModelsCompar.pdf"
print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot trained model performance against hyperparameters ----

ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- plot.perf.hyperparam(outputs = outputs, list.algo = list.algo[2:3], list.taxa = list.taxa)

# Print the plots in a pdf file
file.name <- "PerfvsHyperparam.pdf"
print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)


## ---- Plot variable importance ----

ptm <- proc.time() # to calculate time of simulation

list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)

file.name <- "VarImp.pdf"
print.pdf.plots(list.plots = list.plots, width = 10, height = 10, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

# print directly table with variable importance for each algo and taxa (too complicated to put in a fct)
file.name <- "TableVarImp.pdf"
pdf(paste0(dir.output, info.file.name, file.name), paper = 'special', width = 12, height = 9, onefile = TRUE)
temp.df <- data.frame(matrix(ncol = no.algo*no.taxa, nrow = no.env.fact))
colnames(temp.df) <- c(outer(list.algo, list.taxa, FUN = paste))
rownames(temp.df) <- env.fact
for (j in 1:no.taxa) {
  for (l in 1:no.algo) {
    for (k in 1:no.env.fact) {
      temp.df[env.fact[k],paste(list.algo[l], list.taxa[j])] <- outputs[[l]][[j]][["Variable importance"]][["importance"]][env.fact[k],1]
    }
  }
}
temp.df$mean.imp <- rowMeans(temp.df)
temp.df <- as.matrix(temp.df)
par(mar=c(1,5,15,3)+ 0.2, xaxt = "n")
plot(temp.df, 
     #key = NULL,
     digits = 2, text.cell=list(cex=0.5),
     # axis.col=list(side=3, las=2), 
     axis.row = list(side=2, las=1),
     col = viridis,
     xlab = "",
     ylab = "",
     cex.axis = 0.5,
     srt = 45,
     main = "Variable importance for ML algorithm applied to taxa"
     )
axis(1, at=seq(1:ncol(temp.df)+1), labels = FALSE)
text(seq(1:ncol(temp.df)+1), par("usr")[4] + 0.15, srt = 50, 
     labels = colnames(temp.df), adj= 0, cex = 0.5, xpd = T)

dev.off()

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# For BDM Dataset, only 55 env fact (priority) : 3 sec


## ---- Print accuracy ----

for(l in 1:no.algo){
  for(j in 1:no.taxa){
      # cat("Accuracy on training for",list.taxa[j], "is :",
      #    output[[j]][["Trained model"]][["resample"]][["Accuracy"]], "\n")
      cat("Accuracy on prediction for",list.taxa[j],"with", list.algo[l], "is :",
          outputs[[l]][[j]][["Confusion matrix testing set"]][["overall"]][["Accuracy"]], "\n")
  }
}


# DON'T WORK FROM HERE ####
## ---- Plot PDP ----

ptm <- proc.time() # to calculate time of simulation

# PDP of one model
# list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
#                       list.taxa = list.taxa, env.fact = env.fact)
# 
# file.name <- "PDP.pdf"
# print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

# PDP of all models
# We sub-select taxa and env.fact because it takes a lot of time
list.plots <- plot.pdp(outputs = outputs, list.algo = list.algo,
                       list.taxa = list.taxa[1:2], env.fact = env.fact)

file.name <- "allPDP.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot single predictor ICE ----

ptm <- proc.time() # to calculate time of simulation

# ICE of one model
list.plots <- plot.ice(outputs = outputs, algo = 'rf', list.algo = list.algo,
                       list.taxa = list.taxa[1:2], env.fact = env.fact[1:2])

file.name <- "ICE.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot PDP (multiple predictors) ----

# We just make one example because it's computationally heavy

ptm <- proc.time() # to calculate time of simulation

# Multiple (2) predictors PDP (of one model) for now just for 1 algo and 1 taxa
list.plots <- plot.mult.pred.pdp(outputs = outputs, list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "multpredPDP.pdf"
print.pdf.plots(list.plots = list.plots, width = 17, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

print("hey")
# Stop execution if there is no (prediction on) testing set
# stopifnot(ratio != 1)
if( length(outputs[[1]][[1]]) < 7){
  break
}
print("salut")

## ---- Plot map prediction ----

ptm <- proc.time() # to calculate time of pdf production

map.inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

# make a list with all plots and plot them in a pdf
list.plots <- map.ml.pred.taxa(inputs = map.inputs, outputs = outputs,
                               data.env = data.env, # for now it needs data.env to bind columns X and Y, but maybe delete later
                               list.taxa = list.taxa, list.algo = list.algo,)

file.name <- "ObsvsPred_map.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot env fact vs taxa prediction ----

ptm <- proc.time() # to calculate time of pdf production

# source("plot_functions.r")
list.plots <- response.ml.pred.taxa(outputs = outputs, list.algo = list.algo,
                              list.taxa = list.taxa, env.fact = env.fact)

file.name <- "Resp_EnvFactvsTax.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print("Producing PDF time:")
print(proc.time()-ptm)
