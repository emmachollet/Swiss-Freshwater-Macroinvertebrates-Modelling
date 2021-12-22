## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- TRAINING ANN MODELS ON INVERTEBRATE AND ENVIRONMENTAL DATA ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- December 6, 2021 -- Emma Chollet and Jonas Wydler ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Libraries, data and functions ----

# Load packages
# data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data

# # plots
# if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
# if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") } # to arrange multiple plots on a page
# if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page
# if ( !require("cowplot") ) { install.packages("cowplot"); library("cowplot") } # to arrange multiple plots on a page
# if ( !require("pdp") ) { install.packages("pdp"); library("pdp") } # to plot partial dependance plots
# if ( !require("gt") ) { install.packages("gt"); library("gt") } # to plot nice tables
# if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") } # to plot nice tables
# if(!require(viridis)) {install.packages("viridis", repos="http://cloud.r-project.org"); library(viridis)} # to do even nicer plots
# if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map

# # ml
# if ( !require("caret") ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models
# if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") } # to run GAM ml algorithm
# # if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
# # if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") }
# # if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
# if ( !require("randomForestSRC") ) { install.packages("randomForestSRC"); library("randomForestSRC") } # to run RF and additional features
# if ( !require("keras") ) { devtools::install_github("rstudio/keras"); library("keras") } # to run Neural Networks
# if ( !require("tensorflow") ) { devtools::install_github("rstudio/tensorflow"); library("tensorflow") } # to run Neural Networks
# use_condaenv("r-tensorflow")
#

# ANN

# install.packages("devtools")

# # just needed once
# install.packages("reticulate")
# devtools::install_github("rstudio/reticulate")
library("reticulate")
# install_miniconda()
# 
# install.packages("tensorflow")
library("tensorflow")
# install_tensorflow()
# 
# install.packages("keras")
library("keras")
# install_keras()
# use_condaenv()


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
dir.workspace     <- "../Intermediate results/"
dir.plots.output  <- "../Plots/Models analysis plots/"
dir.models.output <- "../Intermediate results/Trained models/"

file.inv.data     <- "All_occ_data_2020-06-25.dat"
file.inv.BDM.data <- "BDM_occ_data_2020-06-25.dat"
file.env.data     <- "All_environmental_data_2020-06-25.dat"
file.env.BDM.data <- "BDM_environmental_data_2020-06-25.dat"

file.prev         <- "All_prevalence_2020-06-25.dat"
file.prev.BDM     <- "BDM_prevalence_2020-06-25.dat"


# Set if we want to compute for the All or the BDM dataset
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

# set date for file names
d <- Sys.Date()    # e.g. 2021-12-17

# Load functions ####

source("ml_model_functions.r")
source("stat_model_functions.r")
source("plot_functions.r")
source("utilities.r")

# Setup options ####
# set if we want to fit models to whole dataset or perform cross-validation (CV)
CV <- T
dl <- F
#set number of cores
n.cores <-  1

#Settings Stat models 
##Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models, also make sure the cross-validation (CV) flag is
## set correctly
sampsize <- 10 #10000
n.chain  <- 2 #2
comm.corr <- F

#select taxa
all.taxa = T
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa


## ---- DATA WRANGLING  ----

# Select env. factors ####

env.fact <- c("temperature",       # Temp
              "velocity",          # FV
              "A10m",              # A10m
              "cow.density",       # LUD
              "IAR",               # IAR
              "urban.area",        # Urban
              "FRI",               # FRI
              "bFRI",              # bFRI
              "width.variability")#, # WV
#"temperature2",
# "velocity2")
env.fact.full <- c("temperature",       # Temp
                   "velocity",          # FV
                   "A10m",              # A10m
                   "cow.density",       # LUD
                   "IAR",               # IAR
                   "urban.area",        # Urban
                   "FRI",               # FRI
                   "bFRI",              # bFRI
                   "width.variability",#, # WV
                   "temperature2",
                   "velocity2")

#env.fact <- env.fact.full
no.env.fact <- length(env.fact)

# Select taxa ####

cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))

# Select taxa for prediction
if (all.taxa == T){
  list.taxa <- colnames(data.inv)[cind.taxa]
} else if (all.taxa == F){
  
  # 2 taxa
  # list.taxa       <- c("Occurrence.Gammaridae", "Occurrence.Heptageniidae")
  
  # 6 taxa
  # list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55),
  #                            "Occurrence.taxa"] # Select only few taxa
  
  # 22 taxa
  list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25),
                              "Occurrence.taxa"] # Select with prevalence percentage between 25 and 75%
}

no.taxa <- length(list.taxa)

# Summary of prevalence of chosen taxa
for ( i in 1:no.taxa){
  cat("Summary of absence, presence and NA for", list.taxa[i], ":", summary(data.inv[, list.taxa[i]]), "\n")
}

# Select ml algorithms ####

# Select models to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c("#030AE8" = 'glm', # Random Forest
               #"#048504" = 'bam', # Generalized Additive Model using splines
               "#948B8B" = 'gamSpline',#
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # "#A84E05" = 'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
               "#DB1111" = 'svmRadial', # Support Vector Machine
               "#790FBF" = 'rf') # Random Forest

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

# Construct main dataset (with inv and env)
data.full <- data.env[, c("SiteId", "SampId", "X", "Y", env.fact.full)] %>%
  left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data.full)

# drop rows with incomplete influence factors
ind <- !apply(is.na(data.full[,env.fact.full]),1,FUN=any)
ind <- ifelse(is.na(ind),FALSE,ind)
data.full <- data.full[ind,]
print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))
data <- subset(data.full, select = -c(temperature2, velocity2))


#calculate mean and sd for env data for normalisation with data leakage

mean.dl <- apply(select(data.full, all_of(env.fact.full)), 2, function(k){
  mean(k, na.rm = TRUE)
})
sd.dl <- apply(select(data.full, all_of(env.fact.full)), 2, function(k){
  sd(k, na.rm = TRUE)
})


# Split for CV (and save) data ####

# Split for ml purpose (e.g. 80% training, 20% testing)
# ratio <- 0.8 # set a ratio for training dataset (can be 1)
# split.var <- "random" # data splitted according to this variable (default is "random")
# splitted.data <- split.data(data = data, training.ratio = ratio, variable = split.var)

if(CV == T){
  
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
    splits <- split.data(data)
    saveRDS(splits, file = file.name)
  }
}

# Normalize data ####

# Normalize the data (each folds for CV and whole data set else)
if(CV == T){
  
  #center the splits
  centered.splits.tmp <- lapply(splits, FUN = center.data, CV = CV, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl)
  #centered.splits.tmp <- lapply(splits, FUN = center.data, CV = CV)
  #extract necessary information
  centered.splits <- lapply(centered.splits.tmp,"[", 1:2) # only the splits without the mean, sd info
  normalization.data <- lapply(centered.splits.tmp,"[", 3:4) # the mean and sd of the splits
  
  # Normalize the folds but replace '0' and '1' by factors
  centered.splits.factors <- lapply(centered.splits, function(split){
    #split <- centered.splits[[1]] #to test
    return(lapply(split, function(fold){
      #fold <- split[[1]] # to test
      cind.taxa <- which(grepl("Occurrence.",colnames(fold)))
      #Replace "0" and "1" by "absent" and "present" and convert them to factors
      for (i in cind.taxa ) {
        fold[which(fold[,i] == 0),i] <- "absent"
        fold[which(fold[,i] == 1),i] <- "present"
        fold[,i] = as.factor(fold[,i])
      }
      return(fold)
    })
    )
  })
} else {
  
  centered.data <- center.data(list(data), CV = CV, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl)
  centered.data.factors <- centered.data
  normalization.data <- list(mean.dl, sd.dl) #without cv we just use the global mean and sd
  # Replace '0' and '1' by factors
  cind.taxa <- which(grepl("Occurrence.",colnames(centered.data.factors[[1]])))
  #Replace "0" and "1" by "absent" and "present" and convert them to factors
  for (i in cind.taxa ) {
    centered.data.factors[[1]][which(centered.data.factors[[1]][,i] == 0),i] <- "absent"
    centered.data.factors[[1]][which(centered.data.factors[[1]][,i] == 1),i] <- "present"
    centered.data.factors[[1]][,i] = as.factor(centered.data.factors[[1]][,i])
  }
}
### Test centering
if(exists("centered.splits") == T){
  if(mean(centered.splits$Split1$`Training data`$temperature) <= 0.001){
    cat("The data is normalized.")
  }else{
    cat("The data isn't normalized.")
    break()
  }
}else if (exists("centered.data") == T){
  if(mean(centered.data$`Entire dataset`$temperature) <= 0.001){
    cat("The data is normalized.")
  }else{
    cat("The data isn't normalized.")
    break()
  }
}else{
  cat("The data isn't normalized.")
  
}

## TEST ANN WITH KERAS ####

# Prepare train and test sets ####

splitted.data <- as.data.frame(centered.data[[1]])

splitted.data <- as.data.frame(centered.splits[[1]][[1]])

# spot col with too many NAs
too.many.na <- c()
for(i in 1:dim(splitted.data)[2]){
  if(sum(is.na(splitted.data[,i])) > 200){ too.many.na <- c(too.many.na, i)}
}

# remove col with too many NA
splitted.data <- splitted.data[, -too.many.na]

# remove remaining NA
splitted.data.NA <- splitted.data
splitted.data <- na.omit(splitted.data)


list.taxa <- which(grepl("Occurrence.", colnames(splitted.data)))

# # Write a dataframe on dat for Andreas
# filename <- paste(dir.models.output, "TrainingData1", ".dat", sep="")
# write.table(splitted.data$`Training data`, filename,sep="\t",row.names=F,col.names=TRUE)
set.seed(2021)
folds <- groupKFold(splitted.data$SiteId, 3)

train <- splitted.data[folds$Fold1,]
test <- splitted.data[-folds$Fold1,]

# data <- saved.data

Xtrain <- as.matrix(train[ ,env.fact])
Ytrain <- as.matrix(train[ ,list.taxa])

Xtest <- as.matrix(test[ ,env.fact])
Ytest <- as.matrix(test[ ,list.taxa])
# no.fact <- ncol(X)
# no.taxa <- ncol(Y)
# 
# folds.NA <- groupKFold(splitted.data.NA$SiteId, 3)
# 
# train.NA <- splitted.data.NA[folds.NA$Fold1,]
# test.NA <- splitted.data.NA[-folds.NA$Fold1,]
# 
# # data <- saved.data
# 
# Xtrain.NA <- as.matrix(train.NA[ ,env.fact])
# Ytrain.NA <- as.matrix(train.NA[ ,list.taxa])
# 
# Xtest.NA <- as.matrix(test.NA[ ,env.fact])
# Ytest.NA <- as.matrix(test.NA[ ,list.taxa])
# 

# Run the ANN #### 

# Learning rate
learning.rate <- 0.01

# Number of epochs
num.epochs <- 20

# Batch size
batch.size <-  64


# hyper.param <- c(1,15)

list.hyper.param <- list(c(3, 32), c(3, 64), c(5, 32))
names(list.hyper.param) <- lapply(
  lapply(list.hyper.param, FUN = paste, c("L", "U"), sep = ""), # paste number of Layers and Units with "L" and "U"
  FUN = paste, collapse = "") # paste the two parts together
names(list.hyper.param) <- paste("ANN", names(list.hyper.param), sep = "")

build_and_train_model <- function (hyper.param = hyper.param,
                                   Xtrain = Xtrain,
                                   Ytrain = Ytrain,
                                   Xtest = Xtest,
                                   Ytest = Ytest,
                                   learning.rate = learning.rate,
                                   num.epochs = num.epochs,
                                   batch.size = batch.size,
                                   CV = T){

    num.layers <- hyper.param[1]
    num.units.per.layer <- hyper.param[2]
    
    # Most models are so-called 'sequential' models
    model <- keras_model_sequential()
    
    model <- model %>% layer_dense(units = num.units.per.layer,
                                   input_shape = ncol(Xtrain),
                                   activation = "tanh")
    
    if (num.layers>1){
      for (i in 1:(num.layers-1)){
        model <- model %>% layer_dense(units = num.units.per.layer,
                                       activation = "tanh")
      }
    }
    
    # Add the output layer. Note that it uses a sigmoid activation function. Make sure you know why.
    model <- model %>% layer_dense(units = ncol(Ytrain), activation = "sigmoid")
    
    
    
    summary(model)
    
    
    #  Specify the learning rate for stochastic gradient descent
    opt <- optimizer_adam(learning_rate = learning.rate)
    
    
    model %>% compile(optimizer = opt,
                      loss = "binary_crossentropy", #loss_binary_crossentropy(),
                      metrics = list('accuracy')) # categorical_accuracy ?
    
    # Fit the model
    history <- model %>% fit(x = Xtrain,
                             y = Ytrain,
                             epochs = num.epochs,
                             batch_size = batch.size)
    
    
    list.taxa <- colnames(Ytrain)
    
    # Make a list with the outputs of the algorithm for each taxon in list.taxa
    list.outputs <- vector(mode = 'list', length = length(list.taxa))
    names(list.outputs) <- list.taxa
    
    if(CV == T){which.set <- c("training set", "testing set")
    } else {which.set <- c("training set")}
    out <- c("Observation", #1
             "Prediction factors", #2 
             "Prediction probabilities", #3 
             "Likelihood", #4
             "Performance") #5
    output.names <- c("Trained model", "Training history", c(outer(out, which.set, FUN = paste)))
    
    pred <- list()
    likeli <- list()
    
    for (n in 1:length(which.set)) {
      # Predict 
      pred[[n]] <- model %>% predict(if(n == 1){ Xtrain } else { Xtest })
      
      # Compute standardized deviance
      likeli[[n]] <- if(n == 1){ Ytrain } else { Ytest }
      for (j in 1:ncol(likeli[[n]])) {
        for (i in 1:nrow(likeli[[n]])) {
          obs <- likeli[[n]][i,j]
          likeli[[n]][i,j] <- ifelse(obs == 1, pred[[n]][i,j], 1 - pred[[n]][i,j])
        }
      }
    }
    
    for (j in 1:length(list.taxa)) {
      temp.list <- vector(mode = 'list', length = length(output.names))
      names(temp.list) <- output.names
      
      temp.sets <- if(CV == T){ list(cbind(Xtrain, Ytrain[,j]), cbind(Xtest, Ytest[,j])) } else { list(cbind(Xtrain, Ytrain[,j])) }
      
      temp.list[[1]] <- model
      temp.list[[2]] <- history
      
      for(n in 1:length(which.set)){
      pred.prob <- pred[[n]][,j]
      pred.fact <- ifelse(pred.prob > 0.5, "present", "absent")
      temp.list[[paste(out[1],which.set[n])]] <- temp.sets[[n]]
      temp.list[[paste(out[2],which.set[n])]] <-  pred.fact
      temp.list[[paste(out[3],which.set[n])]] <-  pred.prob
      temp.list[[paste(out[4],which.set[n])]] <-  likeli[[n]][,j]
      perf <- -2 * sum(log(likeli[[n]][,j])) / nrow(temp.sets[[n]])
      temp.list[[paste(out[5],which.set[n])]] <- perf
      }
      list.outputs[[j]] <- temp.list
    }
    
    # Return the model and the training history
    return(list.outputs)
}


outputs.hyperparam <- lapply(list.hyper.param, FUN = build_and_train_model, Xtrain = Xtrain,
                             Ytrain = Ytrain,
                             Xtest = Xtest,
                             Ytest = Ytest,
                             learning.rate = learning.rate,
                             num.epochs = num.epochs,
                             batch.size = batch.size,
                             CV = T)

file.name <- paste0(dir.models.output, "test_ANNoutputs.rds")
saveRDS(outputs.hyperparam, file = file.name)



# WORKS UNTIL HERE ####



output.wo.NA <- build_and_train_model(hyper.param = hyper.param, Xtrain = Xtrain, Ytrain = Ytrain, 
                                      Xtest = Xtest, learning.rate = learning.rate, num.epochs = num.epochs,
                                      batch.size = batch.size)
output.wi.NA <- build_and_train_model(hyper.param = hyper.param, Xtrain = Xtrain.NA, Ytrain = Ytrain.NA, 
                                      Xtest = Xtest.NA, learning.rate = learning.rate, num.epochs = num.epochs,
                                      batch.size = batch.size)




image(output.wo.NA$prediction[1:200,])




dim(pred)

image(pred[1:200,])

# mod.loss <- loss_binary_crossentropy(Y,pred)
# k_mean(mod.loss)
# 
# Xtest <- matrix(c(0,1), ncol = 11, nrow = 2)
# pred.test <- model %>% predict(Xtest)
# 
# get_layer(pred)$output
# 
# expl.var <- env.fact
# target.var <- list.taxa
# 
# temp.train <- splitted.data$`Training data`[, c(target.var,expl.var)]
# # temp.train <- na.omit(temp.train)
# 
# temp.test <- splitted.data$`Testing data`[, c(target.var,expl.var)]
# # temp.test <- na.omit(temp.test)
# 
# Xtrain <- as.matrix(temp.train[7:9, expl.var])
# Ytrain <- as.matrix(temp.train[7:9, target.var])
# 
# 
# 
# Xtest <- as.matrix(temp.test[7:9, expl.var])
# Ytest <- as.matrix(temp.test[7:9, target.var])

# One Hot Encoding
# Ytraincat <- to_categorical(as.numeric(Ytrain[,1]) -1)

# Learning rate
learning.rate <- 0.1

# Number of epochs
num.epochs <- 20

# Batch size
batch.size <-  1

# list.hyper.param <- list("hyper1" = c(3, 32), "hyper2" = c(3, 64), "hyper3" = c(5, 32))

# build_and_train_model <- function (hyper.param = hyper.param,
#                                    # no.env.fact = no.env.fact,
#                                    Xtrain = Xtrain,
#                                    Ytrain = Ytrain,
#                                    Xtest = Xtest,
#                                    learning.rate = learning.rate,
#                                    num.epochs = num.epochs,
#                                    batch.size = batch.size,
#                                    print.model.summary = T){

  # hyper.param <- list.hyper.param[[1]]
hyper.param <- c(1,2)
    num.layers <- hyper.param[1]
    num.units.per.layer <- hyper.param[2]
    
    # Most models are so-called 'sequential' models
    model <- keras_model_sequential()

    # Keras makes building neural networks as simple as adding layer upon layer with simple sequential
    # calls to the function "layer_dense". Take a moment to appreciate how easy that makes things.

    # The input layer is the only layer that requires the user to specify its shape. The shape of all
    # subsequent layers is automatically determined based on the output of the preceding layer. Let's
    # use a ReLU activation function in each node in the input and hidden layers.
    model <- model %>% layer_dense(units = num.units.per.layer,
                                   input_shape = ncol(Xtrain), # PROBLEM HERE ####
                                   activation = "tanh")

    # Add the hidden layers. Note this requires just a simple for loop that calls the function "layer_dense"
    # again and again.
    if (num.layers>1){
        for (i in 1:(num.layers-1)){
            model <- model %>% layer_dense(units = num.units.per.layer,
                                           activation = "tanh")
        }
    }

    # Add the output layer. Note that it uses a sigmoid activation function. Make sure you know why.
    model <- model %>% layer_dense(units = ncol(Ytrain), activation = "sigmoid")

    # Print the model description
    # if (print.model.summary){

        summary(model)
    # }

    #  Specify the learning rate for stochastic gradient descent
    opt <- optimizer_adam(learning_rate = learning.rate)


    # WRITE THIS my.loss fct in R with no loop and keras functions (e.g. k_log)
    
    # Compile the model, using binary cross-entropy to define loss. Measure accuracy during training.
    # Note how easy Keras makes this. Did you have to write any functions for loss or for measuring model
    # performance during training? No, Keras takes care of all of this for you.
    model %>% compile(optimizer = opt,
                      loss = loss_binary_crossentropy(),
                      metrics = list('accuracy')) # categorical_accuracy ?

    # Fit the model
    history <- model %>% fit(x = Xtrain,
                             y = Ytrain,
                             epochs = num.epochs,
                             batch_size = batch.size)

    # Predict on testing set
    pred <- model %>% predict(Xtrain)
    
    
    colnames(pred) <- target.var
    
    # Return the model and the training history
#    return(list(model = model, history = history, prediction = pred))
# }




my.loss(Ytest, pred)
loss_binary_crossentropy(Ytest,pred)

test <- build_and_train_model(hyper.param = list.hyper.param[[1]], print.model.summary = T)

ann.outputs <- lapply(list.hyper.param, FUN = build_and_train_model, no.env.fact, Xtrain, Ytrain, Xtest, lr, ne, bs)

file.name <- paste0(dir.models.output, "output_ANN.rds")
saveRDS(ann.outputs, file = file.name, version = 2)


# test loss functons

# my.loss <- function(obs, pred){
#   loss <- 0
#   loss.persite <- NULL
#   for (i in 1:nrow(pred)) {
#     loss.persite.temp <- 0
#     for (j in ncol(pred)) {
#       loss <- loss + ifelse(obs[i,j] == 1, log(pred[i,j]), log(1 - pred[i,j]))
#       loss.persite.temp <- loss.persite.temp + ifelse(obs[i,j] == 1, log(pred[i,j]), log(1 - pred[i,j]))
#     }
#     loss.persite <- c(loss.persite, -1/(ncol(pred))*loss.persite.temp)
#   }
#   loss <- -1/(ncol(pred)*nrow(pred))*loss
#   return(list(loss, loss.persite))
#   # return(loss)
# }

loss.site <- function(obs, pred){
  -mean(ifelse(obs == 1, log(pred), log(1 - pred)))
}

my.loss <- function(obs, pred){
  loss <- rep(NA, nrow(obs))
  for (site in 1:nrow(obs)) {
    loss[site] <- loss.site(obs[site,], pred[site,])
  }
  return(loss)
}

y <- matrix(as.integer(c(0,1,0,0,0,0,NA,NA)), ncol = 4)
p <- matrix(c(0.1,0.1,0.1,0.1,0.5,0.7,0.3,0.5), ncol = 4)

loss_binary_crossentropy(y,p)
my.loss(y,p)

k_mean(loss_binary_crossentropy(Y,pred))

k_mean(loss_binary_crossentropy(y,pred[1:2,4:7]))
k_mean(loss_binary_crossentropy(ytest,pred[1:2,4:7]))

k_mean(loss_binary_crossentropy(Y[1:2,4:7],pred[1:2,4:7]))

ytest <- Y[1:2,4:7]
ytest == y

loss.site(Y[1,], pred[1,])



# # Most models are so-called 'sequential' models
# model <- keras_model_sequential()
# 
# # Keras makes building neural networks as simple as adding layer upon layer with simple sequential
# # calls to the function "layer_dense". Take a moment to appreciate how easy that makes things.
# 
# # The input layer is the only layer that requires the user to specify its shape. The shape of all
# # subsequent layers is automatically determined based on the output of the preceding layer. Let's
# # use a ReLU activation function in each node in the input and hidden layers.
# model <- model %>% layer_dense(units = num.units.per.layer,
#                                input_shape = ncol(Xtrain),
#                                activation = "relu")
# 
# # Add the hidden layers. Note this requires just a simple for loop that calls the function "layer_dense"
# # again and again.
# if (num.layers>1){
#     for (i in 1:(num.layers-1)){
#         model <- model %>% layer_dense(units = num.units.per.layer,
#                                        activation = "relu")
#     }
# }
# 
# # Add the output layer. Note that it uses a sigmoid activation function. Make sure you know why.
# model <- model %>% layer_dense(units = ncol(Ytrain), activation = "sigmoid")
# 
# #  Specify the learning rate for stochastic gradient descent
# opt <- optimizer_adam(lr = lr)
# 
# # Compile the model, using binary cross-entropy to define loss. Measure accuracy during training.
# # Note how easy Keras makes this. Did you have to write any functions for loss or for measuring model
# # performance during training? No, Keras takes care of all of this for you.
# model %>% compile(optimizer = opt,
#                   loss ='binary_crossentropy',
#                   metrics = list('accuracy'))
# 
# # Fit the model
# history <- model %>% fit(x = Xtrain,
#                          y = Ytrain,
#                          epochs = ne,
#                          batch.size = bs)


