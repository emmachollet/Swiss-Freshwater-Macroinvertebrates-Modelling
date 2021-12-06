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

# # just needed once
# install.packages("reticulate")
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

# Check and set working directory
getwd() # show working directory
# setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts") # set the working directory to this folder

# Free workspace
rm(list=ls())
graphics.off()

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
source("stat_model_functions.r")
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

# preproc.data.env <- preProcess(data.env[,c("SiteId", "SampId", env.fact)], method = c("center", "scale")) 

# Select taxa ####
cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))

all.taxa = F
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa

# Select taxa for prediction
if (all.taxa == T){
    list.taxa <- colnames(data.inv)[cind.taxa]
} else if (all.taxa == F){
    # list.taxa       <- c("Occurrence.Gammaridae", "Occurrence.Heptageniidae") # two taxa
    # list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55), # few taxa
    #                            "Occurrence.taxa"] # Select only few taxa
    list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25), # all taxa with intermediate prevalence
                                "Occurrence.taxa"] # Select with prevalence percentage between 25 and 75%
}

# Construct main dataset (with inv and env)
data <- data.env[, c("SiteId", "SampId", "X", "Y", env.fact)] %>%
    left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data)

# Pre-process data ####
# drop rows with incomplete influence factors: 
# inv.names <- colnames(select(data, contains("Occurrence.group."), contains("Occurrence.")))
# env.names <- colnames(select(data, - all_of(inv.names)))
ind <- !apply(is.na(data[,env.fact]),1,FUN=any)
ind <- ifelse(is.na(ind),FALSE,ind)
data <- data[ind,]
print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))

# - - - - 

# Select ml algorithms ####

# Select models to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c("#030AE8" = 'glm', # Random Forest
               "#048504" = 'bam', # Generalized Additive Model using splines
               # "#948B8B" = 'adaboost',
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # "#A84E05" = 'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
               # "#DB1111" = 'svmRadial', # Support Vector Machine
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
file.name <- paste0(dir.workspace,"SplitsForCV_031221.rds")

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

# Normalize the folds
centered.splits <- lapply(splits, FUN = center.splits, cv = T)

# Normalize the folds but replace '0' ans '1' by factors
centered.splits.factors <- lapply(centered.splits, function(split){
    #split <- centered.splits[[1]]
    return(lapply(split, function(fold){
        
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


## TEST ANN WITH KERAS ####

# Prepare training and testing set

splitted.data <- centered.splits[[1]]

temp.train <- splitted.data$`Training data`[, c(list.taxa,env.fact)]
# temp.train <- na.omit(temp.train)

temp.test <- splitted.data$`Testing data`[, c(list.taxa,env.fact)]
# temp.test <- na.omit(temp.test)

Xtrain <- as.matrix(temp.train[, env.fact])
Ytrain <- as.matrix(temp.train[, list.taxa])

Xtest <- as.matrix(temp.test[, env.fact])
Ytest <- as.matrix(temp.test[, list.taxa])

# One Hot Encoding
# Ytraincat <- to_categorical(as.numeric(Ytrain[,1]) -1)

# Learning rate
lr <- 0.01

# Number of epochs
ne <- 100

# Batch size
bs <-  512

num_layers = 1
num_units_per_layer = 2

# Most models are so-called 'sequential' models
model <- keras_model_sequential()

# Keras makes building neural networks as simple as adding layer upon layer with simple sequential
# calls to the function "layer_dense". Take a moment to appreciate how easy that makes things.

# The input layer is the only layer that requires the user to specify its shape. The shape of all
# subsequent layers is automatically determined based on the output of the preceding layer. Let's
# use a ReLU activation function in each node in the input and hidden layers.
model <- model %>% layer_dense(units = num_units_per_layer,
                               input_shape = ncol(Xtrain),
                               activation = "relu")

# Add the hidden layers. Note this requires just a simple for loop that calls the function "layer_dense"
# again and again.
if (num_layers>1){
    for (i in 1:(num_layers-1)){
        model <- model %>% layer_dense(units = num_units_per_layer,
                                       activation = "relu")
    }
}

# Add the output layer. Note that it uses a sigmoid activation function. Make sure you know why.
model <- model %>% layer_dense(units = 1, activation = "sigmoid")

#  Specify the learning rate for stochastic gradient descent
opt <- optimizer_adam(lr = lr)

# Compile the model, using binary cross-entropy to define loss. Measure accuracy during training.
# Note how easy Keras makes this. Did you have to write any functions for loss or for measuring model
# performance during training? No, Keras takes care of all of this for you.
model %>% compile(optimizer = opt,
                  loss ='binary_crossentropy',
                  metrics = list('accuracy'))

# Fit the model
history <- model %>% fit(x = Xtrain,
                         y = Ytrain,
                         epochs = ne,
                         batch_size = bs)


# 
# 
# build_and_train_model <- function (data = temp.train,
#                                    Xtrain = Xtrain,
#                                    Ytrain = Ytrain,
#                                    learning_rate = lr,
#                                    num_epochs = ne,
#                                    batch_size = bs,
#                                    num_layers = 1,
#                                    num_units_per_layer = 2,
#                                    print_model_summary = F){
# 
#     # Most models are so-called 'sequential' models
#     model <- keras_model_sequential()
# 
#     # Keras makes building neural networks as simple as adding layer upon layer with simple sequential
#     # calls to the function "layer_dense". Take a moment to appreciate how easy that makes things.
# 
#     # The input layer is the only layer that requires the user to specify its shape. The shape of all
#     # subsequent layers is automatically determined based on the output of the preceding layer. Let's
#     # use a ReLU activation function in each node in the input and hidden layers.
#     model <- model %>% layer_dense(units = num_units_per_layer,
#                                    input_shape = ncol(Xtrain),
#                                    activation = "relu")
# 
#     # Add the hidden layers. Note this requires just a simple for loop that calls the function "layer_dense"
#     # again and again.
#     if (num_layers>1){
#         for (i in 1:(num_layers-1)){
#             model <- model %>% layer_dense(units = num_units_per_layer,
#                                            activation = "relu")
#         }
#     }
# 
#     # Add the output layer. Note that it uses a sigmoid activation function. Make sure you know why.
#     model <- model %>% layer_dense(units = 1, activation = "sigmoid")
# 
#     # Print the model description
#     if (print_model_summary){
# 
#         summary(model)
#     }
# 
#     #  Specify the learning rate for stochastic gradient descent
#     opt <- optimizer_adam(lr = learning_rate)
# 
#     # Compile the model, using binary cross-entropy to define loss. Measure accuracy during training.
#     # Note how easy Keras makes this. Did you have to write any functions for loss or for measuring model
#     # performance during training? No, Keras takes care of all of this for you.
#     model %>% compile(optimizer = opt,
#                       loss ='binary_crossentropy',
#                       metrics = list('accuracy'))
# 
#     # Fit the model
#     history <- model %>% fit(x = Xtrain,
#                              y = Ytrain,
#                              epochs = num_epochs,
#                              batch_size = batch_size,
#     )
# 
#     # Return the model and the training history
#     return(list(model = model, history = history))
# }
# 
# ann.outputs <- build_and_train_model