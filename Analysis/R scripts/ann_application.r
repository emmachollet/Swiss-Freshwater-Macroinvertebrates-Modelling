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

file.env.data     <- "environmental_data_2020-06-25.dat"
file.inv.data     <- "occ_data_2020-06-25.dat"
file.prev         <- "prevalence_2020-06-25.dat"

# Set if we want to compute for the All or the BDM dataset
BDM <- F

file.prefix <- ifelse(BDM, "BDM_", "All_")

data.env          <- read.delim(paste0(dir.env.data, file.prefix, file.env.data),header=T,sep="\t", stringsAsFactors=T)
data.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.inv.data),header=T,sep="\t", stringsAsFactors=F)
prev.inv          <- read.delim(paste0(dir.inv.data, file.prefix, file.prev),header=T,sep="\t", stringsAsFactors=F)

# set date for file names
d <- Sys.Date()    # e.g. 2021-12-17

# Load functions ####

source("ml_model_functions.r")
source("stat_model_functions.r")
source("plot_functions.r")
source("utilities.r")

# Setup options ####

# Set if we want to fit models to whole dataset or perform cross-validation (CV)
CV <- T # Cross-Validation
dl <- F # Data Leakage

# Set number of cores
n.cores.splits <-  3 # a core for each split, 3 in our case
n.cores.stat.models <- 2 # a core for each stat model 2 in our case (UF0, and CF0)
# Settings Stat models 
# Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models,
# also make sure the cross-validation (CV) flag is set correctly
sampsize <- 1000 #10000 #I think this needs to be an even number for some reason (stan error)
n.chain  <- 2 #2

# Select taxa
all.taxa <- T
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa

# Run the script on the server
server <- F

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

prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, 
                               env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                               BDM = BDM, dl = dl, CV = CV)
data <- prepro.data$data
if(CV){ splits <- prepro.data$splits ; rem.taxa <- prepro.data$rem.taxa }
list.taxa <- prepro.data$list.taxa # list taxa after data pre-processing 
centered.data <- prepro.data$centered.data
centered.data.factors <- prepro.data$centered.data.factors
normalization.data <- prepro.data$normalization.data

remove(prepro.data)

# Select taxa ####

cind.taxa <- which(grepl("Occurrence.",colnames(data)))
list.taxa.full <- colnames(data)[cind.taxa] # list of all taxa, before removing these too unbalanced or with too many NA
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

## ---- APPLY MODELS ----

# Null model

null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]

# Neural Networks #### 

# Learning rate
learning.rate <- 0.01

# Number of epochs
# num.epochs <- 50

# Batch size
batch.size <-  64

grid.hyperparam <- expand.grid(layers = c(3,5), units = c(32, 64), 
                               act.fct = c("tanh", "leakyrelu"), no.epo = c(50,100))
no.hyperparam <- nrow(grid.hyperparam)
list.hyper.param <- vector("list", no.hyperparam)
for (n in 1:no.hyperparam) {
  list.hyper.param[[n]] <- grid.hyperparam[n,]
  names(list.hyper.param)[n] <- paste(paste0(grid.hyperparam[n,], c("L", "U", "FCT", "epo")), collapse = "")
}

# list.ann <- list()
# names(list.ann)
names(list.hyper.param) <- paste("ANN_", names(list.hyper.param), sep = "")

source("ann_model_functions.r")
source("utilities.r")

if(CV){
  ann.outputs <- lapply(centered.data, function(split){
    lapply(list.hyper.param, FUN = build_and_train_model, split = split,
           env.fact = env.fact,
           list.taxa = list.taxa,
           learning.rate = learning.rate,
           # num.epochs = num.epochs,
           batch.size = batch.size,
           CV = CV)
  })
} else {
  ann.outputs <- lapply(list.hyper.param, FUN = build_and_train_model, split = centered.data,
                        env.fact = env.fact,
                        list.taxa = list.taxa,
                        learning.rate = learning.rate,
                        # num.epochs = num.epochs,
                        batch.size = batch.size,
                        CV = CV)
}

if(CV){
  list.splits <- names(ann.outputs)
  no.splits <- length(list.splits)
}

no.ann <- length(list.hyper.param)
list.models <- names(list.hyper.param)
names(list.models) <- rainbow(no.ann) # assign colors
no.models <- length(list.models)

info.file.ann.name <-  paste0("ANN_model_",
                             file.prefix, 
                             no.ann, "ann_",
                             no.taxa, "taxa_", 
                             ifelse(CV, "CV_", "FIT_"),
                             ifelse(dl, "DL_", "no_DL_"))

file.name <- paste0(dir.models.output, info.file.ann.name, ".rds")
cat(file.name)
saveRDS(ann.outputs, file = file.name)

info.file.name <- info.file.ann.name


if(CV){ 
  
  # Merge all CV outputs in one
  outputs.cv <- ann.outputs
  # for (s in list.splits) {
  #   outputs.cv[[s]][[list.stat.mod[1]]] <- stat.outputs.transformed[[1]][[s]]
  #   outputs.cv[[s]][[list.stat.mod[2]]] <- stat.outputs.transformed[[2]][[s]]
  #   outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv[[s]])
  # }
  
  # Make final outputs as list
  outputs <- make.final.outputs.cv(outputs.cv = outputs.cv, list.models = list.models, list.taxa = list.taxa)
  
  # Make final outputs as tables
  df.cv <- make.df.outputs(outputs = outputs.cv, list.models = list.models, 
                           list.taxa = list.taxa, list.splits = list.splits,
                           null.model = null.model, prev.inv = prev.inv, CV = CV)
  df.perf.cv <- df.cv$`Table performance CV`
  df.perf <- df.cv$`Table performance`
  remove(df.cv)
} else {
  # Make final outputs as list
  outputs <-  ann.outputs
  
  # Make final outputs as tables
  df.perf <- make.df.outputs(outputs = outputs, list.models = list.models, 
                             list.taxa = list.taxa, list.splits = list.splits,
                             null.model = null.model, prev.inv = prev.inv, CV = CV)
}

# Table with performance
if(CV){
  list.plots.cv <- plot.df.perf(df.perf = df.perf.cv, list.models = list.models, list.taxa = list.taxa, CV)
  list.plots <- plot.df.perf(df.perf = df.perf, list.models = list.models, list.taxa = list.taxa, CV)
  list.plots <- append(list.plots.cv, list.plots)
} else {
  list.plots <- plot.df.perf(df.perf = df.perf, list.models = list.models, list.taxa = list.taxa, CV)
}

name <- "TablesPerf"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 12, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Performance against prevalence and boxplots
list.plots <- model.comparison(df.perf = df.perf, list.models = list.models, CV = CV)
name <- "ModelsCompar"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 12, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)





ann.outputs.fit <- readRDS(paste0(dir.models.output, "ANN_model_All_3ann_126taxa_FIT_no_DL_.rds"))

summary(data.env$BIOGEO)
# Alpennordflanke  Alpensüdflanke Jura      Mittelland    Zentralalpen 
# 641             146             419       1549          260
select.biogeo <- "Zentralalpen"
select.ann <- names(ann.outputs.fit)[1]
select.samples <- data.env[which(data.env$BIOGEO == select.biogeo), "SampId"]
select.samples <- select.samples[select.samples %in% centered.data[[1]][, "SampId"]]
no.samples <- length(select.samples)
pred.df <- data.frame(matrix(nrow = no.samples, ncol = no.taxa))
colnames(pred.df) <- list.taxa
rownames(pred.df) <- select.samples
for (j in list.taxa) {
  temp.df <- bind_cols(centered.data[[1]], ann.outputs.fit[[select.ann]][[j]][["Prediction probabilities training set"]])
  temp.df <- temp.df[,c("SampId", "present")]
  for (s in select.samples) {
    # Recover sample ID
    pred.df[s,j] <- temp.df[which(temp.df$SampId == s),"present"]
  }
}
pred.df$SampId <- select.samples

library("reshape2")
library("ggplot2")
melted.df <- melt(pred.df, id = "SampId")
title <- paste("Predictions of", select.ann, "\nin region", select.biogeo, "(with", no.samples, "samples)")
p <- ggplot(data = melted.df, aes(x = variable, y = SampId, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0.5, mid ="grey70", low = "#c2141b", high = "#007139",
                       limits = c(0, 1)) +
  scale_x_discrete(limits = list.taxa) +
  labs(title = title, 
       x = "", y = "", fill = "Predicted probability \n of occurence") +
  theme(plot.title = element_text(hjust = 0.5, colour = "black"), 
        axis.title.x = element_text(face="bold", colour="darkgreen", size = 2),
        axis.text.x = element_text(angle=90),
        axis.text.y = element_blank(),
        legend.title = element_text(face="bold", colour="brown", size = 10))  
  # geom_text(aes(x = Taxa, y = variable, label = round(value, 2)),
  #           color = "black", fontface = "bold", size = size.val)

file.name <- paste0("TablesPred_", select.ann, "_", select.biogeo,".pdf")
cat(file.name)
dir.plots.output <- "C:/Users/cholleem/Documents/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/"
pdf(paste0(dir.plots.output, file.name), paper = 'special', width = 20, onefile = TRUE)
print(p)
dev.off()
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


