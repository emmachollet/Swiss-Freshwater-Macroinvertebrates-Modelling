## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- TRAINING ML MODELS ON INVERTEBRATE AND ENVIRONMENTAL DATA ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 21, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Libraries, data and functions ----

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
if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") }
if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm

# install.packages(c('caret', 'skimr', 'RANN', 'randomForest', 'fastAdaboost', 'gbm', 'xgboost', 'caretEnsemble', 'C50', 'earth'))
# install.packages("ggpubr")
# library(dplyr) # to sort, join, merge data
# library(ggplot2) # to do nice plots
# library(caret) # comprehensive framework to build machine learning models
# library(skimr) # to show key descriptive stats
# library(RANN)  # required for knnInpute
# library(e1071) # to do a recursive feature elimination
# library(caretEnsemble) # to do ensemble prediction with caret

# Check and set working directory
getwd() # show working directory
# setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts") # set the working directory to this folder

# Free workspace
rm(list=ls())
graphics.off()

# Define directory and files
dir.env.data      <- "../../Data/Processed data/Environmental data/"
dir.inv.data      <- "../../Data/Processed data/Invertebrate data/"
dir.output        <- "../Plots/"

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
BDM <- TRUE

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

# Load functions
source("model_functions.r")
source("plot_functions.r")

## ---- Construct datasets  ----

# Replace "0" and "1" by "absent" and "present" and convert them to factors
cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))

for (i in cind.taxa ) {
    data.inv[which(data.inv[,i] == 0),i] <- "absent"
    data.inv[which(data.inv[,i] == 1),i] <- "present"
    data.inv[,i] = as.factor(data.inv[,i])
}

# # Convert environmental character columns to factors
# data.env[sapply(data.env, is.character)] <- lapply(data.env[sapply(data.env, is.character)], as.factor)

# Construct main dataset (with inv and env)
data <- data.env %>%
    left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data)

# Select environmental factors for prediction

# Replace "_" and " " by "." in colnames to be consistent
colnames(rank.env) <- gsub("_", ".", colnames(rank.env))
colnames(rank.env) <- gsub(" ", ".", colnames(rank.env))

# Sort the columns in 4 categories: the sample/site information = 3,
# the environmental factors to keep in priority = 2, keep to explore = 1, exclude = 0
info    <- colnames(rank.env)[which(rank.env[1,] == 3)]
prio    <- colnames(rank.env)[which(rank.env[1,] == 2)]
explo   <- colnames(rank.env)[which(rank.env[1,] == 1)]
excl    <- colnames(rank.env)[which(rank.env[1,] == 0)]

a.lot = F     # set to FALSE if it's for exploration (few taxa and env. fact.)
              # set to TRUE for a proper simulation (all important taxa and env. fact.)
# salut
if (a.lot == T){
    
    env.fact <- prio
    # env.fact <- c(prio,explo)
    
} else if (a.lot == F){
    
    env.fact <- c("temperature",
                    "velocity",
                    "cow.density",
                    "IAR",
                    "urban.area",
                    "FRI"#,
                    # "WALD.ANT",
                    # "width.variability",
                    # "bed.modification",
                    # "morphology",
                    # "A.EDO",
                    # "F.EDO",
                    # "ARA.fraction"#,
                    # "agri.land",
                    # "Slope",
                    # "fields",
                    # "saprobic.cond"#,
                    # "normcov.mobile.blocks",
                    # "normcov.coarse.inorganic.sediments",
                    # "normcov.gravel",
                    # "normcov.sand.silt",
                    # "normcov.fine.sediments")
                    )
}

# Remove env. fact. missing in data.env
env.fact <- env.fact[which(env.fact %in% colnames(data.env))]
# env.fact[-which(env.fact %in% colnames(data.env))] # which selected factors are not in data.env ?
# [1] "A.EDO" "F.EDO" are missing

# Select taxa for prediction

if (a.lot == T){
    
    list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25), 
                                "Occurrence.taxa"] # Select with prevalence percentage between 25 and 75%

    
} else if (a.lot == F){
    
    # list.taxa <- c("Occurrence.Gammaridae", "Occurrence.Heptageniidae")
    list.taxa       <- prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55), 
                                "Occurrence.taxa"] # Select only few taxa
}


# Construct training and testing datasets (randomly or according to an env. fact.)
ratio <- 1 # set a ratio for training dataset (can be 1)
split.var <- "temperature" # data splitted according to this variable (default is "random")
splitted.data <- split.data(data = data, training.ratio = ratio, variable = split.var)

# Assemble information to insert in file names
no.env.fact <- length(env.fact)
no.taxa <- length(list.taxa)
percentage.train.set <- ratio * 100
info.file.name <- paste0(file.prefix, d,
                         no.taxa, "taxa_", 
                         no.env.fact, "envfact_", 
                         "trainset", percentage.train.set, 
                         if( ratio != 1) {split.var}, 
                         "_")

# Summary of prevalence of chosen taxa
for ( i in 1:no.taxa){
    cat("Summary of absence and presence for", list.taxa[i], ":", summary(data[, list.taxa[i]]), "\n")
}

## ---- Plot taxa vs env. fact. (data visualization) ----

# Compute plots
list.plots <- plot.data.envvstax(data = data, env.fact = env.fact, list.taxa = list.taxa)

# Print the plots in a pdf file
file.name <- "Data_EnvFactvsTax.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)


## ---- Feature selection (rfe) ----

# We might need a rigorous way to determine the important variables first before feeding them to the ML algorithm.
# RFE works in 3 broad steps:
# Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
# Step 2: Keeping priority to the most important variables, iterate through by building models of given subset sizes, that is, 
#         subgroups of most important predictors determined from step 1. Ranking of the predictors is recalculated in each iteration.
# Step 3: The model performances are compared across different subset sizes to arrive at the optimal number and list of final predictors.

# ptm <- proc.time() # to calculate time of calculation
# 
# ctrl <- rfeControl(functions = rfFuncs, # based on random forest algorithm
#                    method = "repeatedcv", # method of cross-validation
#                    repeats = 2,
#                    verbose = FALSE)
# 
# subsets <- c(1:5, round(no.env.fact/2), round(3/4*no.env.fact)) # select size of subsets of env. fact. to be tested
# 
# # Make a list with the different predictive performance for each taxon in list.taxa
# rfe.result <- vector(mode = 'list', length = no.taxa)
# names(rfe.result) <- list.taxa
# 
# for (j in 1:no.taxa){
#   
#   temp.data <- na.omit(data[, c(list.taxa[j],env.fact)])
#   
#   rfe.result[[j]] <- rfe(x=temp.data[, env.fact], y=temp.data[, list.taxa[j]], # receives the output of the rfeControl() as values
#                          sizes = subsets, # what all model sizes (the number of most important features) the rfe should consider.
#                          rfeControl = ctrl)
#   
# }
# 
# print(paste(no.env.fact, "rfe calculation time:"))
# print(proc.time()-ptm)

# WORKSPACE: BDM, 55 env.fact, 29 taxa : 1:30 hours, saved 12.08.2021

## ---- Apply models ----
ptm <- proc.time() # to calculate time of simulation

# Apply machine learning (ML) models
# Select models to apply (! their packages have to be installed first)
list.algo <- c('glm', # Random Forest
               'rf',
               # 'adaboost',
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # 'xgbDART',
               #'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
               'svmRadial')

# if ( !require("elmNN") ) { install.packages("elmNN"); library("elmNN") }
# if ( !require("arm") ) { install.packages("arm"); library("arm") } 

# Assemble information to insert in file names
no.algo <- length(list.algo)
info.file.name <- paste0(info.file.name, no.algo, "algo_")

# Define the training control
train.control <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 3,                      # number of folds
  # repeats = 1,                   # for repeated k-fold cross-validation 'repeatedcv' only: the number of complete sets of folds to compute
  # savePredictions = 'final',     # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  #summaryFunction = twoClassSummary, # default metric function choosen (accuracy, ROC)
  summaryFunction = stand.dev ,
  selectionFunction = lowest       # we want to minimize the metric
)

# Make a list to store the outputs of each model
outputs <- vector(mode = 'list', length = no.algo)
names(outputs) <- list.algo

# Apply the ML models to list.taxa and env.fact
for(k in 1:no.algo){
  outputs[[k]] <- apply.ml.model(splitted.data = splitted.data, list.taxa = list.taxa,
                                        env.fact = env.fact, algorithm = list.algo[k], train.control = train.control)
}

# "Apply" null model
null.model <- apply.null.model(data = data, list.taxa = list.taxa, prev.inv = prev.inv)

print(paste("Simulation time of different models ", info.file.name))
print(proc.time()-ptm)

# WORKSPACE: BDM, 55 env.fact, 29 taxa, 4 algo : 18 hours, saved 18.08.2021
# BDM, 6 env.fact, 2 taxa, 4 algo : 3 min
# WORKSPACE: BDM, 55 env.fact, 29 taxa, 1 algo (rf) : 15 hours, saved 11.08.2021
# WORKSPACE: BDM, 9 env.fact, 9 taxa, 4 algo : 1 hours, saved 01.09.2021
# Don't forget to reload the functions if we upload old workspace
# source("model_functions.r")
# source("plot_functions.r")

## ---- Plot models comparison ----

# Compare models by resampling
ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- model.comparison(outputs = outputs, null.model = null.model, list.algo = list.algo, list.taxa = list.taxa, prev.inv = prev.inv)

# Print the plots in a pdf file
file.name <- "ModelsCompar.pdf"
print.pdf.plots(list.plots = list.plots, width = 9, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)


## ---- Plot variable importance ----

ptm <- proc.time() # to calculate time of simulation

list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)

file.name <- "VarImp.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

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
          outputs[[l]][[j]][["Confusion matrix"]][["overall"]][["Accuracy"]], "\n")
  }
}

## ---- Plot PDP ----

ptm <- proc.time() # to calculate time of simulation

# PDP of one model
# list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
#                       list.taxa = list.taxa, env.fact = env.fact)

#file.name <- "PDP.pdf"
#print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

# PDP of all models
list.plots <- plot.pdp(outputs = outputs, list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "allPDP.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot single predictor ICE ----

ptm <- proc.time() # to calculate time of simulation

# PDP of one model
list.plots <- plot.ice(outputs = outputs, algo = list.algo[1], list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "ICE.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

## ---- Plot PDP (multiple predictors) ----

# We just make one example because it's computationally heavy

ptm <- proc.time() # to calculate time of simulation

# Multpiple predictors PDP (of one model) for now just for 1 algo and 1 taxa
list.plots <- plot.mult.pred.pdp(outputs = outputs, list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

file.name <- "multpredPDP.pdf"
print.pdf.plots(list.plots = list.plots, width = 17, dir.output = dir.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

print("hey")
# Stop execution if there is no (prediction on) testing set
# stopifnot(ratio != 1)
if( ratio == 1 ){
  exit
}
print("salut")

## ---- Plot map prediction ----

ptm <- proc.time() # to calculate time of pdf production

map.inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

# make a list with all plots and plot them in a pdf
list.plots <- map.ml.pred.taxa(inputs = map.inputs, outputs = outputs,
                               list.taxa = list.taxa, list.algo = list.algo)

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
