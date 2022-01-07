## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- PhD Project ---
## 
##                          --- December 23, 2021 -- Happy Christmas 
##
## --- Emma Chollet, Jonas Wydler, Andreas Scheidegger and Nele Schuwirth ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# next two lines are for jonas, dont run these
# setwd("Q:/Abteilungsprojekte/siam/Jonas Wydler/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/R scripts")
#.libPaths("C:/Program Files/R/R-4.1.1/library")

## ---- PACKAGES, DATA & FCTS ----

# Load libraries ####

if ( !require("parallel") ) { install.packages("parallel"); library("parallel") } # need to run things in parallel

# Data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") } # to sort, join, merge data
if ( !require("splitTools") ) { install.packages("splitTools"); library("splitTools") } # to split the data

# Plots
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
#if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") } # to arrange multiple plots on a page
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page
if ( !require("cowplot") ) { install.packages("cowplot"); library("cowplot") } # to arrange multiple plots on a page
if ( !require("pdp") ) { install.packages("pdp"); library("pdp") } # to plot partial dependance plots
if ( !require("gt") ) { install.packages("gt"); library("gt") } # to plot nice tables
if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") } # to plot nice tables
if ( !require("viridis")) {install.packages("viridis", repos="http://cloud.r-project.org"); library("viridis")} # to do even nicer plots
#if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map

# Stat model
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") } # to read layers for map

# ML algorithms
if ( !require("caret") ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models
if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") } # to run generalized additive model (GAM) algorithm
if ( !require("gam") ) { install.packages("gam"); library("gam") } # to run generalized additive model (GAM) algorithm
# if ( !require("fastAdaboost") ) { install.packages("fastAdaboost"); library("fastAdaboost") } # to run adaboost ml algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") } # to run support vector machine (svm) algorithm
# if ( !require("earth") ) { install.packages("earth"); library("earth") } # to run MARS ml algorithm
if ( !require("randomForest") ) { install.packages("randomForest"); library("randomForest") } # to run random forest (RF)
# if ( !require("randomForestSRC") ) { install.packages("randomForestSRC"); library("randomForestSRC") } # to run RF and additional features

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
all.taxa <- F
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa

# Run the script on the server
server <- T

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

# DISCUSS TAXA DELETING WITH NELE ? ####
# JW: THERE IS ALSO ONE MORE SPECIES THAT GETS DROPPED IN THE CENTERING (ENDS UP AT 126 rather tha 127, move to
#data prep as well)

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

# # Summary of prevalence of chosen taxa
# for ( i in 1:no.taxa){
#     cat("Summary of absence, presence and NA for", list.taxa[i], ":", summary(data.inv[, list.taxa[i]]), "\n")
# }

# Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c( "#030AE8" = 'glm', # Random Forest
               # "#048504" = 'bam', # Generalized Additive Model using splines
                "#948B8B" = 'gamSpline', #
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # "#A84E05" = 'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
                "#DB1111" = 'svmRadial', # Support Vector Machine
               # "#DB1111" = 'svmPoly', # Support Vector Machine
               # "#DB1111" = 'svmLinear', # Support Vector Machine
                "#790FBF" = 'rf') # Random Forest

no.algo <- length(list.algo)
                                            
## ---- APPLY MODELS ----

# Null model

null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]

# Statistical models ####

ptm <- proc.time() # to calculate time of simulation
#file.name <- paste0(dir.models.output, "Stat_model_100iterations_corr_22taxa_CV_no_DL.rds") #to test

comm.corr.options <- c(T,F)
names(comm.corr.options) <- c("CF0", "UF0")
stat.outputs <- mclapply(comm.corr.options, mc.cores = 1, function(comm.corr){
  
  comm.corr <- comm.corr.options[[1]]
  info.file.stat.name <- paste0("Stat_model_",
                                file.prefix,
                                no.taxa.full, "taxa_", 
                                sampsize,"iterations_",
                                ifelse(comm.corr,"CF0_","FF0_"),
                                ifelse(CV, "CV_", "FIT_"),
                                # else(dl, "DL_", "no_DL_"),
                                "no_DL_") # for now with just apply it without DL
  
  file.name <- paste0(dir.models.output, info.file.stat.name, ".rds")
  cat(file.name)
  
  # If the file with the output already exist, just read it
  if (file.exists(file.name) == T){
      
      if(exists("stat.output") == F){
          cat("File with statistical model output already exists, we read it from", file.name, "and save it in object 'stat.output'")
          stat.output <- readRDS(file = file.name)
          
          }
      else{
          cat("List with statistical model output already exists as object 'stat.output' in this environment.")
      }
  } else {
      
      cat("No statistical model output exist yet, we produce it and save it in", file.name)
      
      if(CV == T){
  
        stat.output <- mclapply(centered.splits, mc.cores = n.cores, FUN = stat_mod_cv, CV, comm.corr, sampsize, n.chain)
        
        cat("Saving output of statistical models in", file.name)
        saveRDS(stat.output, file = file.name, version = 2) #version two here is to ensure compatibility across R versions
        
      } else {
        # apply temporary on training data of Split1, later take whole dataset centered etc
  
        stat.output <- stat_mod_cv(data.splits = centered.data, CV, comm.corr, sampsize, n.chain = n.chain)
        cat("Saving output of statistical models in", file.name)
        saveRDS(stat.output, file = file.name, version = 2)
        }
      
  }
  return(stat.output)
})
print(paste("Simulation time of statistical model "))
print(proc.time()-ptm)

#Process output from stat models to fit structure of ml models (makes plotting easier)
stat.outputs.transformed <- transfrom.stat.outputs(CV, stat.outputs)

# Machine Learning models ####

ptm <- proc.time() # to calculate time of simulation

info.file.ml.name <-  paste0("ML_model_",
                             file.prefix, 
                             no.algo, "algo_",
                             no.taxa, "taxa_", 
                             ifelse(CV, "CV_", "FIT_"),
                             ifelse(dl, "DL_", "no_DL_"))


file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
# file.name <- paste0(dir.models.output, "glm_gamSpline_svmRadial_rf_22taxa_FIT.rds")
# file.name <- paste0(dir.models.output, "ML_model_All_4algo_126taxa_FIT_no_DL_.rds")

cat(file.name)

# # If the file with the outputs already exist, just read it
# if (file.exists(file.name) == T ){
#     
#     if(exists("ml.outputs") == F){
#         cat("File with ML outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
#         ml.outputs <- readRDS(file = file.name)
#         }
#     else{
#         cat("List with ML outputs already exists as object 'ml.outputs.fit' in this environment.")
#     }
# 
# } else {
#

# File with lot of taxa already exist, we read it and subselect taxa later

info1 <- paste0(file.prefix, no.algo, "algo_")
info2 <- paste0(ifelse(CV, "CV_", "FIT_"), ifelse(dl, "DL_", "no_DL_"))
list.ml.files <- list.files(path = dir.models.output)
file.name.temp <- paste0(dir.models.output, list.ml.files[which(grepl(info1, list.ml.files) & grepl(info2, list.ml.files))])

if( file.exists(file.name.temp) == T ){
    # & askYesNo(paste("Another file named", file.name.temp, "already exists. Should we read it and subselect taxa later ? If no, we run the algorithms for selected information."), default = T) ){

    cat("Reading", file.name.temp, "and saving it in object 'ml.outputs.fit'")
    ml.outputs <- readRDS(file = file.name.temp)

    } else {

        cat("No ML outputs exist yet, we produce it and save it in", file.name)

    if(CV == T){
        
        if(server == T){
            # Compute three splits in paralel (should be run on the server)
            outputs.cv <- mclapply(centered.data.factors, mc.cores = n.cores.splits,
                                   FUN = apply.ml.model, list.algo, list.taxa, env.fact, prev.inv = prev.inv)
        } else {
            # Compute one split after the other
            outputs.cv <- lapply(centered.data.factors, FUN = apply.ml.model, list.algo, list.taxa, env.fact, CV, prev.inv = prev.inv)
        }
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(outputs.cv, file = file.name, version = 2)
    
        } else {
        
        splitted.data <- list("Training data" =  centered.data.factors[[1]], "Testing data" = data.frame())
        ml.outputs <- apply.ml.model(splitted.data = splitted.data, list.algo = list.algo, list.taxa = list.taxa,
                                      env.fact = env.fact, CV = F, prev.inv = prev.inv)
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(ml.outputs, file = file.name, version = 2)
        
        }
}
# }

print(paste("Simulation time of different models ", info.file.ml.name))
print(proc.time()-ptm)

if(!server){

# Neural Networks ####

# read output from ann_model.r, hopefully ...

# file.name <- paste0(dir.models.output, "test_ANNoutputs.rds")
# ann.outputs <- readRDS(file = file.name)


# Training GAM bam for 6 taxa: 5 hours

# source("ml_model_functions.r")
# source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# Clean outputs ####

# Make mean over splits for final cross validation
if(CV == T){
    
    outputs <- vector(mode = "list", length = length(list.algo))
    names(outputs) <- list.algo
    
    out <- c("Observation", #1
             "Prediction factors", #2 
             "Prediction probabilities", #3 
             "Likelihood", #4
             "Performance", #5
             "Performance splits") #6
    output.names <- paste(out, "testing set")
    
    for (l in 1:no.algo) {
        
        temp.list.taxa <- vector(mode = "list", length = length(list.taxa))
        names(temp.list.taxa) <- list.taxa
        
        for(j in 1:no.taxa){
            
            temp.output <- vector(mode = "list", length = length(output.names))
            names(temp.output) <- output.names
            
            for (m in output.names[c(1,3)]) {
                temp.output[[m]] <- bind_rows(outputs.cv[[1]][[l]][[j]][[m]],
                                              outputs.cv[[2]][[l]][[j]][[m]],
                                              outputs.cv[[3]][[l]][[j]][[m]], .id = "Split")
                }
            for (m in output.names[c(2,4)]) {
                temp.output[[m]] <- c(outputs.cv[[1]][[l]][[j]][[m]],
                                      outputs.cv[[2]][[l]][[j]][[m]],
                                      outputs.cv[[3]][[l]][[j]][[m]])
                }
            
            temp.vect <- vector(mode ="numeric", length = length(outputs.cv)) 
            names(temp.vect) <- names(outputs.cv)
            for (n in 1:length(outputs.cv)) {
              #n = 1
                temp.vect[n] <- outputs.cv[[n]][[l]][[j]][["Performance testing set"]]
            }
            
            temp.output[[5]] <- mean(temp.vect)
            temp.output[[6]] <- temp.vect
            
            temp.list.taxa[[j]] <- temp.output
        }
    
        outputs[[l]] <- temp.list.taxa
    }
    # # Add output of stat model
    # # outputs.cv
    # temp.list.stat <- vector(mode = "list", length = length(stat.cv.res_table$Taxon))
    # names(temp.list.stat) <- stat.cv.res_table$Taxon
    # 
    # for( j in 1:length(stat.cv.res_table$Taxon)){
    #   
    #   temp.vect <- vector(mode ="numeric", length = length(stat.cv.res_table$Taxon))
    #   for (n in 1:length(stat.cv.res_table$Taxon)) {
    #     temp.list.stat[n] <- stat.cv.res_table$performance[[n]]
    #   }
    # }
    # UF0 <- as.list(temp.list.stat)
    # UF0 <- UF0[list.taxa]
    # 
    # outputs[[5]] <- UF0
    # names(outputs)[[5]] <- "UF0"
    
    # # Add output of ANN
    # for (a in 1:length(ann.outputs)) {
    #     outputs[[4 + a]] <- ann.outputs[[a]]
    #     names(outputs)[4 + a] <- names(ann.outputs)[a]
    #     
    # }
}


# Make list and info for all models

list.models <- c(# "Stat", 
                 # "ANN",
                 list.algo)
no.models <- length(list.models)

info.file.name <- paste0(file.prefix, 
                         no.models, "models_",
                         # d, # don't need to include the date
                         no.taxa, "taxa_", 
                         # no.env.fact, "envfact_",
                         ifelse(CV, "CV_", "FIT_"),
                         ifelse(dl, "DL_", "no_DL_"),
                         # "trainset", percentage.train.set, 
                         # if( ratio != 1) {split.var}, 
                         "")


# Select taxa in outputs for analysis

# Make a list to store the outputs of each model
# outputs <- vector(mode = 'list', length = no.models)
# names(outputs) <- list.models

# test <- lapply(outputs.full[l], "[[", list.taxa)
# 
# for (l in list.models) {
#     lind <- which(names(outputs.full[l][[]]) == list.taxa)
#     outputs[l]  <- outputs.full[l][list.taxa]
# }

# AND CHANGE THE PLOTS DEPENDING ON IF THEY NEED ML OR ALL MODELS

## ---- PLOTS ----
#Add stat models to the list of algos
outputs <- append(ml.outputs.fit, stat.outputs.transformed)
list.algo <- append(list.algo, list("#2CA25F" = 'CF0', "#99D8C9" = 'UF0'))
no.algo <- length(list.algo)
# list.algo <- c(list.algo, "blue" = names(ann.outputs)[1],"red" = names(ann.outputs)[2], "grey" = names(ann.outputs)[3])
# no.algo <- length(list.algo)
# 
# # Write information for file names
# # percentage.train.set <- ratio * 100
# info.file.name <- paste0(file.prefix, 
#                          # d, # don't need to include the date
#                          no.taxa, "taxa_", 
#                          # no.env.fact, "envfact_",
#                          no.algo, "algo_",
#                          ifelse(CV, "CV_", "FIT_"),
#                          ifelse(dl, "DL_", "no_DL_"),
#                          # "trainset", percentage.train.set, 
#                          # if( ratio != 1) {split.var}, 
#                          "")

# PROBLEM WITH PERF GAM ####

for (j in 1:no.taxa) {
    for (l in 1:no.algo) {
        perf <- outputs[[l]][[list.taxa[j]]][["Performance training set"]]
        outputs[[l]][[list.taxa[j]]][["Performance training set"]] <- ifelse(perf > 1.5, Inf, perf)
    }
}

# Models comparison ####

# table with perf.

if(CV == T){
    # maybe look at perf. in each split
} else {
    file.name <- "TableFitPerf.pdf"
    pdf(paste0(dir.plots.output, info.file.name, file.name), paper = 'special', width = 12, height = 9, onefile = TRUE)
    
    temp.df <- data.frame(matrix(ncol = no.taxa, nrow = no.algo))
    colnames(temp.df) <- list.taxa
    rownames(temp.df) <- list.algo
    for (j in 1:no.taxa) {
        for (l in 1:no.algo) {
            perf <- outputs[[l]][[list.taxa[j]]][["Performance training set"]]
            temp.df[l,j] <- round(outputs[[l]][[list.taxa[j]]][["Performance training set"]], digits = 2)
        }
    }
    temp.df$mean.perf <- rowMeans(temp.df)
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
         main = "Quality of fit across ML algorithms and taxa"
    )
    axis(1, at=seq(1:ncol(temp.df)+1), labels = FALSE)
    text(seq(1:ncol(temp.df)+1), par("usr")[4] + 0.15, srt = 50, 
         labels = colnames(temp.df), adj= 0, cex = 0.5, xpd = T)
    
    dev.off()
}

# saved.outputs <- outputs

# outputs <- outputs2algo[[1]]

ptm <- proc.time() # to calculate time of pdf production

# TO BE FIXED EARLIER ####
# list.algo <- c(list.algo, "green" = "UF0")
# source("plot_functions.r")

# Compute plots
list.plots <- model.comparison(outputs = outputs, null.model = null.model, 
                               list.algo = list.algo, list.taxa = list.taxa[1:5], prev.inv = prev.inv, CV = CV)

# plot it
## THIS SHOULD BE FIXED ####
# list.plots <- model.comparison.cv(outputs = outputs, outputs.cv = outputs.cv, null.model = null.model, list.algo = list.algo, list.taxa = list.taxa, prev.inv = prev.inv)

name <- "ModelsCompar"
file.name <- paste0(name, ".pdf")

print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# Performance vs hyperparameters ####

ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- plot.perf.hyperparam(outputs = outputs, 
                                   list.algo = list.algo, # GLM algo doesn't have hyperparameters
                                   list.taxa = list.taxa)

name <- "PerfvsHyperparam"

# Print a pdf file
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()

print(paste(file.name, "printing:"))
print(proc.time()-ptm)


# Variables importance ####

ptm <- proc.time() # to calculate time of simulation

list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)

name <- "VarImp"

# Print a pdf file
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 10, height = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# Print a jpeg file
file.name <- paste0(name, ".jpg")
jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
print(list.plots[[1]])
dev.off()

# print directly table with variable importance for each algo and taxa (too complicated to put in a fct)
file.name <- "TableVarImp.pdf"
pdf(paste0(dir.plots.output, info.file.name, file.name), paper = 'special', width = 12, height = 9, onefile = TRUE)
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

# Print accuracy

for(l in 1:no.algo){
    for(j in 1:no.taxa){
        # cat("Accuracy on training for",list.taxa[j], "is :",
        #    output[[j]][["Trained model"]][["resample"]][["Accuracy"]], "\n")
        cat("Accuracy on prediction for",list.taxa[j],"with", list.algo[l], "is :",
            outputs[[l]][[j]][["Confusion matrix testing set"]][["overall"]][["Accuracy"]], "\n")
    }
}

# PDP don't work ####

source("plot_functions.r")

ptm <- proc.time() # to calculate time of simulation

# PDP of one model
list.plots <- plot.pdp(outputs = outputs[[1]], algo = "rf", list.algo = list.algo,
                       list.taxa = list.taxa, env.fact = env.fact)

ptm <- proc.time() # to calculate time of simulation

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
