## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- PhD Project ---
## 
##                          --- December 10, 2021 -- 
##
## --- Emma Chollet, Jonas Wydler, Andreas Scheidegger and Nele Schuwirth ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# next two lines are for jonas, dont run these
# setwd("Q:/Abteilungsprojekte/siam/Jonas Wydler/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/R scripts")
#.libPaths("C:/Program Files/R/R-4.1.1/library")

## ---- PACKAGES, DATA & FCTS ----

# Load libraries ####

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
if ( !require("viridis")) {install.packages("viridis", repos="http://cloud.r-project.org"); library("viridis")} # to do even nicer plots
if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map

# stat model
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") } # to read layers for map

# ml algorithms
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

# set if we want to fit models to whole dataset or perform cross-validation (CV)
CV <- F

# set date for file names
d <- Sys.Date()    # e.g. 2021-12-17

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

no.env.fact <- length(env.fact)

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
                         # no.env.fact, "envfact_",
                         no.algo, "algo_",
                         ifelse(CV, "CV_", "FIT_"),
                         # "trainset", percentage.train.set, 
                         # if( ratio != 1) {split.var}, 
                         "")

# Pre-process data ####

# Construct main dataset (with inv and env)
data <- data.env[, c("SiteId", "SampId", "X", "Y", env.fact)] %>%
    left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data)

# drop rows with incomplete influence factorls
ind <- !apply(is.na(data[,env.fact]),1,FUN=any)
ind <- ifelse(is.na(ind),FALSE,ind)
data <- data[ind,]
print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))

# Split (and save) data ####

# Split for ml purpose (e.g. 80% training, 20% testing)
# ratio <- 0.8 # set a ratio for training dataset (can be 1)
# split.var <- "random" # data splitted according to this variable (default is "random")
# splitted.data <- split.data(data = data, training.ratio = ratio, variable = split.var)

# Split for CV
file.name <- paste0(dir.workspace,"SplitsForCV_101221.rds")

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

### ----Investigate splits----

# mean.train1 <- as.numeric(lapply(centered.splits[[1]][[1]][,env.fact], FUN = mean))
# mean.train2 <- as.numeric(lapply(centered.splits[[2]][[1]][,env.fact], FUN = mean))
# mean.train3 <- as.numeric(lapply(centered.splits[[3]][[1]][,env.fact], FUN = mean))
# 
# sd.train1 <- as.numeric(lapply(centered.splits[[1]][[1]][,env.fact], FUN = sd))
# 
# dif1 <- 1-(mean.train1/mean.train2)
# dif2 <- 1-(mean.train2/mean.train3)
# testest <- as.data.frame(rbind(dif1,dif2))
# 
# colnames(testest) <- env.fact
# testest.table <- lapply(testest, mean, 1)
# 
# 
# map.input <- map.inputs(dir.env.data, data.env)
# 
# 
# g1 <- ggplot()
# g1 <- g1 + geom_sf(data = map.input$ch, fill="#E8E8E8", color="black")
# g1 <- g1 + geom_point(data = splits[[1]][[1]], aes(x=X, y=Y), color = "red")
# g1 <- g1 + geom_point(data = splits[[1]][[2]], aes(x=X, y=Y), color = "blue")
# g1
# 
# 
# g2 <- ggplot()
# g2 <- g2 + geom_sf(data = map.input$ch, fill="#E8E8E8", color="black")
# g2 <- g2 + geom_point(data = splits[[2]][[1]], aes(x=X, y=Y), color = "red")
# g2 <- g2 + geom_point(data = splits[[2]][[2]], aes(x=X, y=Y), color = "blue")
# g2
# 
# g3 <- ggplot()
# g3 <- g3 + geom_sf(data = map.input$ch, fill="#E8E8E8", color="black")
# g3 <- g3 + geom_point(data = splits[[3]][[1]], aes(x=X, y=Y), color = "red")
# g3 <- g3 + geom_point(data = splits[[3]][[2]], aes(x=X, y=Y), color = "blue")
# g3
## ---- APPLY MODELS ----

 # :)

# Statistical models ####

# read in results produced by Jonas
file.name <- paste0(dir.models.output, "Output25112021.rds")

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T ){
    
    if(exists("outputs") == F){
        cat("File with statistical model outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
        stat.outputs <- readRDS(file = file.name)
        }
    else{
        cat("List with statistical model outputs already exists as object 'outputs' in this environment.")
    }
} else {
    
    cat("No statistical model outputs exist yet, we produce it and save it in", file.name)
    
    # #1) No cross validation
    # centered.occ.data <- center.splits(list(inv.occ), cv = F)
    # #a) No cross validation, no comm corr
    # res.1 <- stat_mod_cv(data.splits = centered.occ.data, cv = F, comm.corr = F)
    # 
    # #b) No cross validation, comm corr
    # res.1 <- stat_mod_cv(data.splits = centered.occ.data, cv = F, comm.corr = T)
    # 
    # #2) Cross validation
    # centered.splits <- lapply(splits, FUN = center.splits, cv = T)
    # 
    # #b) Cross validation, no comm corr
    #res.3 <- mclapply(centered.splits, mc.cores = 3, FUN = stat_mod_cv, cv = T, comm.corr = F)
    # 
    # #b) Cross validation, comm corr
    # res.4 <- mclapply(centered.splits, mc.cores = 3, FUN = stat_mod_cv, cv = T, comm.corr = T)
    
    # cat("Saving outputs of statistical model in", file.name)
    # saveRDS(stat.outputs, file = file.name)
    
}

# ML models ####

ptm <- proc.time() # to calculate time of simulation

# "Apply" null model
null.model <- apply.null.model(data = data, list.taxa = list.taxa, prev.inv = prev.inv)

file.name <- paste0(list.algo, collapse = "_")
file.name <- paste0(dir.models.output, file.name,"_", no.taxa, "taxa_", ifelse(CV, "CV_", "FIT_"), d, ".rds")

# file.name <- paste0("Q:/Abteilungsprojekte/siam/Jonas Wydler/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/Intermediate results/Trained models/", no.algo, "MLAlgoTrained.rds")

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T ){
    
    if(exists("outputs") == F){
        cat("File with ML outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
        outputs.CV <- readRDS(file = file.name)
        }
    else{
        cat("List with ML outputs already exists as object 'outputs' in this environment.")
    }

} else {
    
    if(CV == T){
        cat("No ML outputs exist yet, we produce it and save it in", file.name)
        
        # Compute one split after the other
        outputs <- lapply(centered.splits.factors, FUN = apply.ml.model, list.algo, list.taxa, env.fact)
        
        # Compute three splits in paralel (should be run on the server)
        # outputs <- mclapply(centered.splits.factors, mc.cores = 3, FUN = apply.ml.model, list.algo, list.taxa, env.fact)
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(outputs, file = file.name, version = 2)
    
        } else {
        
        # apply temporary on training data of Split1, later take whole dataset centered etc
        splitted.data <- list("Training data" =  centered.splits.factors$Split1$`Training data`, "Testing data" = data.frame())
        
        outputs.fit <- apply.ml.model(splitted.data = splitted.data, list.algo = list.algo, list.taxa = list.taxa,
                                      env.fact = env.fact, CV = F)
        outputs <- outputs.fit
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(outputs, file = file.name, version = 2)
        
        }
}

print(paste("Simulation time of different models ", info.file.name))
print(proc.time()-ptm)

# ANN ####

# read output from ann_model.r, hopefully ...

# Training GAM bam for 6 taxa: 5 hours

# source("ml_model_functions.r")
# source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# make mean over splits for final cross validation

if(CV == T){
    
    outputs.cv <- vector(mode = "list", length = length(list.algo))
    names(outputs.cv) <- list.algo
    
    for (l in 1:no.algo) {
        
        temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
        names(temp.list.st.dev) <- list.taxa
        
        for( j in 1:no.taxa){
            
            temp.vect <- vector(mode ="numeric", length = length(outputs)) 
            for (n in 1:length(outputs)) {
                temp.vect[n] <- outputs[[n]][[l]][[j]][["Performance testing set"]]
            }
            temp.list.st.dev[[j]] <- mean(temp.vect)
        }
    
        outputs.cv[[l]] <- temp.list.st.dev
    }
    # outputs.cv
    
    # plot it
    list.plots <- model.comparison.cv(outputs = outputs, outputs.cv = outputs.cv, null.model = null.model, list.algo = list.algo, list.taxa = list.taxa, prev.inv = prev.inv)
    
    file.name <- "ModelsComparCV.pdf"
    print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

    }
## ---- PLOTS ----

# Models comparison ####

# saved.outputs <- outputs
# outputs <- saved.outputs[[1]]
# outputs <- outputs2algo[[1]]

ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- model.comparison(outputs = outputs, null.model = null.model, list.algo = list.algo, list.taxa = list.taxa, prev.inv = prev.inv)

# Print the plots in a pdf file
file.name <- "ModelsComparFIT.pdf"
print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# Performance vs hyperparameters ####

ptm <- proc.time() # to calculate time of pdf production

# Compute plots
list.plots <- plot.perf.hyperparam(outputs = outputs, 
                                   list.algo = list.algo[2:length(list.algo)], # GLM algo doesn't have hyperparameters
                                   list.taxa = list.taxa)

# Print the plots in a pdf file
file.name <- "PerfvsHyperparam.pdf"
print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)


# Variables importance ####

ptm <- proc.time() # to calculate time of simulation

list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)

file.name <- "VarImp.pdf"
print.pdf.plots(list.plots = list.plots, width = 10, height = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

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
# list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
#                       list.taxa = list.taxa, env.fact = env.fact)
# 
# ptm <- proc.time() # to calculate time of simulation
# 
# # PDP of one model
# # list.plots <- plot.pdp(outputs = outputs, algo = "rf", list.algo = list.algo,
# #                       list.taxa = list.taxa, env.fact = env.fact)
# # 
# # file.name <- "PDP.pdf"
# # print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
# 
# # PDP of all models
# # We sub-select taxa and env.fact because it takes a lot of time
# list.plots <- plot.pdp(outputs = outputs, list.algo = list.algo,
#                        list.taxa = list.taxa, env.fact = env.fact)
# 
# file.name <- "allPDP.pdf"
# print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

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

print("hey")
# Stop execution if there is no (prediction on) testing set
# stopifnot(ratio != 1)
if( length(outputs[[1]][[1]]) < 7){
  break
}
print("salut")

# Map predictions ####

source("plot_functions.r")

ptm <- proc.time() # to calculate time of pdf production

inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

# make a list with all plots and plot them in a pdf
list.plots <- lapply(list.taxa, FUN = map.ml.pred.taxa, inputs, outputs, algo = list.algo[3], list.algo)

file.name <- "ObsvsPred_map.pdf"
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print(paste(file.name, "printing:"))
print(proc.time()-ptm)

# Response shape ####

ptm <- proc.time() # to calculate time of pdf production

# select an algorithm to plot
algo <- list.algo[1]

# make a list with all plots and plot them in a pdf
list.plots <- lapply(list.taxa, FUN = response.ml.pred.taxa, outputs, list.algo, env.fact, algo)


file.name <- paste0(algo, "_Resp_EnvFactvsTax.pdf")
print.pdf.plots(list.plots = list.plots, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

print("Producing PDF time:")
print(proc.time()-ptm)

