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

# Set if we want to fit models to whole dataset or perform cross-validation (CV)
CV <- T # Cross-Validation
dl <- F # Data Leakage

# Set number of cores
n.cores <-  3

# Settings Stat models 
# Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models,
# also make sure the cross-validation (CV) flag is set correctly
sampsize <- 1000 #10000 #I think this needs to be an even number for some reason (stan error)
n.chain  <- 2 #2
comm.corr <- F

# Select taxa
all.taxa = T
# set to FALSE if it's for exploration (very few taxa or only with intermediate prevalence)
# set to TRUE to apply models to all taxa


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

# MOVE THIS TO DATA PREPARATION SCRIPTS ? ####
# JW: THERE IS ALSO ONE MORE SPECIES THAT GETS DROPPED IN THE CENTERING (ENDS UP AT 126 rather tha 127, move to
#data prep as well)

prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, 
                                     env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                                     dl = dl, CV = CV)
data <- prepro.data$data
if(CV){ splits <- prepro.data$splits }
centered.data <- prepro.data$centered.data
centered.data.factors <- prepro.data$centered.data.factors
normalization.data <- prepro.data$normalization.data

# Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c("#030AE8" = 'glm', # Random Forest
               # "#048504" = 'bam', # Generalized Additive Model using splines
               "#948B8B" = 'gamSpline',#
               # 'earth', # MARS: Multivariate Adaptive Regression Splines
               # "#A84E05" = 'elm', # Extreme Learning Machine (Neural Network)
               # 'bayesglm') #, # Bayesian Generalized Linear Model
               "#DB1111" = 'svmRadial', # Support Vector Machine
               # "#DB1111" = 'svmPoly', # Support Vector Machine
               # "#DB1111" = 'svmLinear', # Support Vector Machine
               "#790FBF" = 'rf') # Random Forest

no.algo <- length(list.algo)

# Select taxa ####

cind.taxa <- which(grepl("Occurrence.",colnames(data)))
list.taxa.full <- colnames(data)[cind.taxa]
list.taxa <- list.taxa.full

# Select taxa for prediction
if (all.taxa == F){
    # 2 taxa
    #list.taxa       <- list.taxa.full[list.taxa.full %in% c("Occurrence.Gammaridae", "Occurrence.Heptageniidae")]
    # 
    # 5 taxa
    # list.taxa       <- list.taxa.full[list.taxa.full %in% prev.inv[which(prev.inv[, "Prevalence"] < 0.7 & prev.inv[,"Prevalence"] > 0.55),
    #                            "Occurrence.taxa"]] # Select only few taxa
    
    # 22 taxa
    list.taxa       <- list.taxa.full[list.taxa.full %in% prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25),
                                  "Occurrence.taxa"]] # Select with prevalence percentage between 25 and 75%

}

no.taxa.full <- length(list.taxa.full)
no.taxa <- length(list.taxa)

# # Summary of prevalence of chosen taxa
# for ( i in 1:no.taxa){
#     cat("Summary of absence, presence and NA for", list.taxa[i], ":", summary(data.inv[, list.taxa[i]]), "\n")
# }

# Write information for file names
# percentage.train.set <- ratio * 100
info.file.name <- paste0(file.prefix, 
                         # d, # don't need to include the date
                         no.taxa, "taxa_", 
                         # no.env.fact, "envfact_",
                         no.algo, "algo_",
                         ifelse(CV, "CV_", "FIT_"),
                         ifelse(dl, "DL_", "no_DL_"),
                         # "trainset", percentage.train.set, 
                         # if( ratio != 1) {split.var}, 
                         "")

info.file.stat.name <- paste0("Stat_model_",
                              # d, # don't need to include the date
                              no.taxa.full, "taxa_", 
                              # no.env.fact, "envfact_",
                              sampsize,"iterations_",
                              ifelse(comm.corr,"CF0_","FF0_"),
                              ifelse(CV, "CV_", "FIT_"),
                              ifelse(dl, "DL_", "no_DL_"),
                              # "trainset", percentage.train.set, 
                              # if( ratio != 1) {split.var}, 
                              "")

### ----Investigate splits----

# mean.train1 <- lapply(splits, function(split){
#   split <- splits[[1]]
#   test <- split[[1]][,env.fact]
#   return(mean(split[[1]][,env.fact]))
# })
# 
# 
# mean.train1 <- as.numeric(lapply(centered.splits[[1]][[1]][,env.fact], FUN = mean))
# mean.train2 <- as.numeric(lapply(splits[[2]][[1]][,env.fact], FUN = mean))
# mean.train3 <- as.numeric(lapply(splits[[3]][[1]][,env.fact], FUN = mean))
# 
# dif1 <- 1-(mean.train1/mean.train2)
# dif2 <- 1-(mean.train2/mean.train3)
# 
# splits2 <- split.data(data, 1)
# 
# mean2.train1 <- as.numeric(lapply(splits2[[1]][[1]][,env.fact], FUN = mean))
# mean2.train2 <- as.numeric(lapply(splits2[[2]][[1]][,env.fact], FUN = mean))
# mean2.train3 <- as.numeric(lapply(splits2[[3]][[1]][,env.fact], FUN = mean))
# 
# dif3 <- 1-(mean2.train1/mean2.train2)
# dif4 <- 1-(mean2.train2/mean2.train3)
# 
# 
# splits3 <- split.data(data, 1)
# mean3.train1 <- as.numeric(lapply(splits3[[1]][[1]][,env.fact], FUN = mean))
# mean3.train2 <- as.numeric(lapply(splits3[[2]][[1]][,env.fact], FUN = mean))
# mean3.train3 <- as.numeric(lapply(splits3[[3]][[1]][,env.fact], FUN = mean))
# dif5 <- 1-(mean3.train1/mean3.train2)
# dif6 <- 1-(mean3.train2/mean3.train3)
# 
# testest <- as.data.frame(rbind(dif1,dif2,dif3,dif4, dif5, dif6))
# testest <- as.data.frame(rbind(dif1,dif2))
# 
# colnames(testest) <- env.fact
# testest.table <- lapply(testest, mean, 1)
# 
# 
# split1 <- splits[[1]]
# split2 <- splits[[2]]
# split3 <- splits[[3]]
# 
# map.inputs <- map.inputs(dir.env.data, data.env)
# 
# 
# g1 <- ggplot()
# g1 <- g1 + geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black")
# g1 <- g1 + geom_point(data = splits[[1]], aes(x=X, y=Y), color = "red")
# g1 <- g1 + geom_point(data = splits[[2]], aes(x=X, y=Y), color = "blue")
# g1
# 
# g2 <- ggplot()
# g2 <- g2 + geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black")
# g2 <- g2 + geom_point(data = split2[[1]], aes(x=X, y=Y), color = "red")
# g2 <- g2 + geom_point(data = split2[[2]], aes(x=X, y=Y), color = "blue")
# g2
# 
# g3 <- ggplot()
# g3 <- g3 + geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black")
# g3 <- g3 + geom_point(data = split3[[1]], aes(x=X, y=Y), color = "red")
# g3 <- g3 + geom_point(data = split3[[2]], aes(x=X, y=Y), color = "blue")
# g3

# explore NAs
# explore(data) # to see

# vis_miss(data[,cind.taxa], cluster = FALSE, warn_large_data = FALSE)
# 
# too.many.na <- c()
# for(i in cind.taxa){
#     if(sum(is.na(data[,i])) > 200){ too.many.na <- c(too.many.na, i)}
# }
# colnames(data)[too.many.na]
# # data[complete.cases(data), "SampleId"]
# 
# data.test1 <- data[, -too.many.na]
# data.test2 <- na.omit(data.test1)
# 
# setdiff(data$SampId, data.test2$SampId)
# # setdiff(data.test1$SampId, data.test2$SampId)

# in ALL, loosing 5 taxa and 50 rows

# Explore geographical distribution

# inputs <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)
# 
# variable <- "BIOGEO"
# temp.data.env <- data.env[, c("X","Y", variable)]
# no.rows <- nrow(temp.data.env)
# no.na <- sum(is.na(temp.data.env[,variable]))
# explanation <- "Geographical distribution"
# 
# g <- ggplot()
# g <- g + geom_sf(data = inputs$ch, fill="#E8E8E8", color="black")
# g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
# g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
# g <- g + geom_point(data = temp.data.env, aes(x=X, y=Y, color= temp.data.env[, variable]), size= 3, alpha= 0.6)
# 
# g <- g + theme_void(base_size = 18)
# g <- g + theme(# plot.title = element_text(hjust = 0.5),
#     panel.grid.major = element_line(colour="transparent"),
#     plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
# g <- g + labs(title = paste("Geographic distribution of",variable),
#               subtitle = paste0(no.na, " NAs out of ", no.rows, " samples \n", explanation), colour = variable)
# 
# print(g)

## ---- APPLY MODELS ----

 # :)

# Statistical models ####

ptm <- proc.time() # to calculate time of simulation
#file.name <- paste0(dir.models.output, "Stat_model_100iterations_corr_22taxa_CV_no_DL.rds") #to test
file.name <- paste0(dir.models.output, info.file.stat.name, ".rds")

cat(file.name)

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T){
    
    if(exists("stat.outputs") == F){
        cat("File with statistical model outputs already exists, we read it from", file.name, "and save it in object 'stat.outputs'")
        stat.outputs <- readRDS(file = file.name)
        }
    else{
        cat("List with statistical model outputs already exists as object 'stat.outputs' in this environment.")
    }
} else {
    
    cat("No statistical model outputs exist yet, we produce it and save it in", file.name)
    
    if(CV == T){

      stat.outputs <- mclapply(centered.splits, mc.cores = n.cores, FUN = stat_mod_cv, CV, comm.corr, sampsize, n.chain)
      
      cat("Saving outputs of statistical models in", file.name)
      saveRDS(stat.outputs, file = file.name, version = 2) #version two here is to ensure compatibility across R versions
      
    } else {
      # apply temporary on training data of Split1, later take whole dataset centered etc

      stat.outputs <- stat_mod_cv(data.splits = centered.data, CV, comm.corr, sampsize, n.chain = n.chain)
      cat("Saving outputs of statistical models in", file.name)
      saveRDS(stat.outputs, file = file.name, version = 2)
      names(stat.outputs) <- c("stan.object", "output")
      }
    
}

print(paste("Simulation time of statistical model ", info.file.stat.name))
print(proc.time()-ptm)


## ---- Process output from stat models
# Model performance (FIT)
# FF0
CV <-  F
comm.corr <- F


fit.output <- stat.outputs$output
dev.tmp <- stat.outputs$output$deviance


# #plot traceplots
# # res <- stat.outputs[[1]][[1]]
# # res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
# # 
# # print(traceplot(res,pars=c(names(res)[1:32],"lp__")))
# # 
# # print(traceplot(res))
# 
# #exract neeeded output (std.deviance from the testing set in this case) from the output of the stat models.
# stat_cv_nocorr <- stat.outputs
# stat_cv_nocorr_res <- lapply(stat_cv_nocorr, function(split){
#   #split <- stat_cv_nocorr[[1]]
#   dev_temp <- split[[2]]$deviance
#   dev_temp$std.deviance <- dev_temp$std.deviance
#   dev_temp <- dev_temp[c("Taxon", "Type", "std.deviance")]
#   tmp_dev_test <- subset(dev_temp, Type = "Testing")
#   return(tmp_dev_test)
# })
# 
# print(traceplot(res))
stat_model_name <- "CF0"
#exract neeeded output (std.deviance from the testing set in this case) from the output of the stat models.
stat_cv_nocorr <- stat.outputs
stat_cv_nocorr_res <- lapply(stat_cv_nocorr, function(split){
  #split <- stat_cv_nocorr[[1]]
  dev_temp <- split[[2]]$deviance

  dev_temp <- split[[2]]$deviance
  dev_temp$std.deviance <- dev_temp$std.deviance
  dev_temp <- dev_temp[c("Taxon", "Type", "std.deviance")]
  tmp_dev_test <- subset(dev_temp, Type = "Testing")
  return(tmp_dev_test)
})

#bind rows for all three splits
stat_cv_nocorr_res_table <- stat_cv_nocorr_res[[1]]
for(i in 2:length(stat_cv_nocorr_res)){
  stat_cv_nocorr_res_table <- rbind(stat_cv_nocorr_res_table, stat_cv_nocorr_res[[i]])

}

#calculate mean std.deviance across splits
stat_cv_nocorr_res_table <- as.list(stat_cv_nocorr_res_table %>% group_by(Taxon) %>% summarise(performance = mean(std.deviance, na.rm = T)))

stat.output.list <- vector(mode = "list", length = length(stat.outputs))

for(n in 1:length(stat.outputs)){
  #n = 1
  temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))

  for(j in 1:length(list.taxa)){
    #i = 1
    temp.dat.dev <- subset(stat.outputs[[n]][[2]]$deviance, Type == "Testing")
    
    temp.dat.dev <- subset(temp.dat.dev, Taxon == list.taxa[[j]])
    
    temp.dat.dev$Performance <- as.numeric(temp.dat.dev$std.deviance)
    temp.list.st.dev[[j]] <- temp.dat.dev
    #names(temp.list.st.dev[[2]]) <- "Performance testing set"
    #temp.dat.prop <- subset(stat.outputs$Split1[[2]]$probability, Type == "Testing")
    #temp.dat.prop$Likelyhood <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
    
    #temp.list.st.dev[[1]] <- list("Performance" = temp.dat.dev$Performance)
    
  }
  names(temp.list.st.dev) <-  list.taxa
  stat.output.list[[n]] <- temp.list.st.dev

}
names(stat.output.list) <- c("Split 1", "Split 2", "Split 3")
#start new workflow to combine models
if(CV == T){
  names(stat_model_name) = "#DD1C77"
  list.algo.comb <- stat_model_name
  
  outputs <- vector(mode = "list", length = length(list.algo.comb))
  names(outputs) <- list.algo.comb
  
  no.algo <- length(list.algo.comb)
  for (l in 1:no.algo) {
    temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
    names(temp.list.st.dev) <- list.taxa
    
    for( j in 1:no.taxa){
      temp.vect <- vector(mode ="numeric", length = length(outputs))
      
      #for (n in 1:length(stat.outputs)) {
        #n = 1
        temp.dat <- subset(stat.outputs$Split1[[2]]$deviance, Type == "Testing")
        temp.dat <- subset(temp.dat, Taxon == list.taxa[[j]])
        temp.vect[1] <- temp.dat$std.deviance
        names(temp.vect[1]) <- "Performance testing set"
        #}
        temp.list.st.dev[[j]] <- mean(temp.vect)
        
    }
  }
}



stat.outputs.perf <- stat.outputs$Split1[[2]]$deviance$std.deviance

temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
temp.list.st.dev[[1]]

# Machine Learning models ####

ptm <- proc.time() # to calculate time of simulation

# "Apply" null model
null.model <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)

file.name <- paste0(dir.models.output, info.file.name,"MLAlgoTrained", ".rds")
# file.name <- paste0(dir.models.output, "glm_gamSpline_svmRadial_rf_22taxa_FIT.rds")
cat(file.name)
# file.name <- paste0("Q:/Abteilungsprojekte/siam/Jonas Wydler/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/Intermediate results/Trained models/", no.algo, "MLAlgoTrained.rds")

# If the file with the outputs already exist, just read it
if (file.exists(file.name) == T ){
    
    if(exists("outputs") == F){
        cat("File with ML outputs already exists, we read it from", file.name, "and save it in object 'outputs'")
        outputs <- readRDS(file = file.name)
        }
    else{
        cat("List with ML outputs already exists as object 'outputs' in this environment.")
    }

} else {
    
    if(CV == T){
        cat("No ML outputs exist yet, we produce it and save it in", file.name)
        
        # Compute one split after the other
        outputs.cv <- lapply(centered.splits.factors, FUN = apply.ml.model, list.algo, list.taxa, env.fact, CV)
        
        # Compute three splits in paralel (should be run on the server)
        # outputs <- mclapply(centered.splits.factors, mc.cores = 3, FUN = apply.ml.model, list.algo, list.taxa, env.fact)
        
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(outputs.cv, file = file.name, version = 2)
    
        } else {
        
        # apply temporary on training data of Split1, later take whole dataset centered etc
        splitted.data <- list("Training data" =  centered.data.factors[[1]], "Testing data" = data.frame())
        
        outputs.fit <- apply.ml.model(splitted.data = splitted.data, list.algo = list.algo, list.taxa = list.taxa,
                                      env.fact = env.fact, CV = F, prev.inv = prev.inv)
        outputs <- outputs.fit
        cat("Saving outputs of algorithms in", file.name)
        saveRDS(outputs, file = file.name, version = 2)
        
        }
}

print(paste("Simulation time of different models ", info.file.name))
print(proc.time()-ptm)


# Neural Networks ####

# read output from ann_model.r, hopefully ...

# file.name <- paste0(dir.models.output, "test_ANNoutputs.rds")
# ann.outputs <- readRDS(file = file.name)


# Training GAM bam for 6 taxa: 5 hours

# source("ml_model_functions.r")
# source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# make mean over splits for final cross validation

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
    # temp.list.stat <- vector(mode = "list", length = length(stat_cv_nocorr_res_table$Taxon))
    # names(temp.list.stat) <- stat_cv_nocorr_res_table$Taxon
    # 
    # for( j in 1:length(stat_cv_nocorr_res_table$Taxon)){
    #   
    #   temp.vect <- vector(mode ="numeric", length = length(stat_cv_nocorr_res_table$Taxon))
    #   for (n in 1:length(stat_cv_nocorr_res_table$Taxon)) {
    #     temp.list.stat[n] <- stat_cv_nocorr_res_table$performance[[n]]
    #   }
    # }
    # FF0 <- as.list(temp.list.stat)
    # FF0 <- FF0[list.taxa]
    # 
    # outputs[[5]] <- FF0
    # names(outputs)[[5]] <- "FF0"
    
    # # Add output of ANN
    # for (a in 1:length(ann.outputs)) {
    #     outputs[[4 + a]] <- ann.outputs[[a]]
    #     names(outputs)[4 + a] <- names(ann.outputs)[a]
    #     
    # }
}

## ---- PLOTS ----

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
# list.algo <- c(list.algo, "green" = "FF0")
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

