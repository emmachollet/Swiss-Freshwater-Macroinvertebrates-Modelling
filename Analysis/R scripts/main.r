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
if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") } # to arrange multiple plots on a page
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page
if ( !require("cowplot") ) { install.packages("cowplot"); library("cowplot") } # to arrange multiple plots on a page
if ( !require("pdp") ) { install.packages("pdp"); library("pdp") } # to plot partial dependance plots
if ( !require("gt") ) { install.packages("gt"); library("gt") } # to plot nice tables
if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") } # to plot nice tables
if ( !require("viridis")) {install.packages("viridis", repos="http://cloud.r-project.org"); library("viridis")} # to do even nicer plots
if ( !require("sf") ) { install.packages("sf"); library("sf") } # to read layers for map
if ( !require("scales") ) { install.packages("scales"); library("scales") } # to look at colors
if ( !require("reshape2") ) { install.packages("reshape2"); library("reshape2") } # to reshape dataframes

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
dl <- T # Data Leakage

# Set number of cores
n.cores.splits <-  1 # a core for each split, 3 in our case
n.cores.stat.models <- 1 # a core for each stat model 2 in our case (UF0, and CF0)
# Settings Stat models 
# Set iterations (sampsize), number of chains (n.chain), and correlation flag (comm.corr) for stan models,
# also make sure the cross-validation (CV) flag is set correctly
sampsize <- 2000 #10000 #I think this needs to be an even number for some reason (stan error)
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

# ECR: We need to discuss taxa selection (removing NAs reduce their number) ####
# JW: THERE IS ALSO ONE MORE SPECIES THAT GETS DROPPED IN THE CENTERING (ENDS UP AT 126 rather tha 127, move to
#data prep as well)

prepro.data <- preprocess.data(data.env = data.env, data.inv = data.inv, 
                                 env.fact.full = env.fact.full, dir.workspace = dir.workspace, 
                                 BDM = BDM, dl = dl, CV = CV)
data <- prepro.data$data
if(CV){ splits <- prepro.data$splits ; rem.taxa <- prepro.data$rem.taxa }
list.taxa <- prepro.data$list.taxa # list taxa after data pre-processing
list.taxa <- list.taxa[order(match(list.taxa, prev.inv$Occurrence.taxa))] # reorder taxa by prevalence
centered.data <- prepro.data$centered.data
centered.data.factors <- prepro.data$centered.data.factors
normalization.data <- prepro.data$normalization.data

# check normalization
# mean(centered.data$Split1$`Training data`$temperature)
# mean(centered.data$Split1$`Testing data`$temperature)
# 
# mean(centered.data$Split2$`Training data`$temperature)
# mean(centered.data$Split2$`Testing data`$temperature)
# 

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

# # Summary of prevalence of chosen taxa
# for ( i in list.taxa){
#     cat("Summary of absence, presence and NA for", list.taxa[i], ":", summary(data.inv[, list.taxa[i]]), "\n")
# }

# Select ml algorithms ####

# Select machine learning algorithms to apply (! their packages have to be installed first)
# Already select the colors assigned to each algorithms for the plots
list.algo <- c( "#256801" = 'glm', # Generalized Linear Model
                "#0A1C51" = 'gamSpline', # Deneralized Additive Model
                "#7B1359" = 'svmRadial', # Support Vector Machine
                "#AB4589" = 'rf') # Random Forest

no.algo <- length(list.algo)
                                          
## ---- APPLY MODELS ----

# Null model

null.model.full <- apply.null.model(data = data, list.taxa = list.taxa.full, prev.inv = prev.inv)
null.model <- null.model.full[list.taxa]

# Statistical models ####

ptm <- proc.time() # to calculate time of simulation
#file.name <- paste0(dir.models.output, "Stat_model_100iterations_corr_22taxa_CV_no_DL.rds") #to test

comm.corr.options <- c(T,F)
#comm.corr.options <- c(F)

names(comm.corr.options) <- c("CF0", "UF0")

stat.outputs <- mclapply(comm.corr.options, mc.cores = n.cores.stat.models, function(comm.corr){
  
  #comm.corr <- comm.corr.options[[1]]
  info.file.stat.name <- paste0("Stat_model_",
                                file.prefix,
                                no.taxa.full, "taxa_", 
                                sampsize,"iterations_",
                                ifelse(comm.corr,"CF0_","UF0_"),
                                ifelse(CV, "CV_", "FIT_"),
                                ifelse(dl, "DL_", "no_DL_"))
                                
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
  
        stat.output <- mclapply(centered.data, mc.cores = n.cores.splits, FUN = stat_mod_cv, CV, comm.corr, sampsize, n.chain)
        
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

# Process output from stat models to fit structure of ml models (makes plotting easier)
stat.outputs.transformed <- transfrom.stat.outputs(CV, stat.outputs)

# Make vector statistical models
list.stat.mod <- names(stat.outputs.transformed)
names(list.stat.mod) <- c("#59AB2D","#782B01")
#names(list.stat.mod) <- c("#59AB2D")

# Machine Learning models ####

ptm <- proc.time() # to calculate time of simulation

info.file.ml.name <-  paste0("ML_model_",
                             file.prefix, 
                             no.algo, "algo_",
                             no.taxa, "taxa_", 
                             ifelse(CV, "CV_", "FIT_"),
                             ifelse(dl, "DL_", "no_DL_"))

#JW: is this even still needed?
file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
#file.name.temp <- paste0(dir.models.output, "ML_model_All_2algo_2taxa_CV_no_DL_.rds")
cat(file.name)

# If the file with the outputs already exist we read it, else we run the algorithms
info1 <- paste0(file.prefix, no.algo, "algo_")
info2 <- paste0(ifelse(CV, "CV_", "FIT_"), ifelse(dl, "DL_", "no_DL_"))
list.ml.files <- list.files(path = dir.models.output)
file.name.temp <- paste0(dir.models.output, list.ml.files[which(grepl(info1, list.ml.files) & grepl(info2, list.ml.files))])
#file.name.temp <- paste0(dir.models.output, "ML_model_All_4algo_126taxa_CV_DL_.rds")
if( file.exists(file.name.temp) == T ){
    # & askYesNo(paste("Another file named", file.name.temp, "already exists. Should we read it and subselect taxa later ? If no, we run the algorithms for selected information."), default = T) ){

    cat("Reading", file.name.temp, "and saving it in object 'ml.outputs.fit'")
    if(CV){ ml.outputs.cv <- readRDS(file = file.name.temp)
      #ml.outputs.cv <- readRDS(file = paste0("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss-Freshwater-Macroinvertebrates-Modelling/Analysis/Intermediate results/Trained models/ML_model_All_4algo_126taxa_CV_DL_.rds"))
  
    } else { ml.outputs <- readRDS(file = file.name.temp) }

    } else {

        cat("No ML outputs exist yet, we produce it and save it in", file.name)

    if(CV == T){
        
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
# }

remove(info1, info2, list.ml.files, file.name.temp)

# Write down the number of splits and their names
if(CV){
    list.splits <- names(ml.outputs.cv)
    no.splits <- length(list.splits)
}

# ECR: Problem with some ml performance ####
if(CV){
    for (s in list.splits) {
    print(s)
        for (l in list.algo) {
            print(l)
            list.taxa.temp <- names(ml.outputs.cv[[s]][[l]])
            for (j in list.taxa.temp) {
                    perf <- ml.outputs.cv[[s]][[l]][[j]][["Performance testing set"]]
                    ml.outputs.cv[[s]][[l]][[j]][["Performance testing set"]] <- ifelse(perf > 1.5, Inf, perf)
                }
        }
    }
} else {
    for (l in list.algo) {
        print(l)
        list.taxa.temp <- names(ml.outputs[[l]])
        for (j in list.taxa.temp) {
            perf <- ml.outputs[[l]][[j]][["Performance training set"]]
            ml.outputs[[l]][[j]][["Performance training set"]] <- ifelse(perf > 1.5, Inf, perf)
        }
    }
}

print(paste("Simulation time of different models ", info.file.ml.name))
print(proc.time()-ptm)

if(!server){

# Neural Networks ####

# read output from ann_model.r, hopefully ...

# info1 <- paste0(file.prefix, no.ann, "ann_")
# info2 <- paste0(ifelse(CV, "CV_", "FIT_"), ifelse(dl, "DL_", "no_DL_"))
# list.ann.files <- list.files(path = dir.models.output)
# file.name.temp <- paste0(dir.models.output, list.ann.files[which(grepl(info1, list.ann.files) & grepl(info2, list.ann.files))])

if(CV){
  if(dl){
    file.name <- paste0(dir.models.output, "ANN_model_All_3ann_126taxa_CV_DL_.rds")
    ann.outputs.cv <- readRDS(file = file.name)
  }else{
    file.name <- paste0(dir.models.output, "ANN_model_All_3ann_126taxa_CV_no_DL_.rds")
    ann.outputs.cv <- readRDS(file = file.name)
  }
} else {
    file.name <- paste0(dir.models.output, "ANN_model_All_3ann_126taxa_FIT_no_DL_.rds")
    ann.outputs <- readRDS(file = file.name) 
    }

# Make vector neural networks
list.ann <- if(CV){names(ann.outputs.cv[[1]])} else {names(ann.outputs)}
names(list.ann) <- c("#A1491A", "#C46634", "#DF895A")

# Merge outputs ####

# Make final list of models
list.models <- c(list.algo, list.stat.mod, list.ann)
#list.models <- c(list.algo, list.ann)

no.models <- length(list.models)

# See the final colors for each model
show_col(names(list.models))

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

# Produce final outputs with mean performance across splits
if(CV){ 
    
    # Merge all CV outputs in one
    outputs.cv <- ml.outputs.cv
    for (s in list.splits) {
        #s = "Split2"
        outputs.cv[[s]][[list.stat.mod[1]]] <- stat.outputs.transformed[[1]][[s]]
        outputs.cv[[s]][[list.stat.mod[2]]] <- stat.outputs.transformed[[2]][[s]]
        outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv[[s]])
    }
    
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
    outputs <- append(append(ml.outputs, stat.outputs.transformed), ann.outputs)
    #outputs <- append(appendml.outputs, ann.outputs)
    
    # Make final outputs as tables
    df.perf <- make.df.outputs(outputs = outputs, list.models = list.models, 
                             list.taxa = list.taxa, list.splits = list.splits,
                             null.model = null.model, prev.inv = prev.inv, CV = CV)
}


# Look into the effect of data leakage ####

# save intermediate output to compare dl no dl 
#saveRDS(df.perf, file = paste0(dir.models.output, info.file.name, "df_perf_.rds"), version = 2)
#df.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_7models_126taxa_CV_no_DL_df_perf_.rds"))
#df.perf.dl <- readRDS(file = paste0(dir.models.output, "All_7models_126taxa_CV_DL_df_perf_.rds"))

# df.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_no_DL_df_perf_.rds"))
# df.perf.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_DL_df_perf_.rds"))
# 
# df.perf.no.dl$DL <- F
# df.perf.dl$DL <- T
# df.perf.dl.comb <- rbind(df.perf.no.dl,df.perf.dl)
# 
# #remove infite vlaues to take mean later
# df.perf.dl.comb$rf[!is.finite(df.perf.dl.comb$rf)] <- NA
# df.perf.dl.comb$gamSpline[!is.finite(df.perf.dl.comb$gamSpline)] <- NA
# df.perf.dl.comb$glm[!is.finite(df.perf.dl.comb$glm)] <- NA
# 
# colnames(df.perf.dl.comb)
# 
# df.perf.dl.comb.table <- as.data.frame(df.perf.dl.comb %>% group_by(DL) %>% 
#                                       summarise(glm.mean = mean(glm, na.rm = T),
#                                       glm.sd = sd(glm, na.rm = T),
#                                       
#                                       gamSpline.mean = mean(gamSpline,na.rm = T),
#                                       gamSpline.sd = sd(gamSpline,na.rm = T),
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
#                                       ANN3L32U.mean = mean(ANN3L32U, na.rm = T),
#                                       ANN3L32U.sd = sd(ANN3L32U, na.rm = T),
#                                       
#                                       ANN3L64U.mean = mean(ANN3L64U, na.rm = T),
#                                       ANN3L64U.sd = sd(ANN3L64U, na.rm = T),
#                                       
#                                       ANN5L32U.mean = mean(ANN5L32U, na.rm = T),
#                                       ANN5L32U.sd = sd(ANN5L32U, na.rm = T)
#                                       
#                                        ))
# 
# #saveRDS(df.perf.dl.comb.table, file = paste0(dir.plots.output, "Table_means_dl.rds"), version = 2)
# df.perf.dl.comb.table <- readRDS(file = paste0(dir.plots.output, "Table_means_dl.rds"))
# diff.means.dl <- df.perf.dl.comb.table[1,c("glm.mean", "gamSpline.mean","svmRadial.mean", 
#                          "rf.mean", "CF0.mean", "UF0.mean", "ANN3L32U.mean", "ANN3L64U.mean",
#                          "ANN5L32U.mean")] - 
#   df.perf.dl.comb.table[2,c("glm.mean", "gamSpline.mean","svmRadial.mean", 
#                             "rf.mean", "CF0.mean", "UF0.mean", "ANN3L32U.mean", "ANN3L64U.mean",
#                             "ANN5L32U.mean")]
# 

## ---- PLOTS ----
source("utilities.r")
source("ml_model_functions.r")
source("plot_functions.r")
# rm(list=ls())
# graphics.off()

# # Stat model traceplots ####
# 
# res <- stat.outputs[[1]][[1]][[1]]
# res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
# 
# print(traceplot(res,pars=c(names(res)[134:162],"lp__")))
# 
# print(traceplot(res))

# Models comparison ####

# Table with performance
if(CV){
        list.plots.cv <- plot.df.perf(df.perf = df.perf.cv, list.models = list.models, list.taxa = list.taxa, CV)
        list.plots <- plot.df.perf(df.perf = df.perf, list.models = list.models, list.taxa = list.taxa, CV)
        list.plots <- append(list.plots.cv, list.plots)
        # if(dl){
        # list.plots.dl <- plot.dl.perf(df.perf.dl.comb, list.models = list.models)
        # }
    } else {
        list.plots <- plot.df.perf(df.perf = df.perf, list.models = list.models, list.taxa = list.taxa, CV)
}

name <- "TablesPerf"
file.name <- paste0(name, ".pdf")
# list.plots.dl[[1]]


# Performance against prevalence and boxplots
list.plots <- model.comparison(df.perf = df.perf, list.models = list.models, CV = CV)
name <- "ModelsCompar"
file.name <- paste0(name, ".pdf")
print.pdf.plots(list.plots = list.plots, width = 12, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)

# # Print a jpeg file
# file.name <- paste0(name, ".jpg")
# jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
# print(list.plots[[1]])
# dev.off()

# Plots specifically related to trained models (and not to CV)

if(!CV){
    
    # Performance vs hyperparameters ####
    
    # Compute plots
    list.plots <- plot.perf.hyperparam(outputs = outputs, 
                                       list.algo = list.algo, # GLM algo doesn't have hyperparameters
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
    
    list.list.plots <- lapply(list.taxa, FUN= plot.ice.per.taxa, outputs, list.algo, env.fact, normalization.data)
    
    for (j in list.taxa) {
        taxon <- sub("Occurrence.", "", j)
        file.name <- paste0("ICE_", taxon, ".pdf")
        print.pdf.plots(list.plots = list.list.plots[[j]], width = 20, height = 20, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
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
