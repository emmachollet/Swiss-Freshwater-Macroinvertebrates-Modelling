## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Explore processed data ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 21, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Import libraries and data ----

# Load packages
if ( !require(dplyr) ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require(ggplot2) ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
if ( !require(skimr) ) { install.packages("skimr"); library("skimr") } # to show key descriptive stats
if ( !require(sf) ) { install.packages("sf"); library("sf") } # to read GIS data (shape files)
# if ( !require(caret) ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models
if ( !require(ggpubr) ) { install.packages("ggpubr"); library("ggpubr") } # to make nice arrangement of nice plots
if ( !require(corrplot) ) { install.packages("corrplot"); library("corrplot") } # to make nice arrangement of nice plots


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

file.rank.env     <- "ranking_env_data.csv"
file.env.explan   <- "environmentaldata_documentation_20210728.csv"

# Read data

# Only colnames of env fact sorted by importance
rank.env          <- read.csv(paste(dir.env.data,file.rank.env,sep=""),header=TRUE, sep=";", stringsAsFactors=FALSE)
env.explan        <- read.csv(paste(dir.env.data,file.env.explan,sep=""),header=TRUE, sep=";", stringsAsFactors=FALSE)

# Read inv and env data, set if we want to compute for the All or the BDM dataset
BDM <- T

if( BDM == TRUE){
    
    data.inv      <- read.delim(paste(dir.inv.data,file.inv.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
#    prev.inv      <- read.delim(paste(dir.inv.data,file.prev.BDM,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
    data.env      <- read.delim(paste(dir.env.data,file.env.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
    file.prefix   <- "BDM_"
    
} else {
    
    data.inv          <- read.delim(paste(dir.inv.data,file.inv.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
#    prev.inv          <- read.delim(paste(dir.inv.data,file.prev,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
    data.env          <- read.delim(paste(dir.env.data,file.env.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
    file.prefix       <- "All_"
    
}

## ----  Load functions ----

source("model_functions.r")
source("plot_functions.r")

## ---- Preliminaries to explore environmental data ----

# Convert character columns to factors
data.env[sapply(data.env, is.character)] <- lapply(data.env[sapply(data.env, is.character)], 
                                                   as.factor)
# Make a table with nb of missing values, mean, sd 
# skimmed <- skim(data.env)

# Use the file ranking.data to make a list of useful factors to explore

# replace "_" and " " by "." in colnames to be consistent
colnames(rank.env) <- gsub("_", ".", colnames(rank.env))
colnames(rank.env) <- gsub(" ", ".", colnames(rank.env))
env.explan$column.name <- gsub("_", ".", env.explan$column.name)
env.explan$column.name <- gsub(" ", ".", env.explan$column.name)

# sort the columns in 4 categories: the sample/site information = 3,
# the environmental factors to, keep in priority = 2, keep to explore = 1, exclude = 0
info    <- colnames(rank.env)[which(rank.env[1,] == 3)]
prio    <- colnames(rank.env)[which(rank.env[1,] == 2)]
explo   <- colnames(rank.env)[which(rank.env[1,] == 1)]
excl    <- colnames(rank.env)[which(rank.env[1,] == 0)]

# make a list of env. fact. we want to explore
a.lot = T    

if (a.lot == T){
    
    env.fact <- prio
    
} else if (a.lot == F){
    
    env.fact <- c("temperature",
                  "velocity",
                  "cow.density",
                  "IAR",
                  "urban.area",
                  "FRI",
                    "WALD.ANT",
                    "width.variability",
                    "bed.modification",
                    "morphology",
                    "A.EDO",
                    "F.EDO",
                    "ARA.fraction",
                    "agri.land",
                    "Slope",
                    "fields",
                    "saprobic.cond",
                    "normcov.mobile.blocks",
                    "normcov.coarse.inorganic.sediments",
                    "normcov.gravel",
                    "normcov.sand.silt",
                    "normcov.fine.sediments")
}


if(BDM != TRUE) {
    # remove InS env. fact.
    cind <- c(grep("InS.",env.fact),
              grep("covclass.",env.fact),
              grep("normcov.",env.fact),
              grep("sfract.",env.fact),
              grep("v.",env.fact, fixed = TRUE)) # has to match exactly, to avoid to remove .variability or .vegetation
    env.fact <- env.fact[-cind]
}

# maybe missing factors from data.env
env.fact <- env.fact[which(env.fact %in% colnames(data.env))]

no.env.fact <- length(env.fact)

## ---- Plot correlation matrix ----

ptm <- proc.time() # to calculate time of pdf production

file.name <- paste0(no.env.fact,"envfact_", "CorrMatrix.pdf")

pdf(paste0(dir.output, file.prefix, file.name), paper = 'special', 
    width = 23,
    height = 23,
    onefile = TRUE)

plot.data <- na.omit(data.env[,env.fact])
plot.data <- as.matrix(plot.data)
plot.data <- cor(plot.data)

corrplot(plot.data, 
         method = "number", 
         type="upper", 
         tl.col = "black",
         tl.srt=45)
corrplot.mixed(plot.data, 
               tl.col = "black",
               tl.srt=45,
               tl.pos = "lt",
               order="hclust", 
               lower = 'number', upper = "circle", 
               lower.col = "black", 
               number.cex = .7)


dev.off()    
print("Producing PDF time:")
print(proc.time()-ptm)


## ---- Plot environmental factors on swiss map ----

ptm <- proc.time() # to calculate time of pdf production

file.name <- "EnvFact_OnMap.pdf"

map.inputs.d <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)

map.env.fact(inputs = map.inputs.d, env.fact = env.fact, data.env = data.env, env.explan = env.explan,
             dir.output = dir.output, file.prefix = file.prefix, file.name = file.name)

# 
# file.name <- "EnvFact_OnMap_Morph.pdf"
# map.env.fact(inputs = map.inputs.d, env.fact = env.fact[13:25], data.env = data.env,
#              dir.output = dir.output, file.prefix = file.prefix, file.name = file.name)

# # to do pdf outside of the function
# pdf(paste0(dir.output, file.prefix, file.name), paper = 'special', 
#     width = 11, 
#     onefile = TRUE)
# map.env.fact(inputs = map.inputs.d, env.fact = env.fact, data.env = data.env,
#              dir.output = dir.output, file.prefix = file.prefix, file.name = NA)

# aukasoù
dev.off()

print("Producing PDF time:")
print(proc.time()-ptm)
# Producing All_pdf -> 4 min
# Producing BDM_pdf -> 4 min

data.env[which(is.na(data.env$temperature)), "Canton"]


# plot.title <- "whole environmental dataset"
# 
# pdf(paste0(dir.output,fileName), paper = 'special', width = 11, onefile = TRUE)
# 
# for (k in 1:length(env.fact)){
#     variable <- env.fact[k]
#     
#     # Set up scales
#     k.min <- round(min(data.env[,variable], na.rm = T), 1)
#     k.max <- round(max(data.env[,variable], na.rm = T), 1)
#     k.int <- (k.max - k.min)/5 # ; if.int <- round(if.int)
#     
#     g <- ggplot()
#     g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
#     g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
#     g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
#     
#     g <- g + geom_point(data = data.env, aes(x=X, y=Y, size=data.env[, env.fact[k]]), alpha = 0.35, color="black")
#     g <- g + scale_radius(name = variable, limits = c(k.min, k.max), breaks = seq(k.min, k.max, k.int), range = c(2, 7))
#     
#     g <- g + theme_void(base_size = 15)
#     g <- g + theme(plot.title = element_text(hjust = 0.5),
#                    panel.grid.major = element_line(colour="transparent"),
#                    plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
#     
#     # g <- g + scale_colour_brewer(palette = "Set1")
#     # g <- g + guides(colour = guide_legend(override.aes = list(size=6)), 
#     # ?? don't know for what this is, and it doen't work for me #####
#                     #size = guide_legend(title=labeller(variable)))
#     g <- g + labs(title = paste("Explanatory variables in ",plot.title,": ",variable, sep=""))
#     print(g)
# }
# 
# dev.off()
# 


# ## Same for BDM
# 
# # Construct BDM data set
# data.env.BDM <- data.env %>%
#     filter(MonitoringProgram == "BDM") # filter the rows, for columns use select
# # data[which(data.env$MonitoringProgram =="BDM") , ]
# dim(data.env.BDM)
# 
# # BDM SampId and row indices for additional checks
# BDM.samples <- data.env[which(data.env$MonitoringProgram =="BDM") , "SampId"]
# length(BDM.samples) # should have 886 samples
# BDMind <- which(data.env$SampId %in% BDM.samples)
# 
# env.fact <- c("IAR", "saprobic.cond", "urban.area", "velocity", "temperature")
# 
# fileName <- "BDM_EnvFact_OnMap.pdf"
# plot.title <- "BDM sites only environmental dataset"
# 
# pdf(paste0(dir.output,fileName), paper = 'special', width = 11, onefile = TRUE)
# 
# for (k in 1:length(env.fact)){
#     variable <- env.fact[k]
#     
#     # Set up scales
#     k.min <- round(min(data.env.BDM[,variable], na.rm = T), 1)
#     k.max <- round(max(data.env.BDM[,variable], na.rm = T), 1)
#     k.int <- (k.max - k.min)/5 # ; if.int <- round(if.int)
#     
#     g <- ggplot()
#     g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
#     g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
#     g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
#     
#     g <- g + geom_point(data = data.env.BDM, aes(x=X, y=Y, size=data.env.BDM[, env.fact[k]]), alpha = 0.35, color="black")
#     g <- g + scale_radius(name = variable, limits = c(k.min, k.max), breaks = seq(k.min, k.max, k.int), range = c(2, 7))
#     
#     g <- g + theme_void(base_size = 15)
#     g <- g + theme(plot.title = element_text(hjust = 0.5),
#                    panel.grid.major = element_line(colour="transparent"),
#                    plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
#     
#     # g <- g + scale_colour_brewer(palette = "Set1")
#     # g <- g + guides(colour = guide_legend(override.aes = list(size=6)), 
#     # ?? don't know for what this is, and it doen't work for me #####
#     #size = guide_legend(title=labeller(variable)))
#     g <- g + labs(title = paste("Explanatory env.fact in ",plot.title,": ",variable, sep=""))
#     print(g)
# }
# 
# dev.off()
# 
# 
# ## ---- Plot probability of occurrence ----
# 
# # Change occ.taxa in factors absent-present (for ALL)
# 
# cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
# 
# for (i in cind.taxa ) {
#     data.inv[which(data.inv[,i] == 0),i] <- "absent"
#     data.inv[which(data.inv[,i] == 1),i] <- "present"
#     data.inv[,i] = as.factor(data.inv[,i])
# }
# 
# # cind.taxa <- which(grepl("Occurrence.",colnames(data.inv.BDM)))
# # 
# # for (i in cind.taxa ) {
# #     data.inv.BDM[which(data.inv.BDM[,i] == 0),i] <- "absent"
# #     data.inv.BDM[which(data.inv.BDM[,i] == 1),i] <- "present"
# #     data.inv.BDM[,i] = as.factor(data.inv.BDM[,i])
# # }
# 
# 
# # Choose taxa (in ALL) with intermediate prevalence
# taxa.int.prev <- data.prev[which(data.prev[, "Percentage.prevalence"] < 75 & data.prev[,"Percentage.prevalence"] > 25), "Occurrence.taxa"]
# occ.taxa <- c("Occurrence.Hydraenidae")
# 
# # Construct main dataset (with env.fact and taxa)
# data <- data.env %>%
#     left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
# dim(data)
# 
# # Remove NAs
# taxa.env <- na.omit(data[, c(occ.taxa, env.fact)])
# 
# # Prepare data
# # Step 1: Get row numbers for the training data
# trainRowNumbers <- createDataPartition(taxa.env[,occ.taxa], p=0.8, list=FALSE) # split into training 80% and test 20% datasets
# 
# # Step 2: Create the training  dataset
# trainData <- taxa.env[trainRowNumbers,]
# 
# # Step 3: Create the test dataset
# testData <- taxa.env[-trainRowNumbers,]
# 
# # Store X and Y for later use.
# x = trainData[, env.fact]
# y = trainData[, occ.taxa]
# 
# f <- reformulate(".", occ.taxa)
# 
# # Train the model using rf and predict on the training data itself
# model_rf = train(f, data=trainData, method='rf')
# fitted <- predict(model_rf)
# 
# # Predict on testData
# predicted <- predict(model_rf, testData, type = "prob")
# predicted.fact <- predict(model_rf, testData)
# 
# # Compute the confusion matrix
# confusionMatrix(reference = testData[,occ.taxa], data = predicted.fact, mode='everything', positive='present')
# 
# 
# # Construct dataset to plot
# plot.data <- data[as.numeric(rownames(testData)), c("SiteId", "SampId", "X", "Y", occ.taxa, env.fact)]
# plot.data$pred <- predicted[, "present"]
# 
# # Plot proba of occ on map
# 
# # fileName <- "ALL_ProbPred_OnMap.pdf"
# # pdf(paste0(dir.output,fileName), paper = 'special', width = 11, onefile = TRUE)
# # 
# # # Map geometries
# # g <- ggplot()
# # g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
# # g <- g + geom_point(data = plot.data, aes(X, Y, size = pred, color = testData$Occurrence.Hydraenidae)) #alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
# # # !! defined with precise taxa, not generalized ####
# # 
# # # Configure themes and labels
# # g <- g + labs(title = paste("Probability of occurrence vs observations of", occ.taxa),
# #               subtitle = paste("Random Forest:", paste(env.fact, collapse = " ", sep = " ")), #, "- page", j),
# #               x = "",
# #               y = "",
# #               size = "Probability of\noccurrence",
# #               # alpha = "Posterior",
# #               color = "Observation")
# # 
# # # Configure legends and scales
# # g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
# #                 # alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
# #                 color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
# # # g <- g + scale_y_continuous(breaks=NULL)
# # # g <- g + scale_x_continuous(breaks=NULL)
# # # g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 6))
# # g <- g + scale_color_manual(values=c(absent = "#FF0000", present = "#0077FF"), labels=c("Absence", "Presence"))
# # # g <- g + scale_alpha_manual(values=c("0.65"="0.65", "0.35"="0.35"), labels=c("5th quantile", "95th quantile"))
# # g <- g + scale_shape_identity() # Plot the shape according to the data
# # 
# # print(g)
# # 
# # dev.off()
# 
# 
# ## ---- Speficif plot proba on map for poster ----
# 
# fileName <- "Plot_OccHydraenidae_SEFS12Poster.png"
# # pdf(paste0(dir.output,fileName), paper = 'special', width = 11, onefile = TRUE)
# 
# taxon <- sub("Occurrence.", "", occ.taxa)
# 
# # Map geometries
# g <- ggplot()
# g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
# g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
# g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
# g <- g + geom_point(data = plot.data, aes(X, Y, size = pred, color = testData$Occurrence.Hydraenidae), alpha =0.7) #alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
# # !! defined with precise taxa, not generalized ####
# 
# # Configure themes and labels
# g <- g + theme_void()
# g <- g + theme(plot.title = element_text(size = 16),#, hjust = 0.5),
#                panel.grid.major = element_line(colour="transparent"),
#                plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
#                legend.title = element_text(size=14))
# 
# g <- g + labs(title = paste("Geographic distribution to compare observation\nand model prediction for", taxon),
#               # subtitle = paste("Random Forest:", paste(env.fact, collapse = " ", sep = " ")), #, "- page", j),
#               x = "",
#               y = "",
#               size = "Probability of\noccurrence",
#               color = "Observation")
# 
# # Configure legends and scales
# g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
#                 # alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
#                 color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
# g <- g + scale_color_manual(values=c(absent = "#c2141b", present = "#007139"), labels=c("Absence", "Presence"))
# g <- g + scale_shape_identity() # Plot the shape according to the data
# print(g)
# ggsave(fileName, path = dir.output, width = 10)
# 
# # dev.off()
# 
# 
# ## ---- Plot taxon-specific response to env fact ----
# 
# taxon.label <- sub("_", " ", occ.taxa)
# 
# # tried to automatize the plotting ... 
# # i <- 1
# # fact <- env.fact[i]
# 
# fileName <- "ALL_TaxResp_vsOneFact.pdf"
# pdf(paste0(dir.output,fileName), paper = 'special', width = 11, onefile = TRUE)
# 
# 
# g <- ggplot(data = plot.data, aes(x = saprobic.cond, y = pred, color = testData$Occurrence.Hydraenidae))
# # !! again specific to the name of variables, it's not general ####
# 
# g <- g + geom_point(alpha = 0.25)
# g <- g + theme_bw(base_size=15)
# 
# g <- g + labs(title = taxon.label,
#               y = "Probability of occurrence",
#               color = "Observation")
# g <- g + scale_color_manual(name = "Observation", values=c(absent = "#FF0000", present = "#0077FF"), labels = c("Absence", "Presence"))
# 
# g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
# g <- g + ylim(0,1)
# print(g)
# 
# dev.off()

