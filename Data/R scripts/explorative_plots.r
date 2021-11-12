## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Explore processed data ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- November 02, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Import libraries and data ----

# Load packages
if ( !require(dplyr) ) { install.packages("dplyr"); library("dplyr") } # to sort, join, merge data
if ( !require(ggplot2) ) { install.packages("ggplot2"); library("ggplot2") } # to do nice plots
if ( !require(skimr) ) { install.packages("skimr"); library("skimr") } # to show key descriptive stats
if ( !require(visdat) ) { install.packages("visdat"); library("visdat") } # to visualize missing data in dataset
if ( !require(sf) ) { install.packages("sf"); library("sf") } # to read GIS data (shape files)
if ( !require(caret) ) { install.packages("caret"); library("caret") } # comprehensive framework to build machine learning models
if ( !require(ggpubr) ) { install.packages("ggpubr"); library("ggpubr") } # to make nice arrangement of nice plots
if ( !require(corrplot) ) { install.packages("corrplot"); library("corrplot") } # to make nice arrangement of nice plots
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") } # to arrange multiple plots on a page

# Check and set working directory
getwd() # show working directory

# Free workspace
rm(list=ls())
graphics.off()

# Define directory and files
dir.env.data      <- "../Processed data/Environmental data/"
dir.inv.data      <- "../Processed data/Invertebrate data/"
dir.output        <- "../../Analysis/Plots/Explorative plots/"

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
    data.env      <- read.delim(paste(dir.env.data,file.env.BDM.data,sep=""),header=T,sep="\t", na = c("<Null>", "NA"), stringsAsFactors=T)
    file.prefix   <- "BDM_"
    
} else {
    
    data.inv          <- read.delim(paste(dir.inv.data,file.inv.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)
    data.env          <- read.delim(paste(dir.env.data,file.env.data,sep=""),header=T,sep="\t", na = c("<Null>", "NA"), stringsAsFactors=T)
    file.prefix       <- "All_"
    
}

## ----  Load functions ----

# source("model_functions.r")
source("../../Analysis/R scripts/plot_functions.r")

## ---- Prepare datasets ----

# make list of taxa
cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
list.taxa <- colnames(data.inv)[cind.taxa]

# replace "0" and "1" by "absent" and "present" and convert them to factors
for (i in cind.taxa ) {
    data.inv[which(data.inv[,i] == 0),i] <- "absent"
    data.inv[which(data.inv[,i] == 1),i] <- "present"
    data.inv[,i] = as.factor(data.inv[,i])
}

# construct main dataset
data <- data.env %>%
    left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
dim(data)

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

# print information
cat("We have", length(prio) + length(explo), "environmental factors to explore.\n",
    length(prio), "of them are interesting (but can still be highly correlated)")

# decide to explore the selected (Bogdan's) factors or the list "priority"
select = T    

if (select == F){
    
    env.fact <- prio
    
} else if (select == T){
    
    env.fact <- c("temperature", # Temp
                  "velocity", # FV
                  "A10m", # A10m
                  "cow.density", # LUD
                  "IAR", # IAR
                  "urban.area", # Urban
                  "FRI", # FRI
                  "bFRI", # bFRI
                  "width.variability")# , # WV
                  # "bed.modification",
                  # "morphology",
                  # "A.EDO",
                  # "F.EDO",
                  # "ARA.fraction",
                  # "agri.land",
                  # "Slope",
                  # "fields",
                  # "saprobic.cond",
                  # "normcov.mobile.blocks",
                  # "normcov.coarse.inorganic.sediments",
                  # "normcov.gravel",
                  # "normcov.sand.silt",
                  # "normcov.fine.sediments")
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

# check if env.fact selected are really factors in data.env
env.fact <- env.fact[which(env.fact %in% colnames(data.env))]

no.env.fact <- length(env.fact)
print(no.env.fact)

## ---- Plot vizualisation of duplicates ----

# count nb of sites and samples
n.sites <- length(unique(data.inv$SiteId))
n.samp <- length(unique(data.inv$SampId))
dup <- duplicated(data.inv$SiteId)
sites.dup <- data.inv$SiteId[dup]

data.inv.dup <- data.inv[which(data.inv$SiteId %in% sites.dup), ] # exclude sites with only 1 sampling
data.inv.dup <- data.inv.dup[order(data.inv.dup$SiteId),]

if( BDM == F){ data.inv.dup <- data.inv.dup[data.inv.dup$MonitoringProgram!="BDM",]} #exclude BDM sites


file.name <- paste0(file.prefix, "DuplicatedSamples.pdf")

if ( !file.exists(paste0(dir.output,file.name)) ){ # don't print if file already exists

    pdf(paste0(dir.output, file.name), paper = 'special', 
        width = 17, height = 13,
        onefile = TRUE)
    
    dnrep <- as.data.frame(cbind("nrep"=sort(table(data.inv.dup$SiteId)),"SiteId"=names(sort(table(data.inv.dup$SiteId)))))
    plot(dnrep$nrep)
    
    tab.nrep <- table(dnrep$nrep) # nr. of replicates, e.g. 4 sites with 8 replicates
    
    data.inv.dup <- left_join(data.inv.dup,dnrep,by="SiteId")
    data.inv.dup$nrep <- as.numeric(data.inv.dup$nrep)
    data.inv.dup <- data.inv.dup[order(data.inv.dup$nrep),]
    
    cols <- rainbow(8,start=0.1)
    plot(data.inv.dup$X,data.inv.dup$Y,pch=16,col=cols[data.inv.dup$nrep],main="sites with repeated samplings")
    legend("topleft",legend=sort(unique(data.inv.dup$nrep)),pch=16,col=cols[unique(data.inv.dup$nrep)])
    
    sdup <- unique(data.inv.dup$SiteId)
    plot(1:length(sdup),type="n",ylim=range(data.inv.dup$Year),xlab="Site",ylab="Year")
    cols <- 1:8
    for(i in 1:length(sdup)){
        years.i <- data.inv.dup$Year[data.inv.dup$SiteId==sdup[i]]
        points(rep(i,length(years.i)),years.i,col=cols[1+i%%10])
        segments(x0=i,y0=min(years.i),y1=max(years.i),col=cols[1+i%%10])
    }

    dev.off()
}

## ---- Produce and plot summ stat ----

# produce table with main infos
df.summaries <- skim(data.env) # produce dataframe containing descriptive statistics for each column
df.summaries.prio <- skim(data.env[,env.fact])
View(df.summaries.prio)

# other way to produce a table to visualize missing data
vis_miss(data.env[,env.fact], cluster = FALSE, warn_large_data = FALSE)

# print table with summary stat in a pdf
temp.df <- as.data.frame(df.summaries.prio)
c.ind <- which(lapply(temp.df,is.numeric) == T)
for (i in c.ind) {
    temp.df[,i] <- round(temp.df[,i], digits = 2)
}

# set pdf size depending of number of env. fact
if( BDM == T){
    dim.pdf <- c("width" = 17, "height" = 17)
} else {
    dim.pdf <- c("width" = 17, "height" = 13)
}


file.name <- paste0(no.env.fact,"envfact_", "Summaries.pdf")

pdf(paste0(dir.output, file.prefix, file.name), paper = 'special', 
    width = dim.pdf["width"], height = dim.pdf["height"],
    onefile = TRUE)
grid.table(temp.df)
dev.off()    


## ---- Plot correlation matrix ----

# set pdf size depending of number of env. fact
if( select == T){
    dim.pdf <- c("width" = 13, "height" = 13)
} else {
    dim.pdf <- c("width" = 20, "height" = 20)
}

ptm <- proc.time() # to calculate time of pdf production

file.name <- paste0(no.env.fact,"envfact_", "CorrMatrix.pdf")

pdf(paste0(dir.output, file.prefix, file.name), paper = 'special', 
    width = dim.pdf["width"], height = dim.pdf["height"],
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


## ---- Plot env. fact. on Swiss map ----

ptm <- proc.time() # to calculate time of pdf production

file.name <- paste0(no.env.fact,"envfact_", "OnMap.pdf")

if ( !file.exists(paste0(dir.output, file.prefix, file.name)) ){ # don't print if file already exists
    
    map.inputs.d <- map.inputs(dir.env.data = dir.env.data, data.env = data.env)
    
    map.env.fact(inputs = map.inputs.d, env.fact = env.fact, data.env = data.env, env.explan = env.explan,
                 dir.output = dir.output, file.prefix = file.prefix, file.name = file.name)
}

print("Producing PDF time:")
print(proc.time()-ptm)
# Producing All_pdf -> 1 min
# Producing BDM_pdf -> 2 min

## ---- Plot taxa vs env. fact. (data visualization) ----

ptm <- proc.time() # to calculate time of pdf production

# produce plots in a list
list.plots <- plot.data.envvstax(data = data, env.fact = env.fact, list.taxa = list.taxa)

# print the plots in a pdf file
file.name <- paste0(file.prefix, no.env.fact,"envfact_EnvFactvsTax.pdf")

if ( !file.exists(paste0(dir.output,file.name)) ){ # don't print if file already exists
    print.pdf.plots(list.plots = list.plots, dir.output = dir.output, file.name = file.name)
}

print("Producing PDF time:")
print(proc.time()-ptm)
# Producing All_pdf -> 1 min
# Producing BDM_pdf -> 4 min
