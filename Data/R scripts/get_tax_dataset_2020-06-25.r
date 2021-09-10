## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Get taxa dataset ----
##
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 29, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Directory and file and date  definitions ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls()) # delete all the objects and variables in the workspace

dir.data                <- "../Processed data/Invertebrate data/"
dir.output              <- "../Processed data/Invertebrate data/"
dir.plot                <- "../../Analysis/Plots/Explorative plots/"

file.inv.data           <- "invertebrates_wide_2020-06-25.dat"      # Output from R scripts/convert_invertebrates_data.r       
file.tax                <- "invertebrates_taxonomy_2020-06-25.dat"  # Output from R scripts/convert_invertebrates_data.r

file.FreshEco        <- "taxa_CH_FreshEco.csv"
file.FreshEco.EU     <- "taxa_EU_FreshEco.csv"
dir.database         <- "../Original data/"

d <- "2020-06-25" # date of the exportation of the MIDAT file database

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Packages and functions ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("get.tax.level.r")
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.inv                <- read.delim(paste(dir.data,file.inv.data,sep=""),header=T,sep="\t",stringsAsFactors=FALSE)           
data.tax                <- read.table(paste(dir.data,file.tax,sep=""),header=T,sep="\t", stringsAsFactors=F)

FreshEco.CH                <- read.csv(paste(dir.database,file.FreshEco,sep=""),header=TRUE,sep=";", stringsAsFactors=FALSE)
FreshEco.EU                <- read.csv(paste(dir.database,file.FreshEco.EU,sep=""),header=TRUE,sep=";", stringsAsFactors=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data.inv ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.inv <- data.inv[which(!is.na(data.inv$IBCH.index) | data.inv$MonitoringProgram == "BDM"),]   # remove data from the monitoring data set that have been collected before standard procedures
data.inv <- data.inv[which(data.inv$Year>2009),] # remove data sampled before 2009
data.inv <- data.inv[-which(data.inv$Canton == "Hors frontières"),] # remove data outside of the border (that are only creating NAs in env.fact later)
colnames(data.inv)[which(grepl("Salmo", colnames(data.inv)))]  # remove columns with fish
if(length(which(grepl("Salmo", colnames(data.inv))))>0){
  data.inv <-  data.inv[,-which(grepl("Salmo", colnames(data.inv))) ]
} 

dim(data.inv) # 4458 2083


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get sampling window and remove samples made out of official sampling window (class 3) ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.inv$Sampling.window <- NA # create a new column for sampling window
i.altitude <- match("Altitude",colnames(data.inv))
data.inv <- data.inv[,c(1:i.altitude,ncol(data.inv),(i.altitude+1):(ncol(data.inv)-1))] # rearrange columns

n.obs <- dim(data.inv)[1] # number of observations/lines

# We classify the samples according to their altitude and date for standardization (Stucki et al, 2019):
# Good sampling window: 1
# Buffer time because of specific hydrological situations: 2
# Out of official sampling window: 3

for (j in 1:n.obs) {
  
  if(data.inv$Altitude[j] >= 190 & data.inv$Altitude[j] <= 600){
    if(data.inv$Month[j] == 3){
      data.inv$Sampling.window[j] <- 1
    }
    else if(data.inv$Month[j] == 2 & data.inv$Day[j] >= 16 | data.inv$Month[j] == 4 & data.inv$Day[j] <= 15){
      data.inv$Sampling.window[j] <- 2
    }
    else {
      data.inv$Sampling.window[j] <- 3
    }
  }
  
  else if(data.inv$Altitude[j] >= 601 & data.inv$Altitude[j] <= 1000){
    if(data.inv$Month[j] == 4){
      data.inv$Sampling.window[j] <- 1
    }
    else if(data.inv$Month[j] == 3 & data.inv$Day[j] >= 16 | data.inv$Month[j] == 5 & data.inv$Day[j] <= 15){
      data.inv$Sampling.window[j] <- 2
    }
    else {
      data.inv$Sampling.window[j] <- 3
    }
  }
  
  else if(data.inv$Altitude[j] >= 1001 & data.inv$Altitude[j] <= 1400){
    if(data.inv$Month[j] == 4 & data.inv$Day[j] >= 16 | data.inv$Month[j] == 5 & data.inv$Day[j] <= 15){
      data.inv$Sampling.window[j] <- 1
    }
    else if(data.inv$Month[j] == 4 & data.inv$Day[j] <= 15 | data.inv$Month[j] == 5 & data.inv$Day[j] >= 16){
      data.inv$Sampling.window[j] <- 2
    }
    else {
      data.inv$Sampling.window[j] <- 3
    }
  }
  
  else if(data.inv$Altitude[j] >= 1401 & data.inv$Altitude[j] <= 1800){
    if(data.inv$Month[j] == 5){
      data.inv$Sampling.window[j] <- 1
    }
    else if(data.inv$Month[j] == 4 & data.inv$Day[j] >= 16 | data.inv$Month[j] == 6 & data.inv$Day[j] <= 15){
      data.inv$Sampling.window[j] <- 2
    }
    else {
      data.inv$Sampling.window[j] <- 3
    }
  }
  
  else {
    if(data.inv$Month[j] == 6){
      data.inv$Sampling.window[j] <- 1
    }
    else if(data.inv$Month[j] == 5 & data.inv$Day[j] >= 16 | data.inv$Month[j] == 7 & data.inv$Day[j] <= 15){
      data.inv$Sampling.window[j] <- 2
    }
    else {
      data.inv$Sampling.window[j] <- 3
    }
  }
}

table(data.inv$Sampling.window)
# 1    2    3 
# 2609  472 1377

# Save raw dataset and remove class 3 in main dataset
data.inv.raw <- data.inv # not outputed in a file
data.inv <- data.inv[which(data.inv$Sampling.window != 3),]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomic harmonization  #### 
# -> this was discussed with Astrid, and then corrected June 2017 # xxx pv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Note: identify old names, typos and so on
## Note: the trait database (directly from the internet) should have the correct names, so we harmonize with this

# names in Freshwater Ecology database

# function to capitalize the first letter of a string
Cap <- function(x) {
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
}
# CH database
true.taxa.CH <- FreshEco.CH[,"Taxon"]                 # xxx nis 01.04.2021 corrected, this was EU before
true.taxa.CH <- sub(" Gen.", "", true.taxa.CH)
true.taxa.CH <- sub(" sp.", "", true.taxa.CH)
true.taxa.CH <- sub(" ", "_", true.taxa.CH)
taxa <- NULL
for(i in 1:length(true.taxa.CH)){
  t <- unlist(strsplit(true.taxa.CH[i], " "))[1]
  taxa <- c(taxa, t)
}
true.taxa.CH <- taxa
true.taxa.CH <- sub("]", "", true.taxa.CH)
true.taxa.CH <- sub("[[]", "", true.taxa.CH)
true.taxa.CH <- sub("ord:", "", true.taxa.CH)
true.taxa.CH <- sub("Stamm:", "", true.taxa.CH)
true.taxa.CH <- sub("Kl:", "", true.taxa.CH)
true.taxa.CH <- sub("Ph:", "", true.taxa.CH)
true.taxa.CH <- sub("UKl:", "", true.taxa.CH)
true.taxa.CH <- sub("UOrd:", "", true.taxa.CH)
true.taxa.CH <- tolower(true.taxa.CH)
taxa <- NULL
for (t in true.taxa.CH){
  taxa <- c(taxa,Cap(t))
}
true.taxa.CH <- taxa

# EU database
true.taxa.EU <- FreshEco.EU[,"Taxon"]
true.taxa.EU <- sub(" Gen.", "", true.taxa.EU)
true.taxa.EU <- sub(" sp.", "", true.taxa.EU)
true.taxa.EU <- sub(" ", "_", true.taxa.EU)
taxa <- NULL
for(i in 1:length(true.taxa.EU)){
  t <- unlist(strsplit(true.taxa.EU[i], " "))[1]
  taxa <- c(taxa, t)
}
true.taxa.EU <- taxa
true.taxa.EU <- sub("]", "", true.taxa.EU)
true.taxa.EU <- sub("[[]", "", true.taxa.EU)
true.taxa.EU <- sub("ord:", "", true.taxa.EU)
true.taxa.EU <- sub("Stamm:", "", true.taxa.EU)
true.taxa.EU <- sub("Kl:", "", true.taxa.EU)
true.taxa.EU <- sub("Ph:", "", true.taxa.EU)
true.taxa.EU <- sub("UKl:", "", true.taxa.EU)
true.taxa.EU <- sub("UOrd:", "", true.taxa.EU)
true.taxa.EU <- tolower(true.taxa.EU)
taxa <- NULL
for (t in true.taxa.EU){
  taxa <- c(taxa,Cap(t))
}
true.taxa.EU <- taxa

# names in our monitoring data
my.taxa <- colnames(data.inv)[which(grepl("Occurrence.",colnames(data.inv)))] 
my.taxa <- sub("Occurrence.", "", my.taxa)

# names of our monitoring data that are not in the Freshwater ecology 
### Note: this includes all genus names
not.true.taxa <- sort(my.taxa[-which(my.taxa %in% true.taxa.EU)])
length(not.true.taxa)

# bad species names:
not.true.taxa[grepl("_",not.true.taxa)]  #xxx nis 01.04.2021 to be clarified! 

fix.matrix <- matrix(NA, nrow=length(not.true.taxa), ncol=4)
colnames(fix.matrix) <- c("old_name", "new_name", "notes", "ref")
fix.matrix[,"old_name"] <- not.true.taxa


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomic resolution within programs ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# original data from the MIDAT database
program <- unique(data.inv[,"MonitoringProgram"])
program 
t <- get_tax_level(data.inv=data.inv, program=program, data.tax=data.tax)

cat("program", "n.samples","n.sites","\n") 
for(p in program)
{
  pdt <- data.inv[which(data.inv$MonitoringProgram==p),]
  cat(p, nrow(pdt), length(unique(pdt$SiteId)), "\n") 
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data.tax ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## add a column with the full name of taxa, and their tax resolution to the data.tax
data.tax$taxon <- NA
data.tax$lowest.level <- NA
for(i in 1:nrow(data.tax))
{
  if(data.tax[i,"Species"]!="")
  {
    data.tax$taxon[i] <- paste(data.tax[i,"Genus"],data.tax[i,"Species"],sep="_" )
    data.tax$lowest.level[i] <- "species"
  } else {
    
    if(data.tax[i,"Genus"]!="")
    {
      data.tax$taxon[i] <- paste(data.tax[i,"Genus"])
      data.tax$lowest.level[i] <- "genus"
    } else {
      if(data.tax[i,"Family"]!="")
      {
        data.tax$taxon[i] <- paste(data.tax[i,"Family"])
        data.tax$lowest.level[i] <- "family"
      } else {
        if(data.tax[i,"Order"]!="")
        {
          data.tax$taxon[i] <- paste(data.tax[i,"Order"])
          data.tax$lowest.level[i] <- "order"
        } else{
          if(data.tax[i,"Class"]!="")
          {
            data.tax$taxon[i] <- paste(data.tax[i,"Class"])
            data.tax$lowest.level[i] <- "class"
          } else {
            data.tax$taxon[i] <- paste(data.tax[i,"Phylum"])
            data.tax$lowest.level[i] <- "phylum"
          }
        }
      }
    }
  }
}


# xxx nis 01.04.2021: this is a bit weird, but used to select a minimum taxonomic level and include coarser levels further below:

## add higher orders to the lower level orders in data.tax, in case these lower orders are empty
# Note: the table now contains the specified taxonomic level, or coarser level
ind <- which(data.tax[,"Class"]=="")
data.tax[ind,"Class"] <- data.tax[ind,"Phylum"] 

ind <- which(data.tax[,"Order"]=="")
data.tax[ind,"Order"] <- data.tax[ind,"Class"] 

ind <- which(data.tax[,"Family"]=="")
data.tax[ind,"Family"] <- data.tax[ind,"Order"] 

ind <- which(data.tax[,"Genus"]=="")
data.tax[ind,"Genus"] <- data.tax[ind,"Family"] 

# The column with species, should get the full name (genus_species), or whatever is higher
data.tax[,"Species"] <- data.tax[,"taxon"] 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get dataset and taxa available ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Species level - BDM dataset ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get dataset
data.BDM <- data.inv[which(data.inv$MonitoringProgram=="BDM"),]
nrow(data.BDM)
length(unique(data.BDM$SiteId))

# identify which taxa are present in the data set
occ.taxa.names <- colnames(data.BDM[,grepl("Occurrence", colnames(data.BDM))])
present.ind <- which(colSums(data.BDM[,occ.taxa.names])>0)
present.occ.taxa <- occ.taxa.names[present.ind]
present.taxa <- substring(occ.taxa.names[present.ind],nchar("Occurrence.")+1) 
length(present.taxa)

# identify where the present.taxa are in data.tax
available.ind <- which(data.tax$taxon %in% present.taxa)
data.tax.BDM <- data.tax[available.ind,]
nrow(data.tax.BDM) #350

# get data core (without any taxa observations)
ind <- c(which(grepl("Occurrence", colnames(data.BDM))), which(grepl("Abundance", colnames(data.BDM))))
data.core <- data.BDM[,-ind]
# colnames(data.core)

# get data of taxa that occur
data.occ <- data.BDM[,which(colnames(data.BDM) %in% present.occ.taxa)]
data.BDM <- cbind(data.core, data.occ)
rownames(data.BDM) <- c(1:nrow(data.BDM))
dim(data.BDM)
cat("no of sites: ",length(unique(data.BDM$SiteId)),"\n")
cat("no of samples: ",length(unique(data.BDM$SampId)),"\n")
cat("data fom years: ",paste(sort(unique(data.BDM$Year))),"\n")
cat("no of taxa: ",ncol(data.occ),"\n")

# 1b. Family level - BDM.fam dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get dataset
data.BDM.fam <- data.BDM

# get data core (without any taxa observations)
ind <- c(which(grepl("Occurrence", colnames(data.BDM.fam))), which(grepl("Abundance", colnames(data.BDM.fam))))
data.core <- data.BDM.fam[,-ind]
# colnames(data.core)

## merge finer taxonomic levels to the family level
# identify names of potential families or coarser level
fam <- unique(data.tax[, "Family"])
length(fam)  
# create dataset with taxa merged at family of coarser level
new.data <- matrix(NA, nrow=nrow(data.BDM.fam), ncol=length(fam))
colnames(new.data) <- fam

for (tn in fam){ #tn: taxa name; 
  # tn = fam[1]
  # cat(tn," ")  
  tn.ind <- which(data.tax[,"Family"]==tn)
  inv.ind <- which(colnames(data.BDM.fam) %in% c(paste("Occurrence.",data.tax[tn.ind, "taxon"], sep="")))
  if(length(inv.ind)>1) {
    new.data[,tn] <- ifelse(rowSums(data.BDM.fam[,inv.ind])>0, 1, 0)
  }else{
    if(length(inv.ind)==1){
      new.data[,tn] <- data.BDM.fam[,inv.ind]
    }   
  }  
}

colnames(new.data) <- paste("Occurrence.",colnames(new.data), sep="") 
# identify which taxa are present
present.ind <- which(colSums(new.data,na.rm=T)>0)
new.data <- new.data[,present.ind]
present.taxa <- substring(colnames(new.data),nchar("Occurrence.")+1)
length(present.taxa)

# join with core data
data.BDM.fam <- cbind(data.core, new.data)
dim(data.BDM.fam)

# identify where the present.taxa are in data.tax
available.ind <- which(data.tax$taxon %in% present.taxa)
data.tax.BDM.fam <- data.tax[available.ind,]
nrow(data.tax.BDM.fam)
not.available <- present.taxa[which(!(present.taxa %in% data.tax$taxon))]
not.available # all taxa are available in the taxonomic dictionary, except:Ptilocolepidae #xxx nis 1.4.21 all seem available

cat("no of sites: ",length(unique(data.BDM.fam$SiteId)),"\n")
cat("no of samples: ",length(unique(data.BDM.fam$SampId)),"\n")
cat("data fom years: ",paste(sort(unique(data.BDM.fam$Year))),"\n")
cat("no of taxa: ",length(present.taxa),"\n")


# 2. Family level - all monitoring programs dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get dataset
data.ALL <- data.inv
nrow(data.ALL)
length(unique(data.ALL$SiteId))

# identify repeated sampling
data.ALL$one <- rep(1,nrow(data.ALL))
data.aggr <- aggregate(data.ALL$one, by=list(data.ALL$SiteId), sum)
rep.sites <- data.aggr[which(data.aggr$x>1), "Group.1"]  # names of repeated sites
mon.program <- unique(data.ALL$MonitoringProgram)
sites <- list(NULL)
i <- 0
for(m in mon.program){
  i <- i+1
  sites[[i]] <- unique(data.ALL[which(data.ALL$MonitoringProgram==m),"SiteId"])
  names(sites)[i] <- m
}
no.rep.sites <- NA
for(i in 1:length(mon.program)){
  no.rep.sites[i] <- length(which(rep.sites %in% unlist(sites[mon.program[i]])))
}
names(no.rep.sites) <- mon.program
no.rep.sites
data.ALL <- data.ALL[,-which(colnames(data.ALL)=="one")]

# get data core (without any taxa observations)
ind <- c(which(grepl("Occurrence", colnames(data.inv))), which(grepl("Abundance", colnames(data.ALL))))
data.core <- data.ALL[,-ind]
# colnames(data.core)

## merge finer taxonomic levels to the family level
# identify names of potential families or coarser level
fam <- unique(data.tax[, "Family"])
length(fam)  
# create dataset with taxa merged at family of coarser level
new.data <- matrix(NA, nrow=nrow(data.ALL), ncol=length(fam))
colnames(new.data) <- fam

for (tn in fam){
  tn.ind <- which(data.tax[,"Family"]==tn)
  inv.ind <- which(colnames(data.ALL) %in% c(paste("Occurrence.",data.tax[tn.ind, "taxon"], sep="")))
  if(length(inv.ind)>1) {
    new.data[,tn] <- ifelse(rowSums(data.ALL[,inv.ind])>0, 1, 0)
  }else{
    if(length(inv.ind)==1){
      new.data[,tn] <- data.ALL[,inv.ind]
    }   
  }  
}

colnames(new.data) <- paste("Occurrence.",colnames(new.data), sep="") 
# identify which taxa are present
present.ind <- which(colSums(new.data)>0)
new.data <- new.data[,present.ind]
present.taxa <- substring(colnames(new.data),nchar("Occurrence.")+1)
length(present.taxa)

# join with core data
data.ALL <- cbind(data.core, new.data)
dim(data.ALL) 


# identify where the present.taxa are in data.tax #xxx now on family or higher level nis 1.4.21
## available.ind <- which(data.tax$taxon %in% present.taxa)
available.ind <- match(present.taxa,data.tax$Family)
data.tax.ALL <- data.tax[available.ind,]
nrow(data.tax.ALL)
not.available <- present.taxa[which(!(present.taxa %in% data.tax$Family))] #xxx fam level
not.available # in data.tax$taxon, because they have always been identified further, 

# xxx nis 1.4.21: this is now resolved
# we have to add these to the data.tax as individual rows
##### add taxa to data.tax
#"Haemopidae"   "Haplotaxidae" "Pyralidae"   


# xxx nis 11.02.21 not.available was empty, therefore I included an if clause

if(length(not.available)>0){
  rind <- NULL
  level <- NULL
  for (i in 1:length(not.available)){
    if(any(data.tax$Phylum==not.available[i])) {
      rind[i] <- which(data.tax$Phylum==not.available[i])[1]
      level[i] <-"Phylum"
    }else{
      if(any(data.tax$Class==not.available[i])) {
        rind[i] <- which(data.tax$Class==not.available[i])[1]
        level[i] <-"Class"
      }else{
        if(any(data.tax$Order==not.available[i])) {
          rind[i] <- which(data.tax$Order==not.available[i])[1]
          level[i] <-"Order"
        }else{
          if(any(data.tax$Family==not.available[i])) {
            rind[i] <- which(data.tax$Family==not.available[i])[1]
            level[i] <-"Family"
          }
        }}}}
  rind
  level
  
  
  # make matrix for additional rows
  rows <- matrix(NA, nrow=length(rind), ncol=ncol(data.tax.ALL))
  colnames(rows) <- colnames(data.tax.ALL)
  
  for (i in 1:length(rind)){
    if (level[i]=="Phylum"){
      rows[i,"Phylum"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Class"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Order"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Family"] <- data.tax[rind[i], "Phylum" ]                         
      rows[i,"Genus"] <- data.tax[rind[i], "Phylum" ] 
      rows[i,"Species"] <- data.tax[rind[i], "Phylum" ]                            
      rows[i,"taxon"] <- data.tax[rind[i], "Phylum" ]    
      rows[i,"lowest.level"] <- "Phylum"
    }else{if (level[i]=="Class"){
      rows[i,"Phylum"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Class"] <- data.tax[rind[i], "Class"]
      rows[i,"Order"] <- data.tax[rind[i], "Class"]
      rows[i,"Family"] <- data.tax[rind[i], "Class" ]                         
      rows[i,"Genus"] <- data.tax[rind[i], "Class" ] 
      rows[i,"Species"] <- data.tax[rind[i], "Class" ]                            
      rows[i,"taxon"] <- data.tax[rind[i], "Class" ]    
      rows[i,"lowest.level"] <- "Class"     
    }else{if (level[i]=="Order"){
      rows[i,"Phylum"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Class"] <- data.tax[rind[i], "Class"]
      rows[i,"Order"] <- data.tax[rind[i], "Order"]
      rows[i,"Family"] <- data.tax[rind[i], "Order" ]                         
      rows[i,"Genus"] <- data.tax[rind[i], "Order" ] 
      rows[i,"Species"] <- data.tax[rind[i], "Order" ]                            
      rows[i,"taxon"] <- data.tax[rind[i], "Order" ]    
      rows[i,"lowest.level"] <- "Order"      
    }else{if (level[i]=="Family"){
      rows[i,"Phylum"] <- data.tax[rind[i], "Phylum"]
      rows[i,"Class"] <- data.tax[rind[i], "Class"]
      rows[i,"Order"] <- data.tax[rind[i], "Order"]
      rows[i,"Family"] <- data.tax[rind[i], "Family" ]                         
      rows[i,"Genus"] <- data.tax[rind[i], "Family" ] 
      rows[i,"Species"] <- data.tax[rind[i], "Family" ]                            
      rows[i,"taxon"] <- data.tax[rind[i], "Family" ]    
      rows[i,"lowest.level"] <- "Family"     
    }}}}}
  dim(data.tax.ALL)  
  data.tax.ALL <- rbind(data.tax.ALL, rows)  
}


dim(data.tax.ALL)
table(data.tax.ALL$lowest.level)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify mixed taxonomic levels  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Species level - BDM dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
warning.list.BDM <- list(NULL)
i <- 0

tax.levels <- c("genus", "family", "order", "class", "phylum")
names(tax.levels) <- c("Genus", "Family", "Order", "Class", "Phylum")

nsamp <- length(data.BDM$SampId)
threshold <- 0.05
n.threshold <- floor(nsamp*threshold)

for (j in 1:length(tax.levels)){
  tax.ind <- which(data.tax.BDM$lowest.level==tax.levels[j])
  names <- unique(data.tax.BDM[tax.ind,names(tax.levels)[j]])
  for (name in names){
    ind <- which(data.tax.BDM[,names(tax.levels)[j]]==name)
    if(length(ind)>1){
      i <- i+1
      warning.list.BDM[[i]] <- list(
        "case" = name,
        "tax.level" = table(data.tax.BDM[ind,"lowest.level"]),
        "taxa" = NULL
      )
      names(warning.list.BDM)[i] <- name 
      nocc <- colSums(data.BDM[,paste("Occurrence.",data.tax.BDM[ind,"taxon"], sep="")])
      warning.list.BDM[[i]]$taxa <- nocc
      warning.list.BDM[[i]]$keeplower <- FALSE
      if(sum(nocc>=n.threshold)>0){
        if(sum(names(which(nocc>=n.threshold))!=name)>0){  # a lower level taxon exists that is above the threshold
          warning.list.BDM[[i]]$keeplower <- TRUE
        }
      }
      names(warning.list.BDM[[i]]$taxa) <- data.tax.BDM[ind,"taxon"]
    }                                 
  }
}
length(warning.list.BDM)

# write to .txt file
sink(file = paste(dir.output,"mixed_taxonomy_BDM_",d,".txt"), append = FALSE, type = c("output", "message"), split = TRUE)
warning.list.BDM
sink()

# 1b. Taxa level - BDM.fam dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
warning.list.BDM.fam <- list(NULL)
i <- 0

tax.levels <- c("family", "order", "class", "phylum")
names(tax.levels) <- c("Family", "Order", "Class", "Phylum")

nsamp <- length(data.BDM.fam$SampId)
threshold <- 0.05
n.threshold <- floor(nsamp*threshold)

for (j in 1:length(tax.levels)){
  tax.ind <- which(data.tax.BDM.fam$lowest.level==tax.levels[j])
  names <- unique(data.tax.BDM.fam[tax.ind,names(tax.levels)[j]])
  for (name in names){
    ind <- which(data.tax.BDM.fam[,names(tax.levels)[j]]==name)
    if(length(ind)>1){
      cat(name, "\n")
      i <- i+1
      warning.list.BDM.fam[[i]] <- list(
        "case" = name,
        "tax.level" = table(data.tax.BDM.fam[ind,"lowest.level"]),
        "taxa" = NULL
      )
      names(warning.list.BDM.fam)[i] <- name
      nocc <- colSums(data.BDM.fam[,paste("Occurrence.",data.tax.BDM.fam[ind,"taxon"], sep="")])
      warning.list.BDM.fam[[i]]$taxa <- nocc
      warning.list.BDM.fam[[i]]$keeplower <- FALSE
      if(sum(nocc>=n.threshold)>0){
        if(sum(names(which(nocc>=n.threshold))!=name)>0){  # a lower level taxon exists that is above the threshold
          warning.list.BDM.fam[[i]]$keeplower <- TRUE
        }
      }
      names(warning.list.BDM.fam[[i]]$taxa) <- data.tax.BDM.fam[ind,"taxon"]
    }
  }
}
length(warning.list.BDM.fam)

# write to .txt file
sink(file = paste(dir.output,"mixed_taxonomy_BDM_fam_",d,".txt"), append = FALSE, type = c("output", "message"), split = TRUE)
warning.list.BDM.fam
sink()

# 2. Family level - all monitoring programs dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
warning.list.ALL <- list(NULL)
i <- 0

tax.levels <- c("family", "order", "class", "phylum")
names(tax.levels) <- c("Family", "Order", "Class", "Phylum")

nsamp <- length(data.ALL$SampId)
threshold <- 0.05
n.threshold <- floor(nsamp*threshold)

for (j in 1:length(tax.levels)){
  tax.ind <- which(data.tax.ALL$lowest.level==tax.levels[j])
  names <- unique(data.tax.ALL[tax.ind,names(tax.levels)[j]])
  for (name in names){
    ind <- which(data.tax.ALL[,names(tax.levels)[j]]==name)
    if(length(ind)>1){
      if(any(data.tax.ALL[ind,"lowest.level"] %in% c("order", "class", "phylum")))
      {
        i <- i+1
        warning.list.ALL[[i]] <- list(
          "case" = name,
          "tax.level" = table(data.tax.ALL[ind,"lowest.level"]),
          "taxa" = NULL)
        
        names(warning.list.ALL)[i] <- name 
        nocc <- colSums(data.ALL[,paste("Occurrence.",data.tax.ALL[ind,"Family"], sep="")])
        warning.list.ALL[[i]]$taxa <- nocc
        warning.list.ALL[[i]]$keeplower <- FALSE
        if(sum(nocc>=n.threshold)>0){
          if(sum(names(which(nocc>=n.threshold))!=name)>0){  # a lower level taxon exists that is above the threshold
            warning.list.ALL[[i]]$keeplower <- TRUE
          }
        }
        names(warning.list.ALL[[i]]$taxa) <- data.tax.ALL[ind,"taxon"]
      }                                 
    }
  }
}
length(warning.list.ALL)

# write to .txt file
sink(file = paste(dir.output,"mixed_taxonomy_ALL_",d,".txt"),
     append = FALSE, type = c("output", "message"), split = TRUE)
warning.list.ALL
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct final dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Species level - BDM dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some checks
dim(data.BDM)         


nrow(data.tax.BDM)         # all of them should be present in the taxonomic dictionary
table(data.tax.BDM$lowest.level)

c1 <- NULL # taxa where lower level taxa should be kept
c2 <- NULL # taxa where lower level taxa should be aggregated

for(i in 1:length(warning.list.BDM)){
  if(warning.list.BDM[[i]]$keeplower){
    c1 <- c(c1,names(warning.list.BDM[i]))
    names(c1)[length(c1)] <- data.tax[which(data.tax$taxon==names(warning.list.BDM[i]))  ,"lowest.level"]
  } else {
    c2 <- c(c2,names(warning.list.BDM[i]))
    names(c2)[length(c2)] <- data.tax[which(data.tax$taxon==names(warning.list.BDM[i]))  ,"lowest.level"]
  }
}

# case 1: coarser taxonomic level only recorded at a limited number of sites, relative to the finer taxonomic level recordings
# solution: set sites where this coarser taxonomic level is recorded to NA for all finer taxonomic level taxa, if their observation was 0.
# reason: finer level taxa would still be present if they were recorded as 1, if they were recorded as 0 it now becomes NA, implying that they could or could not be there.
### 10% version
# c1 <- c("Baetis", "Ecdyonurus", "Epeorus", "Rhithrogena","Habroleptoides","Habrophlebia",
#         "Chloroperla", "Leuctra", "Nemoura", "Protonemura", 
#         "Brachyptera", "Drusus", "Philopotamus", "Rhyacophila", "Chloroperlidae", "Perlodidae",
#         "Taeniopterygidae", "Limnephilidae" ,"Philopotamidae", "Polycentropodidae", 
#         "Psychomyiidae","Rhyacophilidae","Sericostomatidae" )#, "Planorbidae")
# names(c1) <- c("genus","genus","genus","genus","genus","genus","genus","genus",
#                "genus","genus","genus","genus","genus","genus","family","family","family","family",
#                "family","family","family","family","family")#, "family)
### 5% version
# c1 <- c("Baetis", "Ecdyonurus","Electrogena", "Epeorus", "Rhithrogena","Habroleptoides","Habrophlebia",
#         "Chloroperla","Siphonoperla", "Leuctra", "Nemoura", "Protonemura", "Perla",
#         "Brachyptera",  "Hydropsyche", "Drusus", "Philopotamus", "Rhyacophila", 
#         "Chloroperlidae", "Perlodidae","Perlidae","Taeniopterygidae", 
#         "Limnephilidae" ,"Philopotamidae", "Polycentropodidae", 
#         "Psychomyiidae","Rhyacophilidae","Sericostomatidae","Glossosomatidae" )#, "Planorbidae")
# names(c1) <- c("genus","genus","genus","genus","genus","genus","genus",
#                "genus","genus","genus","genus","genus","genus",
#                "genus","genus","genus","genus","genus",
#                "family","family","family","family",
#                "family","family","family",
#                "family","family","family","family")#, "family)
length(c1) # 41
for (i in 1:length(c1)){
  cat(c1[i], "\n")
  taxa <- names(warning.list.BDM[[which(names(warning.list.BDM)== c1[i])]]$taxa)
  names(taxa) <- data.tax.BDM[which(data.tax.BDM$taxon %in% taxa),"lowest.level"]
  if(names(c1[i])=="genus"){
    taxa.to.keep <- taxa[which(names(taxa)=="species")]
    taxa.to.remove <- taxa[-which(names(taxa)=="species")]
    remove.ind <- which(colnames(data.BDM) %in% paste("Occurrence", taxa.to.remove, sep="."))
  }
  if(names(c1[i])=="family"){
    taxa.to.keep <- taxa[which(names(taxa) %in% c("species", "genus"))]
    taxa.to.remove <- taxa[which(names(taxa)=="family")]
    remove.ind <- which(colnames(data.BDM) %in% paste("Occurrence", taxa.to.remove, sep="."))
  }
  if (length(remove.ind)>0){
    if(length(remove.ind)>1){
      rind <- which(rowSums(data.BDM[, remove.ind])==1)
    }else{
      rind <- which(data.BDM[, remove.ind]==1)
    }
    cind <- which(colnames(data.BDM) %in% paste("Occurrence", taxa.to.keep, sep="."))
    data.BDM[rind,cind] <- ifelse(data.BDM[rind,cind]==1,1,NA)
    # remove coarser level taxon
    data.BDM <- data.BDM[,-remove.ind]
  }
}

# case 2: coarser taxonomic level is recorded most of the time, or overall this group of taxa is very rare
# solution: merge finer level taxa into coarser level taxa.
# reason: We would not gain much of the finer level detail, because of high uncertainty on their absence observations, or because anyway there are very few observations
### 10% version
# c2 <- c("Capnia", "Electrogena", "Siphonoperla", "Amphinemura", "Dinocras", "Perla", "Dictyogenus", "Isoperla", "Perlodes",
#         "Rhabdiopteryx","Ernodes", "Agapetus", "Silo", "Hydropsyche", "Stactobia","Mystacides", "Metanoea",
#         "Wormaldia", "Plectrocnemia", "Polycentropus", "Lype","Tinodes",  "Perlidae", "Ancylus",
#         "Goeridae", "Hydroptilidae", "Lepidostomatidae","Glossosomatidae" )
# names(c2) <- c("genus","genus","genus","genus","genus","genus","genus","genus","genus","genus","genus","genus",
#                "genus","genus","genus","genus","genus","genus","genus","genus","genus","genus",
#                 "genus",  "genus", "family","family","family","family")
### 5% version
# c2 <- c("Capnia", "Amphinemura", "Dinocras",  "Dictyogenus", "Isoperla", "Perlodes",
#         "Rhabdiopteryx","Ernodes", "Agapetus", "Silo", "Stactobia","Mystacides", "Metanoea",
#         "Wormaldia", "Plectrocnemia", "Polycentropus", "Lype","Tinodes",  
#         "Goeridae", "Hydroptilidae", "Lepidostomatidae")
# names(c2) <- c("genus","genus","genus","genus","genus","genus",
#                "genus","genus","genus", "genus","genus","genus",
#                "genus","genus","genus","genus","genus","genus",
#                "family",  "family", "family")
length(c2) # 17
for (c in c2){
  taxa <- names(warning.list.BDM[[which(names(warning.list.BDM)== c)]]$taxa)
  cind <- which(colnames(data.BDM) %in% paste("Occurrence", taxa, sep="."))
  print(paste(c,", taxa to merge: ", length(cind), sep=""))
  # merge taxa together
  if(length(cind)>1){
    for (i in 1:nrow(data.BDM)){
      ifelse(sum(data.BDM[i, cind], na.rm=TRUE) >0, data.BDM[i, paste("Occurrence.group", c, sep=".")] <- 1,data.BDM[i, paste("Occurrence.group", c, sep=".")] <- 0)
    }
  }else{
    data.BDM[, paste("Occurrence.group", c, sep=".")] <-  data.BDM[, cind]
  }
  # remove old taxa
  data.BDM <- data.BDM[,-cind]
}

# some checks
dim(data.BDM)       
# colnames(data.BDM)
start <- min(which(grepl("Occurrence",colnames(data.BDM))))
for(j in start:ncol(data.BDM)){
  for (i in 1:nrow(data.BDM)){
    if(!(is.na(data.BDM[i,j])) & data.BDM[i,j]>1) print(colnames(data.BDM)[j]) 
  }
}


## save final data
filename <- paste(dir.output, "BDM_occ_data_", d,".dat", sep="")
write.table(data.BDM, filename, sep="\t", col.names=TRUE, row.names=FALSE)

# some checks
check <- read.table(filename, header = T, sep="\t",stringsAsFactors=F)
dim(check)

cat("no of taxa data.BDM:", sum(grepl("Occurrence.",colnames(data.BDM))))

# 1b. BDM.fam dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## note: no taxonomic mismatches (mismatches in the BDM occur at family or more detailed taxonomic levels)

# some checks
dim(data.BDM.fam)      
#colnames(data.BDM.fam) 
start <- min(which(grepl("Occurrence",colnames(data.BDM.fam)))) ### to check if there are any values greater than 1 in the dataset
for(j in start:ncol(data.BDM.fam)){
  for (i in 1:nrow(data.BDM.fam)){
    if(!(is.na(data.BDM.fam[i,j])) & data.BDM.fam[i,j]>1) print(colnames(data.BDM.fam)[j]) 
  }
}

## save final data
filename <- paste(dir.output, "BDM_fam_occ_data_", d,".dat", sep="")
write.table(data.BDM.fam, filename, sep="\t", col.names=TRUE, row.names=FALSE)

# some checks
check <- read.table(filename, header = T, sep="\t",stringsAsFactors=F)
dim(check)

cat("no of taxa data.BDM.fam:", sum(grepl("Occurrence.",colnames(data.BDM.fam))))

# 2. Family level - all monitoring programs dataset  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some checks
dim(data.ALL)              
length(present.taxa)         
nrow(data.tax.ALL)           # all of them should be present in the taxonomic dictionary
table(data.tax.ALL$lowest.level)


c1 <- NULL # taxa where lower level taxa should be kept
c2 <- NULL # taxa where lower level taxa should be aggregated

for(i in 1:length(warning.list.ALL)){
  if(warning.list.ALL[[i]]$keeplower){
    c1 <- c(c1,names(warning.list.ALL[i]))
    names(c1)[length(c1)] <- data.tax[which(data.tax$taxon==names(warning.list.ALL[i]))  ,"lowest.level"]
  } else {
    c2 <- c(c2,names(warning.list.ALL[i]))
    names(c2)[length(c2)] <- data.tax[which(data.tax$taxon==names(warning.list.ALL[i]))  ,"lowest.level"]
  }
}

# case 1: coarser taxonomic level only recorded at a limited number of sites, relative to the finer taxonomic level recordings
# solution: set sites where this coarser taxonomic level is recorded to NA for all finer taxonomic level taxa, if their observation was 0.
# reason: finer level taxa would still be present if they were recorded as 1, if they were recorded as 0 it now becomes NA, implying that they could or could not be there.
# c1 <- c("Diptera", "Odonata","Oligochaeta", "Plecoptera", "Hirudinea", "Platyhelminthes")
for (c in c1){
  rind <- which(data.ALL[, paste("Occurrence", c, sep=".")]==1)
  taxa <- names(warning.list.ALL[[which(names(warning.list.ALL)== c)]]$taxa)
  taxa <- taxa[-which(taxa==c)]
  cind <- which(colnames(data.ALL) %in% paste("Occurrence", taxa, sep="."))
  data.ALL[rind,cind] <- ifelse(data.ALL[rind,cind]==1,1,NA)
  # remove coarser level taxon
  data.ALL <- data.ALL[,-which(colnames(data.ALL)==paste("Occurrence", c, sep="."))]
}

# case 2: coarser taxonomic level is recorded most of the time, or overall this group of taxa is very rare
# solution: merge finer level taxa into coarser level taxa.
# reason: We would not gain much of the finer level detail, because of high uncertainty on their absence observations, or because anyway there are very few observations

# c2 <- c("Lepidoptera", "Hemiptera", "Cnidaria","Nematomorpha", "Porifera")
for (c in c2){
  taxa <- names(warning.list.ALL[[which(names(warning.list.ALL)== c)]]$taxa)
  cind <- which(colnames(data.ALL) %in% paste("Occurrence", taxa, sep="."))
  # merge taxa together
  for (i in 1:nrow(data.ALL)){
    ifelse(sum(data.ALL[i, cind], na.rm=TRUE) >0, data.ALL[i, paste("Occurrence.group", c, sep=".")] <- 1,data.ALL[i, paste("Occurrence.group", c, sep=".")] <- 0)
  }
  # remove old taxa
  data.ALL <- data.ALL[,-cind]
}

# some checks
dim(data.ALL)           
#colnames(data.ALL) 
start <- min(which(grepl("Occurrence",colnames(data.ALL)))) 
for(j in start:ncol(data.ALL)){
  for (i in 1:nrow(data.ALL)){
    if(!(is.na(data.ALL[i,j])) & data.ALL[i,j]>1) print(colnames(data.ALL)[j]) 
  }
}

## save final data
filename <- paste(dir.output, "All_occ_data_", d,".dat", sep="")
write.table(data.ALL, filename, sep="\t", col.names=TRUE, row.names=FALSE)


# some checks
check <- read.table(filename, header = T, sep="\t",stringsAsFactors=F)
dim(check)

cat("no of taxa data.ALL:", sum(grepl("Occurrence.",colnames(data.ALL))))

#sort(colnames(data.ALL)[grepl("Occurrence.",colnames(data.ALL))])

# # case 1
# rind <- c( 71, 1328, 1341, 1712, 1765, 1767, 2057, 2063, 2127, 2133, 2145, 2154, 2261, 2273, 2283, 2301) # rows where Diptera where
# taxa <- c("Anthomyiidae","Athericidae","Blephariceridae", "Ceratopogonidae", "Chaoboridae", "Chironomidae", "Culicidae") # within the Diptera
# check[rind, paste("Occurrence",taxa, sep=".")]
# case 2
# length(which(check$Occurrence.group.Porifera==1))  # should be 6 sites
# colnames(check)[which(grepl("Porifera", colnames(check)))]  # should only be the Porifera group, no more Porifera on its own

# 3. Write table with sites/samples information to get env data from rs  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c.ind <- grep("Occurrence.",colnames(data.ALL))

data.forenv <- data.ALL[,-c.ind]

filename <- paste(dir.output, "SitesData_for_RS_", d,".dat", sep="")
write.table(data.forenv, filename, sep="\t", col.names=TRUE, row.names=FALSE)


# 4. Plot prevalence ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Construct dataframe with prevalence and taxonomic level for ALL

cind.taxa <- which(grepl("Occurrence.",colnames(data.ALL)))

Occurrence.taxa <- colnames(data.ALL)[cind.taxa]

data.prev <- data.frame(Occurrence.taxa)
data.prev$Prevalence <- NA
data.prev$Taxonomic.level <- NA
data.prev$Missing.values <- NA

for(i in cind.taxa){
  
  n.missobs <- sum(is.na(data.ALL[,i]))
  n.obs <- sum(!is.na(data.ALL[,i]))
  
  if( n.missobs > 100){ # too many missing observation
    
    cat("sum NA for", colnames(data.ALL)[i], "is", n.missobs, "too many missing obs !! \n")
    
    # sum na for Occurrence.Lumbriculidae is 3863 
    # sum na for Occurrence.Lumbricidae is 3691 
    # sum na for Occurrence.Naididae is 3873 
    # sum na for Occurrence.Tubificidae is 3852
    
  } else if ( 100 > n.missobs && n.missobs > 0 ) { # some missing information, we have to be careful
    
    cat("sum NA for", colnames(data.ALL)[i], "is", n.missobs, "\n")
  }
  
  #fill missing values column
  data.prev[which(data.prev[,1] == colnames(data.ALL)[i]), "Missing.values"] <- n.missobs
  
  # fill prevalence column
  prev <- sum(na.omit(data.ALL[,i]))/n.obs # calculate the of prevalence of this taxon
  data.prev[which(data.prev[,1] == colnames(data.ALL)[i]), "Prevalence"] <- prev
  
  # fill taxonomic level column
  rind <- which(data.tax[,"taxon"] == sub("Occurrence.","", colnames(data.ALL)[i])) # look for taxon in data.tax
  
  if (length(rind) != 0) {
    
    data.prev[which(data.prev[,1] == colnames(data.ALL)[i]), "Taxonomic.level"] <- data.tax[rind,"lowest.level"]
    
  } else {
    cat(sub("Occurrence.","", colnames(data.ALL)[i]), "is not in taxa homogenization \n")
  }
}

# reorder prevalence dataframe in decreasing order
data.prev <- arrange(data.prev, desc(Prevalence))


pdf(paste(dir.plot, "All_Prevalence.pdf"), height = 10, width = 13)
ggplot(data=data.prev, aes( x = 1:dim(data.prev)[1], y = Prevalence, colour = Taxonomic.level))+ 
  geom_point() +
  labs(title = "Prevalence of taxa", subtitle = "test", caption = "other test", 
       x = "Taxa (ordered by decreasing number of occurrences)", y = "Prevalence",
       color = "Taxonomic level") + #, tag = "A") +
  scale_x_discrete(limits = sub("Occurrence.","", data.prev[, "Occurrence.taxa"])) +
  theme(axis.text.x = element_text(angle=90))

dev.off()

## save the table for later use
filename <- paste(dir.output, "All_prevalence_", d,".dat", sep="")
write.table(data.prev, filename, sep="\t", col.names=TRUE, row.names=FALSE)

# Construct dataframe with prevalence and taxonomic level for BDM

cind.taxa <- which(grepl("Occurrence.",colnames(data.BDM)))

Occurrence.taxa <- colnames(data.BDM)[cind.taxa]

data.prev.BDM <- data.frame(Occurrence.taxa)
data.prev.BDM$Prevalence <- NA
data.prev.BDM$Taxonomic.level <- NA
data.prev.BDM$Missing.values <- NA


for(i in cind.taxa){
  
  n.missobs <- sum(is.na(data.BDM[,i]))
  n.obs <- sum(!is.na(data.BDM[,i]))
  
  if( n.missobs > 100){ # too many missing observation
    
    cat("sum NA for", colnames(data.BDM)[i], "is", n.missobs, "too many missing obs !! \n")
    
    # sum na for Occurrence.Lumbriculidae is 3863 
    # sum na for Occurrence.Lumbricidae is 3691 
    # sum na for Occurrence.Naididae is 3873 
    # sum na for Occurrence.Tubificidae is 3852
    
  } else if ( 100 > n.missobs && n.missobs > 0 ) { # some missing information, we have to be careful
    
    cat("sum NA for", colnames(data.BDM)[i], "is", n.missobs, "\n")
  }
  
  #fill missing values column
  data.prev.BDM[which(data.prev.BDM[,1] == colnames(data.BDM)[i]), "Missing.values"] <- n.missobs
  
  # fill prevalence column
  prev <- sum(na.omit(data.BDM[,i]))/n.obs # calculate the of prevalence of this taxon
  data.prev.BDM[which(data.prev.BDM[,1] == colnames(data.BDM)[i]), "Prevalence"] <- prev
  
  # fill taxonomic leel column
  rind <- which(data.tax[,"taxon"] == sub("Occurrence.","", colnames(data.BDM)[i])) # look for taxon in data.tax
  
  if (length(rind) != 0) {
    data.prev.BDM[which(data.prev.BDM[,1] == colnames(data.BDM)[i]), "Taxonomic.level"] <- data.tax[rind,"lowest.level"]
  } else {
    cat(sub("Occurrence.","", colnames(data.BDM)[i]), "is not in taxa homogenization \n")
  }
}

# reorder prevalence dataframe in decreasing order
data.prev.BDM <- arrange(data.prev.BDM, desc(Prevalence))


pdf(paste(dir.plot, "BDM_Prevalence.pdf"), height = 10, width = 13)
ggplot(data=data.prev.BDM, aes( x = 1:dim(data.prev.BDM)[1], y = Prevalence, colour = Taxonomic.level))+ 
  geom_point() +
  labs(title = "Prevalence of taxa", subtitle = "test", caption = "other test", 
       x = "Taxa (ordered by decreasing number of occurrences)", y = "Prevalence",
       color = "Taxonomic level") + #, tag = "A") +
  scale_x_discrete(limits = sub("Occurrence.","", data.prev.BDM[, "Occurrence.taxa"])) +
  theme(axis.text.x = element_text(angle=90))

dev.off()

## save the table for later use
filename <- paste(dir.output, "BDM_prevalence_", d,".dat", sep="")
write.table(data.prev.BDM, filename, sep="\t", col.names=TRUE, row.names=FALSE)


