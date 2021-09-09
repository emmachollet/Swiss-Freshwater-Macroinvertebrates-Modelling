## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Convert river invertebrates data ----
##
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 29, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ---- Libraries and global definitions ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set working directory:
rm(list = ls()) # delete all the objects and variables in the workspace

# directory and file definitions:
dir.orig       <- "../Original data/"

file.orig     <- "MIDAT-data-export-20200625_ALL_2_transm_2.csv"

dir.proc       <- "../Processed data/Invertebrate data/"
file.data.proc <- "invertebrates"
file.taxa.proc <- "invertebrates_taxonomy"
file.meta.proc <- "invertebrates_metadata"

d <- "2020-06-25"  # date when the MIDAT dataset was exported, needed to date output files

# utilities:
source("utilities.r")


# ---- Read invertebrates data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read data:
data <- read.csv(paste(dir.orig,file.orig,sep=""),header=TRUE,sep=";",
                 stringsAsFactors=FALSE,strip.white = TRUE)
data <- tail(data, -2) # delete first 2 lines (concepts names in German and in French)


# replace column names for Swiss Coordinates:
colnames(data)[colnames(data)=="Coord.X..m."] <- "X"
colnames(data)[colnames(data)=="Coord.Y..m."] <- "Y"

data$Abundance <- as.numeric(data$Abundance)
data$Day <- as.numeric(data$Day)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$X <- as.numeric(data$X)
data$Y <- as.numeric(data$Y)

# eliminate data (taxonomic levels) with "Abundance" zero:
data <- data[data$Abundance>0,]

# eliminate data from before 2009:
data <- data[data$Year>"2009",]

# ---- Add SiteId and SampId ----
# ~~~~~~~~~~~~~~~~~~~~

data <- cbind(SiteId=paste("CSCF_",data$Station.nb.,sep=""),
              SampId=paste("CSCF_",data$Station.nb.,"_",
                           data$Year,"-",
                           ifelse(data$Month<10,paste("0",data$Month,sep=""),data$Month),"-",
                           ifelse(data$Day<10,paste("0",data$Day,sep=""),data$Day),
                           sep=""),
              MonitoringProgram=NA,
              data)


# ---- Identify monitoring program ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data$MonitoringProgram[data$Project=="Projet NAWA - Observation nationale de la qualité des eaux de surface"] <- "NAWA TREND"
data$MonitoringProgram[data$Project=="NAWA SPEZ"] <- "NAWA SPEZ"

data$MonitoringProgram[data$Project=="BDM, Indicateur Z9, insectes aquatiques"] <- "BDM"

data[is.na(data)] <- "Other programs"

# check for new Monitoring programs that have to be added:
ind.na <- which(is.na(data$MonitoringProgram))
if(length(ind.na)>0) {
  cat("for ",length(ind.na),"rows the monitoring program is missing")
  data.na <- data[ind.na,]
}


# ---- Correct spelling and taxonomic errors ----
# In unclear cases we chose the version used in the freshwaterecology database
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data$Family <- sub("Cordulegasteridae","Cordulegastridae",data$Family)
data$Family <- sub("Limnaeidae","Lymnaeidae",data$Family)

ind <- which(data$Genus=="Drusus" & data$Species=="monticolus");     data$Species[ind]<-"monticola"  # not clear 
ind <- which(data$Genus=="Alainites" & data$Species=="muticus");     data$Genus[ind]<-"Baetis"       # synonymous
ind <- which(data$Genus=="Nigrobaetis" & data$Species=="niger");     data$Genus[ind]<-"Baetis"       # not clear 
ind <- which(data$Genus=="Dictyogenus" & data$Species=="alpinus");   data$Species[ind]<-"alpinum"   
ind <- which(data$Genus=="Haitia");                                  data$Genus[ind]<-"Physella"     # synonymous
ind <- which(data$Genus=="Physella"& data$Species=="heterostropha"); data$Species[ind]<-"acuta" 
ind <- which(data$Genus=="Radix" & data$Species=="ovata");           data$Species[ind]<-"balthica"   # synonymous   
ind <- which(data$Genus=="Gammarus" & data$Species=="roeseli");      data$Species[ind]<-"roeselii"   # not clear
ind <- which(data$Genus=="Atherix" & data$Species=="marginata");     data$Genus[ind]<-"Ibisia"       # synonymous   
ind <- which(data$Genus=="Lestes" & data$Species=="viridis");        data$Genus[ind]<-"Chalcolestes" # synonymous  

ind <- which(data$Family=="Ancylidae");  data$Genus[ind]<-"Ancylus"; data$Family[ind]<-"Planorbidae"  # can lead to problems with taxonomic homogenization!
ind <- which(data$Genus=="Baetis" & data$Species=="pentaphlebodes"); data$Species[ind]<-"nexus"       # to match with traitdatabase, pentaphlebodes seems the newer name
ind <- which(data$Family=="Helodidae");                              data$Family[ind]<-"Scirtidae" 
ind <- which(data$Genus=="Helodes");                                 data$Genus[ind]<-"Elodes"        # synonymous
ind <- which(data$Genus=="Halesus" & data$Species=="tesselatus");    data$Species[ind]<-"tessellatus" # not clear
ind <- which(data$Family=="Hirudidae");                              data$Family[ind]<-"Hirudinidae" 
ind <- which(data$Genus=="Ephemerella" & data$Species=="ignita");    data$Genus[ind]<-"Serratella"    # not clear
ind <- which(data$Genus=="Perla" & data$Species=="burmeisteriana");  data$Species[ind]<-"abdominalis"

#temp <- data[ind,]

# clarify:
# Ancylidae change to Ancylus 
# Baetis pentaphlebodes change to Baetis nexus?, check with Nadine or Astrid  , yes change to nexus to match with FWE even though penta... seems newer
# Helodidae change to Scirtidae?, yes
# Helodes should be Elodes? yes
# Halesus tesselatus change to Halesus tessellatus yes, but unclear
# Hirudidae (in our data) and Hirudinidae in freshwater ecology, yes, change to Hirudinidae
# Ephemerella ignita change to Serratella ignita? yes but unclear
# Perla burmeisteriana change to Perla abdominalis 
# Tubificida is correct (order level)

# ---- Separate abundance and abundance class ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.abundanceclass <- function(a) { return(ifelse(a<=100,ifelse(a<=11,a,11),ifelse(a<=1000,101,1001))) }

classes <- c(0:11,101,1001)
data$AbundanceClass <- NA
i.abundance <- match("Abundance",colnames(data))
data <- data[,c(1:i.abundance,ncol(data),(i.abundance+1):(ncol(data)-1))] # rearrange columns
samples <- paste(data$SiteId,data$Year,data$Month,data$Day,sep="_")
for ( s in unique(samples) ){
  ind <- samples == s
  ind.not.class <- is.na(match(data$Abundance[ind],classes))
  if ( sum(ind.not.class) == 0 )  # only classes 
  {
    data$AbundanceClass[ind] <- data$Abundance[ind]
    data$Abundance[ind] <- ifelse(data$Abundance[ind]<=10,data$Abundance[ind],NA)
  }
  else
  {
    data$AbundanceClass[ind] <- get.abundanceclass(data$Abundance[ind])
    if ( max(data$Abundance[ind][ind.not.class]) <= 100 )  # some small abundances, classes for large abundances
    {
      data$Abundance[ind] <- ifelse(data$Abundance[ind]!=11 & data$Abundance[ind]<=100,data$Abundance[ind],NA)
    }
  }
}

# ---- Construct taxon name ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data$Taxon <- ifelse((nchar(data$Species)>0 & !is.na(data$Species)),
                     paste(data$Genus,data$Species,sep="_"),
                     ifelse((nchar(data$Genus)>0 & !is.na(data$Genus)),
                            data$Genus,
                            ifelse((nchar(data$Family)>0 & !is.na(data$Family)),
                                   data$Family,
                                   ifelse((nchar(data$Order)>0 & !is.na(data$Order)),
                                          data$Order,
                                          ifelse((nchar(data$Class)>0 & !is.na(data$Class)),
                                                 data$Class,
                                                 data$Phylum)))))
i.species <- match("Species",colnames(data))
data <- data[,c(1:i.species,ncol(data),(i.species+1):(ncol(data)-1))]


# ---- Get vector with lowest level taxa and their level ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LowestLevelTaxon <- data$Taxon

names(LowestLevelTaxon) <- ifelse((nchar(data$Species)>0 & !is.na(data$Species)),
                                  "Species",
                                  ifelse((nchar(data$Genus)>0 & !is.na(data$Genus)),
                                         "Genus",
                                         ifelse((nchar(data$Family)>0 & !is.na(data$Family)),
                                                "Family",
                                                ifelse((nchar(data$Order)>0 & !is.na(data$Order)),
                                                       "Order",
                                                       ifelse((nchar(data$Class)>0 & !is.na(data$Class)),
                                                              "Class",
                                                              "Phylum")))))

lltaxa <- unique(LowestLevelTaxon)
ind <- match(lltaxa,LowestLevelTaxon)
names(lltaxa) <- names(LowestLevelTaxon)[ind]


# ---- Transform factors to character ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fac <- sapply(data,is.factor)
data[fac] <- lapply(data[fac],as.character)


# ---- Aggregate duplicated entries ----   # happens due to "species complexes"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  # we decide if we aggregate to genus or to species level

n <- nrow(data)
dat <- paste(data$SampId,data$Taxon,sep="_")
dup <- unique(dat[duplicated(dat)])
if ( length(dup) > 0 ){
  cat(length(dup),"duplicated entries found\n")
  for ( i in 1:length(dup) )
  {
    ind <- which(paste(data$SampId,data$Taxon,sep="_")==dup[i])  
    if ( length(ind) < 2 ) stop("error in dealing with duplicates")
    
    # # equivalent to if else below, but seems to be equally slow
    # data$Abundance[ind[1]] <- ifelse(sum(is.na(data$Abundance[ind])) == 0,sum(as.numeric(data$Abundance[ind])),NA)
    # data$AbundanceClass[ind[1]] <- ifelse(sum(is.na(data$Abundance[ind])) == 0,get.abundanceclass(data$Abundance[ind[1]]), max(data$AbundanceClass[ind]))
    
    if ( sum(is.na(data$Abundance[ind])) == 0 )  # sum abundances:
    {
      data$Abundance[ind[1]] <- sum(as.numeric(data$Abundance[ind]))
      data$AbundanceClass[ind[1]] <- get.abundanceclass(data$Abundance[ind[1]])
    }
    else  # use largest abundance class:
    {
      data$Abundance[ind[1]] <- NA
      data$AbundanceClass[ind[1]] <- max(data$AbundanceClass[ind])
    }
    data <- data[-ind[-1],]   # eliminate duplicates
    if(i%%100==0) cat(i," duplicates eliminated\n")
  }
  cat(n-nrow(data),"duplicated species records aggregated\n")
}


# --- Add Longitude and Latitude in addition to Swiss coordinates ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coord <- CH.to.WGS(x=as.numeric(data$X),y=as.numeric(data$Y))
data <- cbind(data,coord)


# --- resolve missing entries for taxonomy  -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# missing information will be copied from other rows with same taxon name
# cause: Hydrochidae and Helophoridae same as Hydrophilidae? 
# xxx nis 22.03.2021

ind.m <- which(data$Phylum=="")

if(length(ind.m)>0){
  cat("check taxa names for ",data[ind.m,"Taxon.IBCH"])
  
  for(i in ind.m)
  {
    tax <- data[i,"Taxon.IBCH"]
    ind.tax <- setdiff(which(data$Taxon.IBCH==tax),i)
    if(sum(ind.tax)>0 ) {
      data[i,c("Phylum",	"Class",	"Order",	"Family")] <-  
        data[ind.tax[1],c("Phylum",	"Class",	"Order",	"Family")]
    }
  }
}

# ---- Extract taxonomic classifications ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

taxonomy <- sort(unique(paste(data$Phylum,data$Class,data$Order,data$Family,data$Genus,data$Species,sep="\t")))
Phylum   <- unique(data$Phylum);  Phylum  <- Phylum[Phylum!=""]; Phylum <- Phylum[!is.na(Phylum)]
Class    <- unique(data$Class);   Class   <- Class[Class!=""];   Class <- Class[!is.na(Class)]
Order    <- unique(data$Order);   Order   <- Order[Order!=""];   Order <- Order[!is.na(Order)]
Family   <- unique(data$Family);  Family  <- Family[Family!=""]; Family <- Family[!is.na(Family)]
Genus    <- unique(data$Genus);   Genus   <- Genus[Genus!=""];   Genus <- Genus[!is.na(Genus)]
Species  <- unique(data$Species); Species <- Species[Species!=""]; Species <- Species[!is.na(Species)]


# ---- List summary of data ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Summary of invertebrate data processed:\n")
cat("Number of records: ",nrow(data),"\n")
cat("Number of stations:",length(unique(data$Station.nb.)),"\n")
cat("Time span:         ",min(data$Year),"-",max(data$Year),"\n")
cat("Altitude range:    ",trunc(min(as.numeric(data$Altitude),na.rm=TRUE)),"-",
    trunc(max(as.numeric(data$Altitude),na.rm=TRUE)),"m\n")
cat("Phylum:            ",length(Phylum),"\n")
cat("Class:             ",length(Class),"\n")
cat("Order:             ",length(Order),"\n")
cat("Family:            ",length(Family),"\n")
cat("Genus:             ",length(Genus),"\n")
cat("Species:           ",length(Species),"\n")
no.program.id <- sum(is.na(data$MonitoringProgram))
if ( no.program.id > 0 ) cat("*** no monitoring program identified for",no.program.id,"taxa\n")

cat("Summary of invertebrate data processed:\n",
"Number of records: ",nrow(data),"\n",
"Number of stations:",length(unique(data$Station.nb.)),"\n",
"Time span:         ",min(data$Year),"-",max(data$Year),"\n",
"Altitude range:    ",trunc(min(as.numeric(data$Altitude),na.rm=TRUE)),"-",
    trunc(max(as.numeric(data$Altitude),na.rm=TRUE)),"m\n",
"Phylum:            ",length(Phylum),"\n",
"Class:             ",length(Class),"\n",
"Order:             ",length(Order),"\n",
"Family:            ",length(Family),"\n",
"Genus:             ",length(Genus),"\n",
"Species:           ",length(Species),"\n", file = paste(dir.proc,"summary_processed_WideInvData","_",d,".txt"))


# ---- Produce wide format ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

drop <- c("Sample.number","alternative.abundance",
          "Sample.ID..database.","GEWISS.no.", "Linked.stations", "Main.stations",# xxx ecr 29.06.21: additional columns to drop added
          "Stage","Phylum","Class","Order","Family","Genus","Species",
          "Taxon.on.laboratory.protocol","Taxon.ID","Determining.group", 
          "comment","Taxon.IBCH","Indicator.group", # xxx nis 02.03.21: additional columns to drop added
          "Number.of.IBCH.taxa", "Rating.IBCH.index.", "Rating.SPEAR.index.",
          "MakroIndex", "Rating.MakroIndex.", "Project", "Principal.Institution") # xxx ecr 29.06.21: additional columns to drop added

data.copy <- data
data.copy$Abundance <- ifelse(is.na(data.copy$Abundance),-1,data.copy$Abundance) # encode true NAs by -1
data.copy$Occurrence <- rep(1,nrow(data.copy))
data.wide <- reshape(data=data.copy,
                     timevar="Taxon",
                     v.names=c("Occurrence","Abundance","AbundanceClass"),
                     idvar="SampId",
                     drop=drop,
                     direction="wide")
ind.occurrence     <- which(substr(colnames(data.wide),1,11)=="Occurrence.")
ind.abundance      <- which(substr(colnames(data.wide),1,10)=="Abundance.")
ind.abundanceclass <- which(substr(colnames(data.wide),1,15)=="AbundanceClass.")
for ( i in ind.occurrence ){
  data.wide[,i] <- ifelse(is.na(data.wide[,i]),0,data.wide[,i]) # replace "new" NAs by 0
}
for ( i in ind.abundance ){
  data.wide[,i] <- ifelse(is.na(data.wide[,i]),0,data.wide[,i]) # replace "new" NAs by 0
  data.wide[,i] <- ifelse(data.wide[,i]==-1,NA,data.wide[,i])   # restore true NAs
}
for ( i in ind.abundanceclass ){
  data.wide[,i] <- ifelse(is.na(data.wide[,i]),0,data.wide[,i]) # replace "new" NAs by 0
}
offset <- min(c(ind.occurrence,ind.abundance,ind.abundanceclass))-1
data.wide <- data.wide[,c(1:offset,ind.occurrence,ind.abundance,ind.abundanceclass)]


# ---- Produce wide format where lower level taxa are added to the higher levels ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

levelnames <-  c("Phylum","Class","Order","Family","Genus","Species")
ind.occurrence     <- which(substr(colnames(data.wide),1,11)=="Occurrence.")
data.wide.taxhom <- data.wide[,c(1:max(ind.occurrence))]

for(j in ind.occurrence){
  taxon <- substring(colnames(data.wide)[j],12)
  
  rind <- match(taxon,data$Taxon)
  
  level <- names(lltaxa)[match(taxon,lltaxa)]
  higherlevels <- levelnames[1:(match(level,levelnames)-1)]
  hltaxa <- as.character(data[rind,higherlevels])
  
  ind.hl.lltaxa <- match(hltaxa,lltaxa)
  
  if(sum(!is.na(ind.hl.lltaxa))>0)
  {
    hllltaxa <- lltaxa[ind.hl.lltaxa[!is.na(ind.hl.lltaxa)]]
    cind <- match(hllltaxa, substring(colnames(data.wide.taxhom),12))
    for(k in cind) 
    {
      data.wide.taxhom[,k] <- apply(data.wide.taxhom[,c(j,k)],1,max,na.rm=TRUE)
    }
  }
}

# ---- Collect metadata ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

metadata <- list()

metadata$Title          <- "River Invertebrate Data Switzerland"

metadata$Description    <- "Swiss river invertebrate sampling data from national and cantonal authorities"

metadata$Variable       <- "Invertebrates"

metadata$SystemType     <- "River"

metadata$StartingDate   <- min(data$Year)

metadata$Scope          <- "Observations"

metadata$DataFormat     <- "CSV"

metadata$Owner          <- "CSCF - Centre Suisse de Cartographie  de la Faune, Passage Max.-de-Meuron 6, 2000 Neuchâtel, http://www.cscf.ch"

metadata$DataCurator    <- "Nele Schuwirth <nele.schuwirth@eawag.ch>"

metadata$PointOfContact <- "Nele Schuwirth <nele.schuwirth@eawag.ch>"

metadata$Creator        <- "CSCF - Centre Suisse de Cartographie  de la Faune, Passage Max.-de-Meuron 6, 2000 Neuchâtel, http://www.cscf.ch"

metadata$AccessRights   <- "Private; Permission by Owner needed"

metadata$spatial        <- paste("{\"type\": \"MultiPoint\", \"coordinates\": [ ",
                                 paste(unique(paste("[",
                                                    data$Longitude,", ",data$Latitude,
                                                    "]",
                                                    sep="")),
                                       collapse=", "),
                                 " ] }",sep="")


# ---- Write processed data and metadata ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(data,paste(dir.proc,file.data.proc,"_",d,".dat",sep=""),sep="\t",
            col.names=TRUE,row.names=FALSE,quote=FALSE,na="")
write.table(data.wide,paste(dir.proc,file.data.proc,"_wide_",d,".dat",sep=""),sep="\t",
            col.names=TRUE,row.names=FALSE,quote=FALSE,na="")
write.table(data.wide.taxhom,paste(dir.proc,file.data.proc,"_wide_occ_homogenizedtaxonomy_",d,".dat",sep=""),sep="\t",
            col.names=TRUE,row.names=FALSE,quote=FALSE,na="")
write("Phylum\tClass\tOrder\tFamily\tGenus\tSpecies",paste(dir.proc,file.taxa.proc,"_",d,".dat",sep=""))
write(taxonomy,paste(dir.proc,file.taxa.proc,"_",d,".dat",sep=""),append=TRUE)
write.list(metadata,paste(dir.proc,file.meta.proc,"_",d,".txt",sep=""))

