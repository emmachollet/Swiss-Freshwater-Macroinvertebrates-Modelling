## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Conversion of Ecomorphology data collected during ----
##                BDM invertebrate monitoring 
##
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 29, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Based on Conversion_BDM_ecomorphology.R 06 11 2017 Paper 1 Peter Vermeiren

rm(list=ls())  # free workspace etc.
graphics.off()

# only for now, delete when file is ready
# setwd("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Data/R scripts")

# ===============================
# directory and file definitions:
# ===============================

dir.data.ecom                <- "../Original data/"
file.BDM.ecom                <- "BDM_EPT_Kopfdaten_vonChristianStickelberger_20200811.csv"                        

dir.data.inv                 <- "../Processed data/Invertebrate data/"
file.inv                     <- "invertebrates_wide_2020-06-25.dat"

dir.GIS.data                 <- "../Original data/"
file.environment.data        <- "GIS_20210526.csv" 

dir.output                   <- "../Processed data/Environmental data/"

# ============================
# read data:
# ============================

data.BDM.ecom                <- read.csv(paste(dir.data.ecom,file.BDM.ecom,sep=""),header=TRUE,sep=";",stringsAsFactors=FALSE)
data.inv                     <- read.delim(paste(dir.data.inv,file.inv,sep=""),header=T,sep="\t",stringsAsFactors=FALSE)
data.environment             <- read.csv(paste(dir.GIS.data,file.environment.data,sep=""),header=TRUE, sep=";", stringsAsFactors=FALSE)
# ============================
# Preparation:
# ============================

# ecomorphology
colnames(data.BDM.ecom)[which(names(data.BDM.ecom) == "Jahr")] <- "Year"

# make a new column "Station.nb." in data.BDM.ecom 
data.BDM.ecom <- cbind(Station.nb.=data.BDM.ecom$aIdstao,data.BDM.ecom) # create a new column with Station nb
# data.BDM.ecom <- cbind(SiteId=paste("CSCF_",data.BDM.ecom$Station.nb.,sep=""),
#                        SampId=paste("CSCF_",data.BDM.ecom$Station.nb.,"_",
#                                     data.BDM.ecom$Year,sep=""),
#                        data.BDM.ecom)
# lines needed for analysis of missing sites

# invertebrates
data.inv$X <- as.character(data.inv$X)
data.inv$Y <- as.character(data.inv$Y)
data.BDM.inv <- subset(data.inv, data.inv[,"MonitoringProgram"]=="BDM")
cind <- match(c("SiteId", "X", "Y", "Year", "SampId", "Station.nb."),colnames(data.BDM.inv) )
data.BDM.inv <- data.BDM.inv[,cind]

station.inv <- unique(data.BDM.inv$Station.nb.)
station.ecom <- unique(data.BDM.ecom$Station.nb.)
station.missing <- setdiff(station.inv,station.ecom)

length(station.inv) # 494
length(station.ecom) # 491

print(station.missing) # 3 stations missing in data.BDM.ecom

for (i in station.missing){
 cat("Row", which(data.BDM.inv[,"Station.nb."] == i), "for", i, "\n")
}

# # looking for ecom observations missing in invertebrates file
# sample.ecom <- data.BDM.ecom$SampId
# sample.inv <- substr(data.BDM.inv$SampId, start = 0, stop = 16)
# sample.missing <- setdiff(sample.ecom,sample.inv)
# 
# length(sample.missing) # 99 truncated SampleID in ecom that are not in inv
# print(sample.missing) # only "CSCF_668179_2017" is named CSCF_668179_bis_2017-06-20 in inv
# # all the others are from 2019 and we don't have inv data from this year
# # print(unique(data.BDM.inv$Year))

# merge the data
data <- merge(data.BDM.inv, data.BDM.ecom, by=c("Year", "Station.nb."), all.x=T)
  # left_join(data.BDM.inv, data.BDM.ecom, by=c("Station.nb.", "Year"))


# 3 site does not have ecomorphology
length(unique(data$SiteId))
length(unique(data$aIdstao))
na.ind <- which(is.na(data$aIdstao))
# na.ind
# data[na.ind, "SiteId"]
# # add data for this site
# data[na.ind, "aIdstao"] <- 686171
# data[na.ind, "aIdbeob"] <- NA
# data[na.ind, "Planungsjahr"] <- 2010
# data[na.ind, "Bemerkungen"] <- ""
# data[na.ind, "Eindolung"] <- "nein"
# data[na.ind, "Breite"] <- 8
# data[na.ind, "Tiefe"] <- 0.2
# data[na.ind, "Abschnittslaenge"] <- 80
# data[na.ind, "Wasserspiegel"] <- "ausgepraegt"
# data[na.ind, "Sohlenverbauung"] <- "keine"
# data[na.ind, "Nat_Abstuerze"] <- "ja"
# data[na.ind, "Tiefevarabilitaet"] <- "ausgepraegt"
# data[na.ind, "Material"] <- ""
# data[na.ind, "Boschung_links"] <- "keine"
# data[na.ind, "Boschung_rechts"] <- "keine"
# data[na.ind, "Durchlaessigkeit_links"] <- "" 
# data[na.ind, "Durchlaessigkeit_rechts"] <- "" 
# data[na.ind, "Uferbreite_links"] <- 16
# data[na.ind, "Uferbreite_rechts"] <- 16
# data[na.ind, "Beschaffenheit_links"] <- "gewaessergerecht"
# data[na.ind, "Beschaffenheit_rechts"] <- "gewaessergerecht"
# data[na.ind, "Absturze"] <- "ja"
# data[na.ind, "Absturztyp"] <- "naturlich"
# data[na.ind, "Absturzmaterial"] <- "naturlich"
# data[na.ind, "Absturzhoehe"] <- 50.0
# data[na.ind, "Bauwerke"] <- "nein"
# data[na.ind, "Baumaterial"] <- ""
# data[na.ind, "Bauhoehe"] <- NA
# data[na.ind, "Praesenz._Trubung"] <- "kein"
# data[na.ind, "Praesenz_Abfaelle"] <- ""
# data[na.ind, "Praesenz_Bewuchs"] <- "kein"
# data[na.ind, "Praesenz_Eisensulfid"] <- "kein"
# data[na.ind, "Praesenz_Feststoffe"] <- "kein"
# data[na.ind, "Praesenz_Geruch"] <- "kein"
# data[na.ind, "Praesenz_Kolmation"] <- "kein"
# data[na.ind, "Praesenz_Schaum"] <- "kein"
# data[na.ind, "Praesenz_Schlamm"] <- "kein"
# data[na.ind, "Praesenz_Verfaerbung"] <- "kein"
# data[na.ind, "Ursache._Trubung"] <- "" 
# data[na.ind, "Ursache_Bewuchs"] <- "" 
# data[na.ind, "Ursache_Eisensulfid"] <- "" 
# data[na.ind, "Ursache_Geruch"] <- "" 
# data[na.ind, "Ursache_Kolmation"] <- "" 
# data[na.ind, "Ursache_Schaum"] <- "" 
# data[na.ind, "Ursache_Schlamm"] <- "" 
# data[na.ind, "Ursache_Verfaerbung"] <- "" 
# data[na.ind, "Bem._Trubung"] <- "" 
# data[na.ind, "Bem_Abfaelle"] <- "" 
# data[na.ind, "Bem_Bewuchs"] <- "" 
# data[na.ind, "Bem_Eisensulfid"] <- "" 
# data[na.ind, "Bem_Geruch"] <- "" 
# data[na.ind, "Bem_Schaum"] <- "" 
# data[na.ind, "Bem_Schlamm"] <- "" 
# data[na.ind, "Bem_Verfaerbung"] <- "" 
# data[na.ind, "Algen"] <- "mittel"
# data[na.ind, "Moose"] <- "< 10%"
# data[na.ind, "Makrophyten"] <- "< 10%"
# data[na.ind, "MSK_Nr"] <- NA
# data[na.ind, "MSK_Klasse"] <- NA
# data[na.ind, "Feld"] <- NA
# 


# # aIdstao or aIdbeob, which unique aIdstao are sampled 2ce, and when
# length(unique(data$aIdstao))
# length(unique(data$aIdbeob))
# nrow(data)
# data$count<- rep(1,nrow(data))
# data.repeat <- aggregate(data$count, by=list(data$aIdstao), FUN="sum")
# colnames(data.repeat)[which(colnames(data.repeat)=="Group.1")] <- "aIdstao"
# colnames(data.repeat)[which(colnames(data.repeat)=="x")] <- "sum"
# ind <- which(data.repeat$sum >1) 
# stations.repeated <- data.repeat[ind,"aIdstao"]
# rind <- which(data$aIdstao %in% as.character(stations.repeated))
# length(rind)
# unique(data[rind,"Jahr"])
# length(which(data[rind,"Jahr"]=="2015"))
# # remove the ones that are repeatedly sampled in 2015
# data <- data[-rind,]

# ============================
# convert column names and contents
# ============================

# EINDOL
colnames(data)[which(colnames(data) == "Eindolung")] <- "EINDOL"
data$EINDOL <- replace(data$EINDOL, data$EINDOL=="nein",0)
data$EINDOL <- replace(data$EINDOL, data$EINDOL=="Ja",1)
which(is.na(data$EINDOL))
table(data$EINDOL)

# GSBREITE
colnames(data)[which(colnames(data) == "Breite")] <- "GSBREITE"
which(is.na(data$GSBREITE))
summary(data$GSBREITE)
## 1 NA

# GSTIEFE
colnames(data)[which(colnames(data) == "Tiefe")] <- "GSTIEFE"
which(is.na(data$GSTIEFE))  
summary(data$GSTIEFE)
## 114 NA

# Abschnittslaenge = length of section under ecomorphological description
summary(data$Abschnittslaenge)
which(is.na(data$Abschnittslaenge)) 
## 2 NA

# BREITENVAR 
colnames(data)[which(colnames(data) == "Wasserspiegel")] <- "BREITENVAR"
## there are 5 categories of Breitenvar in Rosi`s data, one is 0, but there is 1 5, likely a mistake
unique(data.environment$BREITENVAR)
## the one with the 5 is CSCF_AAR010_BE
data.environment[which(data.environment$BREITENVAR==5),"SiteId"]
unique(data$BREITENVAR) 
data$BREITENVAR <- replace(data$BREITENVAR, data$BREITENVAR=="ausgepraegt",1)
data$BREITENVAR <- replace(data$BREITENVAR, data$BREITENVAR=="eingeschraenkt",2)
data$BREITENVAR <- replace(data$BREITENVAR, data$BREITENVAR=="keine",3)
which(is.na(data$BREITENVAR))
which(data.environment$BREITENVAR==5)
table(data$BREITENVAR)

# SOHLVER
colnames(data)[which(colnames(data) == "Sohlenverbauung")] <- "SOHLVER"
unique(data.environment$SOHLVER) 
unique(data$SOHLVER)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="keine",1)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="< 10%",2)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="10-30%",3)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="30-60%",4)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="> 60%",5)
data$SOHLVER <- replace(data$SOHLVER, data$SOHLVER=="100%",6)
which(is.na(data$SOHLVER))
table(data$SOHLVER)

# VNATABST
colnames(data)[which(colnames(data) == "Nat_Abstuerze")] <- "VNATABST"
unique(data.environment$VNATABST)  
unique(data$VNATABST)
data$VNATABST <- replace(data$VNATABST, data$VNATABST=="nein",0)
data$VNATABST <- replace(data$VNATABST, data$VNATABST=="ja",1)
## 1 empty site and 1 NA
which(data$VNATABST=="")                                                                  
which(is.na(data$VNATABST))
table(data$VNATABST)

# TIEFENVAR
colnames(data)[which(colnames(data) == "Tiefevarabilitaet")] <- "TIEFENVAR"
unique(data.environment$TIEFENVAR) 
unique(data$TIEFENVAR)
data$TIEFENVAR <- replace(data$TIEFENVAR, data$TIEFENVAR=="ausgepraegt",1)
# data$TIEFENVAR <- replace(data$TIEFENVAR, data$TIEFENVAR=="maessig",2)
data$TIEFENVAR <- replace(data$TIEFENVAR, data$TIEFENVAR=="eingeschraenkt",2)
data$TIEFENVAR <- replace(data$TIEFENVAR, data$TIEFENVAR=="keine",3)
## 1 NA
which(is.na(data$TIEFENVAR))

# SOHLMAT
colnames(data)[which(colnames(data) == "Material")] <- "SOHLMAT"
unique(data.environment$SOHLMAT)
unique(data$SOHLMAT) 
data$SOHLMAT <- replace(data$SOHLMAT, data$SOHLMAT=="Steine",1)
data$SOHLMAT <- replace(data$SOHLMAT, data$SOHLMAT=="Holz",2)
data$SOHLMAT <- replace(data$SOHLMAT, data$SOHLMAT=="Beton",3)
data$SOHLMAT <- replace(data$SOHLMAT, data$SOHLMAT=="undurchlaessig",4)
data$SOHLMAT <- replace(data$SOHLMAT, data$SOHLMAT=="andere (dicht)",5)
## lots of empties -> interestingly, they all have 1 as sohlver (=no sohlver)
no.sohlmat <- which(data$SOHLMAT=="")
length(no.sohlmat)
data[no.sohlmat,"SOHLVER"]
## 1 NA
which(is.na(data$SOHLMAT))
table(data$SOHLMAT)

# LBUKVER
colnames(data)[which(colnames(data) == "Boschung_links")] <- "LBUKVER"
unique(data.environment$LBUKVER) 
unique(data$LBUKVER)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="keine",1)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="< 10%",2)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="10-30%",3)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="30-60%",4)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="> 60%",5)
data$LBUKVER <- replace(data$LBUKVER, data$LBUKVER=="100%",6)
## 1 NA
which(is.na(data$LBUKVER))
table(data$LBUKVER)

# RBUKVER
colnames(data)[which(colnames(data) == "Boschung_rechts")] <- "RBUKVER"
unique(data.environment$RBUKVER)
unique(data$RBUKVER) 
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="keine",1)
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="< 10%",2)
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="10-30%",3)
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="30-60%",4)
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="> 60%",5)
data$RBUKVER <- replace(data$RBUKVER, data$RBUKVER=="100%",6)
## 1 NA
which(is.na(data$RBUKVER))
table(data$RBUKVER)

# LBUKMAT
colnames(data)[which(colnames(data) == "Durchlaessigkeit_links")] <- "LBUKMAT"
unique(data.environment$LBUKMAT)                                                         
unique(data$LBUKMAT) 
data$LBUKMAT <- replace(data$LBUKMAT, data$LBUKMAT=="durchlaessig",1)  ## note: we only have 2 classes here, we coded them as 1 and 4
data$LBUKMAT <- replace(data$LBUKMAT, data$LBUKMAT=="undurchlaessig",4)
which(is.na(data$LBUKMAT))
table(data$LBUKMAT)

# RBUKMAT
colnames(data)[which(colnames(data) == "Durchlaessigkeit_rechts")] <- "RBUKMAT"
unique(data.environment$RBUKMAT)                                                          
unique(data$RBUKMAT) 
data$RBUKMAT <- replace(data$RBUKMAT, data$RBUKMAT=="durchlaessig",1)  ## note: we only have 2 classes here, we coded them as 1 and 4
data$RBUKMAT <- replace(data$RBUKMAT, data$RBUKMAT=="undurchlaessig",4)
which(is.na(data$RBUKMAT))
table(data$RBUKMAT)

# LUFBEBRE
colnames(data)[which(colnames(data) == "Uferbreite_links")] <- "LUFBEBRE"
summary(data$LUFBEBRE) 
## 1 NA
which(is.na(data$LUFBEBRE))

# RUFBEBRE
colnames(data)[which(colnames(data) == "Uferbreite_rechts")] <- "RUFBEBRE"
summary(data$RUFBEBRE) 
## 1 NA
which(is.na(data$RUFBEBRE))

# LUFBEBEW
colnames(data)[which(colnames(data) == "Beschaffenheit_links")] <- "LUFBEBEW"
unique(data.environment$LUFBEBEW)
unique(data$LUFBEBEW)
data$LUFBEBEW <- replace(data$LUFBEBEW, data$LUFBEBEW=="gewaessergerecht",1)
data$LUFBEBEW <- replace(data$LUFBEBEW, data$LUFBEBEW=="gewaesserfremd",2)
data$LUFBEBEW <- replace(data$LUFBEBEW, data$LUFBEBEW=="kunstlich",3)
## 1 NA
which(is.na(data$LUFBEBEW))
table(data$LUFBEBEW)

# RUFBEBEW
colnames(data)[which(colnames(data) == "Beschaffenheit_rechts")] <- "RUFBEBEW"
unique(data.environment$RUFBEBEW)
unique(data$RUFBEBEW)
data$RUFBEBEW <- replace(data$RUFBEBEW, data$RUFBEBEW=="gewaessergerecht",1)
data$RUFBEBEW <- replace(data$RUFBEBEW, data$RUFBEBEW=="gewaesserfremd",2)
data$RUFBEBEW <- replace(data$RUFBEBEW, data$RUFBEBEW=="kunstlich",3)
## 1NA
which(is.na(data$RUFBEBEW))
table(data$RUFBEBEW)

########################### below are new columns, I kept the German column names, but classified similar to previous variables ################################

# Absturze    
unique(data$Absturze)
data$Absturze <- replace(data$Absturze, data$Absturze=="nein",0)
data$Absturze <- replace(data$Absturze, data$Absturze=="ja",1)
## 1NA
which(is.na(data$Absturze))
table(data$Absturze)

# Absturztyp
unique(data$Absturztyp)
## seems to be the category "unbekannt"
data$Absturztyp <- replace(data$Absturztyp, data$Absturztyp=="",0) 
data$Absturztyp <- replace(data$Absturztyp, data$Absturztyp=="naturlich",1)
data$Absturztyp <- replace(data$Absturztyp, data$Absturztyp=="kunstlich",2)
## 1NA
which(is.na(data$Absturztyp))
table(data$Absturztyp)

# Absturzhoehe
summary(data$Absturzhoehe)
## many NA

# Bauwerke
unique(data$Bauwerke)
data$Bauwerke <- replace(data$Bauwerke, data$Bauwerke=="nein",0)
data$Bauwerke <- replace(data$Bauwerke, data$Bauwerke=="ja",1)
## 1NA
which(is.na(data$Bauwerke))
## many empties, not sure what this means
which(data$Bauwerke=="")  
table(data$Bauwerke=="") 

# Baumaterial
unique(data$Baumaterial)                                                         ############ I did not change these material codes

# Bauhoehe                
unique(data$Bauhoehe) 

# Praesenz._Trubung 
unique(data$Praesenz._Trubung)                                             ############ I could not find a classification for this, I kept as is

# Praesenz_Abfaelle
unique(data$Praesenz_Abfaelle)                                             ############ I could not find a classification for this, I kept as is

# Praesenz_Bewuchs
unique(data$Praesenz_Bewuchs)                                               ############ I could not find a classification for this, I kept as is

# Praesenz_Eisensulfid
unique(data$Praesenz_Eisensulfid)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Feststoffe"  
unique(data$Praesenz_Feststoffe)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Geruch"  
unique(data$Praesenz_Geruch)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Kolmation"
unique(data$Praesenz_Kolmation)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Schaum"  
unique(data$Praesenz_Schaum)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Schlamm"  
unique(data$Praesenz_Schlamm)                                         ############ I could not find a classification for this, I kept as is

# "Praesenz_Verfaerbung"  
unique(data$Praesenz_Verfaerbung)                                         ############ I could not find a classification for this, I kept as is


# Ursache._Trubung"     "Ursache_Bewuchs"      "Ursache_Eisensulfid"  "Ursache_Geruch"       "Ursache_Kolmation"    "Ursache_Schaum"    "Ursache_Schlamm"      "Ursache_Verfaerbung"  "Bem._Trubung"         
# "Bem_Abfaelle"         "Bem_Bewuchs"          "Bem_Eisensulfid" "Bem_Geruch"           "Bem_Schaum"           "Bem_Schlamm"          "Bem_Verfaerbung"



# BEWALGEN
colnames(data)[which(colnames(data) == "Algen")] <- "BEWALGEN"
unique(data$BEWALGEN)
data$BEWALGEN <- replace(data$BEWALGEN, data$BEWALGEN=="< 10%",1)
data$BEWALGEN <- replace(data$BEWALGEN, data$BEWALGEN=="mittel",2)
data$BEWALGEN <- replace(data$BEWALGEN, data$BEWALGEN=="> 50%",3)
## 1 NA
which(is.na(data$BEWALGEN))

# Moose
unique(data$Moose)
data$Moose <- replace(data$Moose, data$Moose=="< 10%",1) 
data$Moose <- replace(data$Moose, data$Moose=="mittel",2)
data$Moose <- replace(data$Moose, data$Moose=="> 50%",3)
## 1 NA
which(is.na(data$Moose))

# BEWMAKRO
colnames(data)[which(colnames(data) == "Makrophyten")] <- "BEWMAKRO"
unique(data$BEWMAKRO)
data$BEWMAKRO <- replace(data$BEWMAKRO, data$BEWMAKRO=="< 10%",1)
data$BEWMAKRO <- replace(data$BEWMAKRO, data$BEWMAKRO=="mittel",2)
data$BEWMAKRO <- replace(data$BEWMAKRO, data$BEWMAKRO=="> 50%",3)
## 1 NA
which(is.na(data$BEWMAKRO))

# MSK_Nr
summary(data$MSK_Nr)

# OEKOMKLASSE_1
colnames(data)[which(colnames(data) == "MSK_Klasse")] <- "OEKOMKLASSE_1"
table(data$OEKOMKLASSE_1)

# Feld
unique(data$Feld)

# TOTHOLZ?-> not part of ecomorph assessment


# ============================
# put EM for ecomorphology in front of all the columnnames
# ============================
colnames(data) <- paste("InS", colnames(data),sep="_")


# ============================
# save the dataset
# ============================


date <- "2020-08-11" # Date at which the BDM Kopfdaten data was extracted in the file used here

write.table(data,paste(dir.output,"BDM_EcoMorph_data_HintermannWeber_",date,".dat", sep=""),sep="\t",row.names=T,col.names=NA)




