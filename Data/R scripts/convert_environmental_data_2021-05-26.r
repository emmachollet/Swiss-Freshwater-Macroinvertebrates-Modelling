## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Construct environmental data ----
##
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- July 30, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#here 
rm(list=ls())  # free workspace etc.
graphics.off()
cat("\14")

# packages and functions  ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }      
if ( !require("ecoval") ) { install.packages("ecoval"); library("ecoval") }
msk.morphol <- msk.morphol.1998.create(language="EnglishNodes")
# plot(msk.morphol)
# temp <- ecoval.dictionaries.default  # to look at the dictionary

# directory and file definitions ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.env.data      <- "../Processed data/Environmental data/"
dir.inv.data      <- "../Processed data/Invertebrate data/"
dir.GIS.data      <- "../Original data/"

# Compiled dataset from rs (FRI and bFRI factors integrated with Bogdan Caradima's python scripts)
file.environment        <- "SitesData_for_RS_2020-06-25_result.dat"  

# Data extracted by Bogdan
# now integrated in GIS file from rs
# file.Bogdan             <-  "site_variablesBogdan.csv" 

# Substrate data extracted by Nele (with thanks to Karin for entering the data)
# extracted with "/R scripts/convert_habitatdata_BDM_2020-08-11.r"
file.sub.bdm               <- "BDM_Habitat_substrate_coverage_2020-08-11.dat" # extracted by nis
file.sub.bdm.fract         <- "BDM_Habitat_samplesubstratefraction_2020-08-11.dat" # extracted by nis

# Ecomorphology for BDM from Hintermann and Weber (Nicolas Martinez, Christian Stickelberger)
# extracted with "/R scripts/convert_ecomorph_BDM_2020-08-11.r"
file.InS                 <- "BDM_EcoMorph_data_HintermannWeber_2020-08-11.dat" #xxx to be updated for new sites, contact Nicolas.

# invertebrate data
file.inv.data           <- "All_occ_data_2020-06-25.dat"
file.inv.BDM.data       <- "BDM_occ_data_2020-06-25.dat"

# Only colnames of env fact sorted by importance (need to be updated if there is a change in workflow)
file.rank.env     <- "ranking_env_data.csv"


# read data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Rosi`s compiled dataset
data.environment        <- read.delim(paste(dir.GIS.data,file.environment,sep=""), na = c("<Null>", " ", "NA"), # directly replace problems by Na
                                      header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Substrate data extracted by Nele
data.sub             <- read.delim(paste(dir.env.data,file.sub.bdm,sep=""),header=TRUE,sep="\t", stringsAsFactors=FALSE)
data.sub.fract       <- read.delim(paste(dir.env.data,file.sub.bdm.fract,sep=""),header=TRUE,sep="\t", stringsAsFactors=FALSE)

# Ecomorphology for BDM from Hintermann and Weber (Nicolas Martinez)
data.InS                 <- read.delim(paste(dir.env.data,file.InS,sep=""),header=T,sep="\t", stringsAsFactors=FALSE)  

# invertebrate data
data.inv                <- read.delim(paste(dir.inv.data,file.inv.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE) 
data.inv.BDM            <- read.delim(paste(dir.inv.data,file.inv.BDM.data,sep=""),header=T,sep="\t", stringsAsFactors=FALSE) 

dim(data.inv) # 3081  150

# Only colnames of env fact sorted by importance
rank.env          <- read.csv(paste(dir.env.data,file.rank.env,sep=""),header=TRUE, sep=";", stringsAsFactors=FALSE)

#check if all sites are in the environmental data set:
sites.inv <- unique(data.inv$SiteId)
sites.env <- unique(data.environment$SiteId)
sites.missing <- setdiff(sites.inv,sites.env) #13 sites missing (probably problems with coordinates so excluded by Rosi)

# ind <- match(sites.missing,data.inv[,"SiteId"])
# data.inv[ind,"MonitoringProgram"]

# data preparation ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## invertebrate dataset

cat("number of sites:", length(unique(data.inv$SiteId)),"\n")
cat("number of samples:", length(data.inv$SampId),"\n")

## environmental dataset

dim(data.environment) # 3081  174
length(unique(data.environment$SampId))
length(unique(data.environment$SiteId))


c.ind <- grep("Occurrence.",colnames(data.environment)) 

if (length(c.ind) > 0){
  data.environment <- data.environment[,-c.ind]  #xxx ecr 8.6.21 remove columns with occ taxa because of misunderstanding with rs,
                                                # but for next time don't need it, it's corrected in the workflow and GIS file is already corrected
                                                # because of this there are also useless columns like "OBJECTID...1" but shouldn't be a problem
}
  
data.environment <- data.environment[,-which(colnames(data.environment)=="SiteId")] #xxx nis 7.4.21 remove column to avoid duplicated columns when joining by SampID

## combine env and inv, in order to add a SampleId to env and only retain those samples that have invertebrate observations

data <- data.inv[,c("SiteId", "SampId")]
data <- left_join(data, data.environment, by=c("SampId"))   #xxx nis: changed 7.4.21
dim(data)

# ind.dif <- which(is.na(data$MonitoringProgram.x == data$MonitoringProgram.y)) # check if some columns duplicated, should be empty

## <  add substrate data ##############################################################################

substrata <- c("mobile_blocks","coarse_inorganic_sediments", "gravel",
               "sand_silt", "fine_sediments", "algae_clay", 
               "natural_artifical_surfaces", "moss", "hydrophytes",  "helophytes", "coarse_organic_matter" )
normcov.substrata <- paste("normcov", substrata, sep="_")
covclass.substrata <- paste("covclass", substrata, sep="_")
colnames(data.sub.fract)[4:14] <- paste("sfract", colnames(data.sub.fract)[4:14], sep="_")
velocities <- c("v_5", "v_5_25", "v_25_75", "v_75_150", "v_150")

# Files data.sub and data.sub.fract only have the year (not the exact date),
# so we first create SiteId and then a truncated SampId with only the year
# to finally join it to data

colnames(data.sub)[which(names(data.sub) == "Jahr")] <- "Year"
data.sub <- cbind(SiteId=paste("CSCF_",data.sub$aIdstao,sep=""), data.sub)
data.sub <- cbind(TruncSampId=paste(data.sub$SiteId, data.sub$Year, sep="_"), data.sub)
# data.sub$SampId <- paste(data.sub$SiteId, data.sub$Date, sep="_")
# data.sub.fract$SampId <- paste(data.sub.fract$SiteId, data.sub.fract$Date, sep="_")

colnames(data.sub.fract)[which(names(data.sub.fract) == "Jahr")] <- "Year"
data.sub.fract <- cbind(SiteId=paste("CSCF_",data.sub.fract$aIdstao,sep=""), data.sub.fract)
data.sub.fract <- cbind(TruncSampId=paste(data.sub.fract$SiteId, data.sub.fract$Year, sep="_"), data.sub.fract)

data <- cbind(TruncSampId=paste(data$SiteId, data$Year, sep="_"), data)

data.substrate <- left_join(data.sub[,c("TruncSampId",covclass.substrata,normcov.substrata)],
                            data.sub.fract[,c("TruncSampId",paste0("sfract_",substrata),velocities)],by="TruncSampId")

data <- left_join(data,data.substrate,by="TruncSampId")

data <- data[,-which(colnames(data)=="TruncSampId")] # remove column "TruncSampId" created just for this step

## < add ecomorphology ##############################################################################

# rename columns in the data with morphology so that they start with "CH_"
morph.names <- c("OEKOMKLASSE", "GSBREITE", "EINDOL", "VNATABST",                           
                 "BREITENVAR", "TIEFENVAR", "SOHLVER",
                 "SOHLMAT", "LBUKVER", "RBUKVER", 
                 "LBUKMAT", "RBUKMAT", "LUFBEBRE",          
                 "RUFBEBRE", "LUFERBER", "RUFERBER",
                 "LUFBEBEW","RUFBEBEW", "BEWALGEN",
                 "BEWMAKRO", "TOTHOLZ")
ind <- which(colnames(data) %in% morph.names)
colnames(data)[ind] <- paste("CH", colnames(data)[ind], sep="_")

# add morphology collected during the invertebrate sampling, those columns start with "InS_"
data.InS <- data.InS[,-c(1:6, 8:10)]
colnames(data.InS)[which(colnames(data.InS) == "InS_SampId")] <- "SampId"
data <- left_join(data,data.InS, by="SampId")

## final cleaning and checks

cat("no of unique sites:",length(unique(data$SiteId)),"\n") # 2386 sites
dim(data)   # 3081  290

# note: ideally we should have environmental data for all monitoring sites, but some sites were sampled outside of Switzerland


# calculate environmental factors ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## < mean max temp in summer, maximal morning water temperature in summer ##############

data$A_EZG <- as.numeric(data$A_EZG)   # in m2
data$DEM_EZG <- as.numeric(data$DEM_EZG)  # in m

# data$temp.max.sum <- 17.836 + 1.858 * log10(data$A_EZG/1e6) -  0.005 * data$DEM_EZG  # mean maximum temperature in summer
data$temperature <- 18.4 + 2.77 * log10(data$A_EZG/1e6) -  0.006 * data$DEM_EZG        # mean morning temperature in summer
data[which(data$temperature<0), "temperature"] <- 0  # there are no negative temperatures now, constrained to 0
# cat("missing temperature values for MP:",data[which(is.na(data$temperature)),"MonitoringProgram"],"\n")


## < discharge ##########################################################################

# data$AS_M_S <- as.numeric(data$AS_M_S)
# summary(data$AS_M_S)
# plot(data$AS_M_S, data$Q, xlim=c(0,350),ylim=c(0,350))
data$Discharge <- as.numeric(data$Discharge)


# cat("missing Discharge values for MP:",unique(data[which(is.na(data$Discharge)),"MonitoringProgram"]),"\n")

## < wastewater fraction ################################################################


data$ARA_fraction <- as.numeric(data$ARA)/(data$Discharge*31536000)  

# !! discute with rosi about sites with really small discharge ####
data$ARA_fraction[which(data$ARA_fraction > 200)] <- NA


## < pesticides (agriculture) ###########################################################
  
# Pesticides derived from pesticide landuse (weighted landuse data) 

# preparation                                                                                                 
data$RAPS_ANT <- as.numeric(data$RAPS_ANT)/100   
data$GEM_ANT  <- as.numeric(data$GEM_ANT)/100   
data$OBST_ANT <- as.numeric(data$OBST_ANT)/100    
data$REB_ANT  <- as.numeric(data$REB_ANT)/100    
data$KART_ANT <- as.numeric(data$KART_ANT)/100    
data$GETR_ANT <- as.numeric(data$GETR_ANT)/100   
data$HUEF_ANT <- as.numeric(data$HUEF_ANT)/100    
data$RUEB_ANT <- as.numeric(data$RUEB_ANT)/100    
data$MAIS_ANT <- as.numeric(data$MAIS_ANT)/100    

data$IAR <- 1.83*data$RAPS_ANT +
  2.66*data$GEM_ANT +
  3.10*data$OBST_ANT +
  0.37*data$REB_ANT +
  0.44*data$KART_ANT +
  0.03*data$GETR_ANT +
  0.38*data$HUEF_ANT +
  0.07*data$RUEB_ANT + 
  0.01*data$MAIS_ANT

summary(data$IAR)


## < urban areas [range:0-1] ###########################################################
data$urban_area <- as.numeric(data$SIED_ANT)/100


## < saprobic condition (agriculture + wastewater) #####################################

# preparation
data$cow_density <- as.numeric(data$GVE)/data$A_EZG*1000000  # GVE: cow units, A_EZG: m2
summary(data$cow_density)

data$agri_land <-  as.numeric(data$ACK_ANT) + as.numeric(data$OBST_ANT) + as.numeric(data$REB_ANT)
summary(data$agri_land)

data$fields <- as.numeric(data$NAWHW_ANT) + as.numeric(data$AJWW_ANT)


## 4 models with  interpolation points incl. xeno and using max or mean of the 5 chemical variables =2-dis, DOC, NH4, NO3, PO4
## we keep the model using the mean and including agri land, WW fract and lifestock densities
# data$saproby.max_xeno = 1.287 + 0.0248*data$agri_land + 5.957*data$ARA_fract
# summary(data$saproby.max_xeno)
# data[which(data$saproby.max_xeno>4), "saproby.max_xeno"] <- 4
# data$saproby.max_xeno2 = 1.151 + 0.0220*data$agri_land + 5.716*data$ARA_fract + 0.00727*data$cow_density 
# summary(data$saproby.max_xeno2)
# data[which(data$saproby.max_xeno2>4), "saproby.max_xeno2"] <- 4
# data$saproby.mean_xeno = 0.870 + 0.0208*data$agri_land + 4.648*data$ARA_fract
# summary(data$saproby.mean_xeno)
# data[which(data$saproby.mean_xeno>4), "saproby.mean_xeno"] <- 4
# data$saproby.mean_xeno2 = 0.746 + 0.0182*data$agri_land + 4.427*data$ARA_fract + 0.00668*data$cow_density
# summary(data$saproby.mean_xeno2)
# data[which(data$saproby.mean_xeno2>4), "saproby.mean_xeno2"] <- 4


data$saprobic_cond = 0.746 + 0.0182*data$agri_land + 4.427*data$ARA_fract + 0.00668*data$cow_density
summary(data$saprobic_cond)

data[which(data$saprobic_cond>4), "saprobic_cond"] <- 4


# < heavy metals ###########################################################################################

# consider once we get data from Zurich
# vineyards for Cu (also orchards and organic farming)
# urban and traffic areas


# < organic pollutants (persistant) #################################################################################

# could consider this at a later stage, but
# from experience this does not work too well 


# < river morphology ###########################################################################################

# identify Samples from the BDM
BDM.sites <- data.inv[ which(data.inv$MonitoringProgram =="BDM") , "SampId"]
length(BDM.sites) # should have 886 samples
BDMind <- which(data$SampId %in% BDM.sites)

ind <- !is.na(data[BDMind, "InS_EINDOL"]) # remove new sites for avoiding problems on 28.04.2021
BDMind <- BDMind[ind]



# make new columns for those variables that are needed for the ecoval package to calculate morphology variables

# If available we use data from BDM sampling (after checking there are no data missing for all BDM sites)
# If not, we use CH_ data provided by rs

## EINDOL
data$EINDOL <- as.numeric(data$CH_EINDOL)
length(which(is.na(data[BDMind, "InS_EINDOL"])))      # none are missing from bdm sites
data[BDMind,"EINDOL"] <- data[BDMind,"InS_EINDOL"] 
summary(data$EINDOL)                                  # should be 0 or 1

# GSBREITE
data$GSBREITE <- as.numeric(data$Channel_Width)
length(which(is.na(data[BDMind, "InS_GSBREITE"])))        # none are missing from bdm sites
data[BDMind,"GSBREITE"] <- data[BDMind,"InS_GSBREITE"]
summary(data$GSBREITE)

## BREITENVAR
data$BREITENVAR <- as.numeric(data$CH_BREITENVAR)
length(which(is.na(data[BDMind, "InS_BREITENVAR"])))   # none are missing from bdm sites
data[BDMind,"BREITENVAR"] <- data[BDMind,"InS_BREITENVAR"] 
summary(data$BREITENVAR)                               # should be 1-3
data[which(data$BREITENVAR>3),"BREITENVAR"] <- 3       # We checked this site, it is a straight channel => breitenvar = 3
data[which(data$BREITENVAR<1),"BREITENVAR"] <- 3       # this is ok, we assume that 0 means 3 (no width variability)

# SOHLMAT
data$SOHLMAT <- as.numeric(data$CH_SOHLMAT)
length(which(is.na(data[BDMind, "InS_SOHLMAT"])))      # this is ok, if SOHLVER = 1 (= no sohlver), then there is no SOHLMAT (=NA)
data[BDMind[which(is.na(data[BDMind,"InS_SOHLMAT"]))],"InS_SOHLMAT"] <- 0   # 0 means no sohlmaterial because no sohlverbauung
data[BDMind,"SOHLMAT"] <- data[BDMind,"InS_SOHLMAT"] 
summary(data$SOHLMAT)                                  # should be 0-5
data[which(data$SOHLMAT>5),"SOHLMAT"] <- 5             # we have to check the NAWA morphology, this would likely dissapear then (indeed, it does, see paper 2)

# SOHLVER
data$SOHLVER <- as.numeric(data$CH_SOHLVER)
length(which(is.na(data[BDMind, "InS_SOHLVER"])))      # none are missing from bdm sites
data[BDMind,"SOHLVER"] <- data[BDMind,"InS_SOHLVER"] 
summary(data$SOHLVER)                                  # should be 1-6
data[which(data$SOHLVER<1),"SOHLVER"] <- 1             # this is ok, we assume that 0 means no sohlverbauung, which would be coded as 1

# LBUKMAT
data$LBUKMAT <- as.numeric(data$CH_LBUKMAT)
length(which(is.na(data[BDMind, "InS_LBUKMAT"])))      # many missing, but these have no LBUKVER (LBUKVER = 1)
# this is correct,  0 is coded as impermeable, this makes it clear that something is wrong there. However, as there is no LBUKVER this is not taken further in the calculations. So, no problem
data[BDMind[which(is.na(data[BDMind,"InS_LBUKMAT"])&!is.na(data[BDMind,"BREITENVAR"]))],"InS_LBUKMAT"] <- 0  # 0 means no LBUKVER because no LBUKVER #xxx
# but only for sites with morph info (i.e e.g. with BREITENVAR info)
data[BDMind,"LBUKMAT"] <- data[BDMind,"InS_LBUKMAT"]

# LBUKVER
data$LBUKVER <- as.numeric(data$CH_LBUKVER)
length(which(is.na(data[BDMind, "InS_LBUKVER"])))     # none are missing from bdm sites
data[BDMind,"LBUKVER"] <- data[BDMind,"InS_LBUKVER"] 
data[which(data$LBUKVER<1),"LBUKVER"]  <- 1           # 0 and NA are equal to 1 (class indicating no Lbukver, based on ecoval.dictionaries.default)
data[which(is.na(data$LBUKVER)&!is.na(data$EINDOL)),"LBUKVER"]  <- 1  #xxx nis 09.08.2021
# but only for sites with morph info (i.e e.g. with EINDOL info)
summary(data$LBUKVER)                                 # should be 1-6

# RBUKMAT
data$RBUKMAT <- as.numeric(data$CH_RBUKMAT)
length(which(is.na(data[BDMind, "InS_RBUKMAT"])))     # many missing, but these have no RBUKVER (RBUKVER = 1)
data[BDMind[which(is.na(data[BDMind,"InS_RBUKMAT"])&!is.na(data[BDMind,"BREITENVAR"]))],"InS_RBUKMAT"] <- 0  # 0 means no LBUKVER because no LBUKVER xxx
# but only for sites with morph info (i.e e.g. with BREITENVAR info)
data[BDMind,"RBUKMAT"] <- data[BDMind,"InS_RBUKMAT"]

# RBUKVER
data$RBUKVER <- as.numeric(data$CH_RBUKVER)
length(which(is.na(data[BDMind, "InS_RBUKVER"])))     # none are missing from bdm sites
data[BDMind,"RBUKVER"] <- data[BDMind,"InS_RBUKVER"] 
data[which(data$RBUKVER<1),"RBUKVER"]  <- 1           # 0 and NA are equal to 1 (class indicating no Lbukver, based on ecoval.dictionaries.default)
data[which(is.na(data$RBUKVER)&!is.na(data$EINDOL)),"RBUKVER"]  <- 1  #xxx nis 09.08.2021
# but only for sites with morph info (i.e e.g. with EINDOL info)
summary(data$RBUKVER)                                 # should be 1-6

# LUFBEBRE
data$LUFBEBRE <- as.numeric(data$CH_LUFBEBRE)
length(which(is.na(data[BDMind, "InS_LUFBEBRE"])))   # none are missing from bdm sites
data[BDMind,"LUFBEBRE"] <- data[BDMind,"InS_LUFBEBRE"] 
summary(data$LUFBEBRE)

# RUFBEBRE
data$RUFBEBRE <- as.numeric(data$CH_RUFBEBRE)   
length(which(is.na(data[BDMind, "InS_RUFBEBRE"])))   # none are missing from bdm sites
data[BDMind,"RUFBEBRE"] <- data[BDMind,"InS_RUFBEBRE"] 
summary(data$RUFBEBRE)

# LUFBEBEW
data$LUFBEBEW <- as.numeric(data$CH_LUFBEBEW)
length(which(is.na(data[BDMind, "InS_LUFBEBEW"])))   # none are missing from bdm sites
length(which(data[BDMind, "InS_LUFBEBEW"]<1))        # none are 0
data[BDMind,"LUFBEBEW"] <- data[BDMind,"InS_LUFBEBEW"] 
summary(data$LUFBEBEW)                                # should be 1-3

# RUFBEBEW
data$RUFBEBEW <- as.numeric(data$CH_RUFBEBEW)
length(which(is.na(data[BDMind, "InS_RUFBEBEW"])))    # none are missing from bdm sites
length(which(data[BDMind, "InS_RUFBEBEW"]<1))         # none are 0
data[BDMind,"RUFBEBEW"] <- data[BDMind,"InS_RUFBEBEW"] 
summary(data$RUFBEBEW)                                # should be 1-3

# looking at the LUFBEBEW and RUFBEBEW problem
ind <- unique(which(data$LUFBEBEW==0 | data$RUFBEBEW==0 ))
temp2 <- data[ind, c("SiteId", "SampId", "MonitoringProgram", "Latitude", "Longitude","X", "Y",
            "GSBREITE", "EINDOL", "LBUKVER", "RBUKVER", "LBUKMAT" ,
            "RBUKMAT", "LUFBEBEW", "RUFBEBEW" , "LUFBEBRE", "RUFBEBRE")]
# #### We assumed that 0 means 3 (artificial), we agreed that this is the best option 
data$LUFBEBEW <- replace(data$LUFBEBEW, list= c(which(data$LUFBEBEW==0)), 3)  
data$RUFBEBEW <- replace(data$RUFBEBEW, list= c(which(data$RUFBEBEW==0)), 3)  

# BEWMAKRO
data$BEWMAKRO <- as.numeric(data$CH_BEWMAKRO)
length(which(is.na(data[BDMind, "InS_BEWMAKRO"])))        # none are missing from bdm sites
data[BDMind,"BEWMAKRO"] <- data[BDMind,"InS_BEWMAKRO"]

data[which(is.na(data$BEWMAKRO) & !is.na(data$BREITENVAR)),"BEWMAKRO"]  <- 0        # ?? NA means 0, is it right also for CH_BEWMAKRO ? Or only fo InS ?
data[which(is.na(data$InS_BEWMAKRO)),"InS_BEWMAKRO"]  <- 0 
summary(data$BEWMAKRO)

# calculate and attach values morphology msk
val <- evaluate(msk.morphol,attrib=data)   # errors related to Sohlenverbgrad_Prozent are ok, this is not needed
# plot(msk.morphol,u=val[1:10,],main=data[1:10,"SiteId"])
dim(data)
dim(val)            

colnames(val)[which(colnames(val)=="Good morphological state")] <- "ecomorphology"
colnames(val)[which(colnames(val)=="High width variability")] <- "width variability"
colnames(val)[which(colnames(val)=="No bed modification")] <- "bed modification"
colnames(val)[which(colnames(val)=="No bank modification")] <- "bank modification"
colnames(val)[which(colnames(val)=="Near natural riparian zone")] <- "riparian zone"



data <- cbind(data,val)
dim(data)                    


# little checks
length(which(is.na(data[BDMind, "ecomorphology"])))       # none are missing from bdm sites 
length(which(is.na(data[BDMind, "width variability"])))   # none are missing from bdm sites
length(which(is.na(data[BDMind, "bed modification"])))    # none are missing from bdm sites
length(which(is.na(data[BDMind, "bank modification"])))   # none are missing from bdm sites
length(which(is.na(data[BDMind, "riparian zone"])))       # none are missing from bdm sites

## < flow velocity ###########################################################################################
 
# Note: simple velocity model using S0, n, Q, width
# assumption: n= 0.08
# lots of missing mqn_Jahr, hence I used Bogdan`s calculations (Q in m3/s)

data$Slope <- as.numeric(data$Slope)

# # calculate velocity
# # velocity(m/s), slope (dimensionless), Q(m3/s), width(m) , n=s/m(1/3)
# data$velocity <- (sqrt(data$Slope/100)/0.08)^(3/5)*(data$Discharge/data$GSBREITE)^(2/5)
# if(length(which(data$velocity>5))>0) data <- data[-which(data$velocity>5), ]   ### remove 1 site, 4 samples with extreme discharge values (is a disconnected river section)
# summary(data$velocity)

# Mannings n: 
data[,"n.Man"] <- 0.03 +
  0.03*data[,"width variability"] +
  0.03*(data[,"normcov_gravel"] +
          data[,"normcov_coarse_inorganic_sediments"] +
          data[,"normcov_mobile_blocks"]) +
  0.05*((data[,"BEWMAKRO"]-1)/2) 

# # estimate makrophyte cover
# BEWMAKRO_EST <- (1-(data[,"FRI"]/100)) * 
#   ifelse(data[,"Discharge"]<10 & data[,"Slope"]<2,1,0)
# 
# data[,"n.Man2"] <- 0.03 +
#   0.03*data[,"width variability"] +
#   0.03*(data[,"normcov_gravel"] +
#           data[,"normcov_coarse_inorganic_sediments"] +
#           data[,"normcov_mobile_blocks"]) +
#   0.05*BEWMAKRO_EST


# flow velocity:

v.new <- (sqrt(data[,"Slope"]/100)/data[,"n.Man"])^(3/5) *
  (data[,"Discharge"]/data[,"GSBREITE"])^(2/5)
# v.new2 <- (sqrt(data[,"Slope"]/100)/data[,"n.Man2"])^(3/5) *
#   (data[,"Discharge"]/data[,"GSBREITE"])^(2/5)
v.old <- (sqrt(data[,"Slope"]/100)/0.08)^(3/5) *
  (data[,"Discharge"]/data[,"GSBREITE"])^(2/5)
# make 2 new ones where you fill NA values with the old.vel
na.ind <- which(is.na(v.new))
v.new.nafill <- v.new
v.new.nafill[na.ind] <- v.old[na.ind]
# na.ind <-  which(is.na(v.new2))
# v.new2.nafill <- v.new2
# v.new2.nafill[na.ind] <- v.old[na.ind]
# # plot(data$v.old,data$v.new,pch=16,col=rgb(0,0,0,0.2));abline(a=0,b=1)
# # plot(data$v.old,data$v.new2,pch=16,col=rgb(0,0,0,0.2));abline(a=0,b=1)

data$velocity <- v.new.nafill
summary(data$velocity)

## < riparian zone degradation ################################################################################
# note: 0 riparian zone, means bad riparian zone
# 1 for riparian zone, means best status of the riparian zone => 
# rip.zone.deg [0:none - 1 complete degradation] = 1 - rip zone [0: bad state - 1: good state]
data$rip.zone.deg <- 1 - data[,"riparian zone"]

## < Bed degradation ##########################################################################################

# note: 0 for bed modification means bad bed modification,
# 1 for bed modification means best status of the river bed => 
# bed degradataion [0: no - 1 complete degradation] = 1 - bed modification [0: bad state - 1: good state]
data$bed.deg <- 1 - data[,"bed modification"]



# check data$Note for any outstanding issues ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SiteId = identifier of the monitoring program
# SiteID = identifier of Rosi


# same SiteId, but actually a different location (50m apart) => Rosi gave 2 different SiteIDs to them
# data[ which(grepl("Same SiteId", data$Note)) , "SampId"]
# problem.samples <- c("CSCF_S294_AG_2011-03-24", "CSCF_S294_AG_2016-03-10")
# setdiff(data[problem.samples[1], ], data[problem.samples[2], ]) #xxx ecr 9.6.21 to try to understand
# main.var <- c("temperature", "saprobic_cond", "IAR", "rip.zone.deg", "SOHLVER")
# data[which(data$SampId == problem.samples),main.var]
# data <- data[-which(data$SampId %in% problem.samples)[c(2,4)],] # I decided to remove the repeated site, It looks like it really is the same site
# 
# # same SiteID, but 2 different SiteIds  = the same location, but different Id
# data[ which(grepl("Same Location", data$Note)) , "SampId"]
# Rosi.Id <- unique(data[ which(grepl("Same Location", data$Note)) , "Site_ID"]) # should be 8 sites
# data[which(data$Site_ID== Rosi.Id[2]) , "SampId"]


# column names to be removed

to.remove <- c("Site_ID",   "Note",    "rid" ,            "fmeas",            "tmeas" ,           "gwn_ID",           
"gwn_ID_ARA" ,     "OM_FID",           "Line_ID",         "EZG_in",           "EZG_ARASplit",     "EZG_gwn",          
"FGTID",           "Main_Length",      "Quell_ID",        "Total_Length",    
#"regimetyp",        "regimenr",  "abflussvar",      
"OBJECTID_g",       "herkunft_r",      "ID",               "a",                "r_2",              
"std",             "up_len",           "z1",              "z2",               "OM_ABSNR",         "OM_GWLNR",         
"OM_FMEAS",        "OM_TMEAS",         "ERSTERFASSUNG",   "DATENHERKUNFT",    "IDHERKUNFT",       "LETZTEAENDERUNG",
"LETZTERBEARBEITER", "OBJECTID_G_1",   "EZG_NR_1",        "EZGNR_ARA",        "gwn_ID_ARA_1",     "h1",              
"h2",              "OBJECTID_G_12",    "FLOZ_1",          "EZG_NR_12",        "AS_M_S",           "SEEN_ID",          
"See_EZG",         "DEM_Lake",         "Hypeak_Site_Dist","InS_Bemerkungen",  "InS_MSK_Nr",
"InS_OEKOMKLASSE_1", "InS_Feld") 
#"n.Man",

dim(data)
data <- data[,-which(colnames(data) %in% to.remove)]
dim(data)

# column names to be renamed

colnames(data)[which(colnames(data)=="ecomorphology")] <- "morphology" # not needed with the newest ecoval package


# move morphology columns together
CH_ind <- which(grepl("CH_",colnames(data)))
CH_ind <- CH_ind[1:(length(CH_ind)-2)]  # to remove two colums with DACH_

colnames(data)[CH_ind]

InS_ind <-  which(grepl("InS_",colnames(data)))
morph_ind <- which(colnames(data) %in% sub("InS_", "", colnames(data)[InS_ind]))

data <- data[,c( c(1:ncol(data))[-c(CH_ind, InS_ind, morph_ind)],
                  CH_ind,
                  InS_ind, 
                  morph_ind)
             ]
## ---- Remove columns or rows with too many NA ----

# Lots of Null values, that are missing values (should be NA)
# as.numeric(as.character(data.env$mqn.Jahr))


## ---- Construct env dataset with main factors ----

# If the rank.file below doesn't match with the colnames/env factors of the data produced by this workflow
# write env data set now and produce manually a new rank.file with new colnames/env factors sorted in 0, 1, 2 or 3 levels


# d <- "2020-06-25"  # date of the MIDAT invertebrate datafile
# 
# # write All env data set
# filename <- paste(dir.env.data,"All_environmental_data_",d,".dat", sep="")
# write.table(data, filename,sep="\t",row.names=F,col.names=TRUE)


# filename <- paste(dir.env.data,"raw_environmental_data_",d,".dat", sep="")
# write.table(data, filename,sep="\t",row.names=F,col.names=TRUE)


# replace "_" and " " by "." in colnames to be consistent (also with colnames of invertebrate data)
# we could also change upper or lower case letter but hard to be consistent
colnames(rank.env) <- gsub("_", ".", colnames(rank.env))
colnames(rank.env) <- gsub(" ", ".", colnames(rank.env))

colnames(data) <- gsub("_", ".", colnames(data))
colnames(data) <- gsub(" ", ".", colnames(data))


# sort the columns in 4 categories: the sample/site information = 3,
# the environmental factors to, keep in priority = 2, keep to explore = 1, exclude = 0
info    <- colnames(rank.env)[which(rank.env[1,] == 3)]
prio    <- colnames(rank.env)[which(rank.env[1,] == 2)]
explo   <- colnames(rank.env)[which(rank.env[1,] == 1)]
excl    <- colnames(rank.env)[which(rank.env[1,] == 0)]

# check if colnames of env dataset are the same as the one in ranking file
if ( length(setdiff(colnames(data), colnames(rank.env))) != 0){
  
  print(paste("Colnames of new env data don't match with colnames of", file.rank.env, "Need to update rank.file."), sep="")
  
  data.info <- data[,which(colnames(data) %in% info)]
  data <- data[,-which(colnames(data) %in% excl)]
  
} else {
  
  data.info <- data[, info]
  data <- data[, c(info, prio, explo)]
  
}

# # transform character (that are not info) columns in factors
# temp <- data
# temp[sapply(temp, is.character)] <- lapply(temp[sapply(temp, is.character)], as.factor)
# # temp[,which(colnames(temp) %in% info)] <- as.character(data[,which(colnames(data) %in% info)])
# data <- temp
# inutile de le faire maintenant, chaque fois qu'on relira les fichiers les facteurs seront des characters


# save the environmental data and list of different env fact ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- "2020-06-25"  # date of the MIDAT invertebrate datafile

# write All env data set
filename <- paste(dir.env.data,"All_environmental_data_",d,".dat", sep="")
write.table(data, filename,sep="\t",row.names=F,col.names=TRUE)

# write BDM env data set

## combine env and BDM inv, in order to add a SampleId to env and only retain those samples that have invertebrate observations
temp.data.env <- data[,-which(colnames(data)=="SiteId")] # remove column to avoid duplicated columns when joining by SampID
data.BDM <- data.inv.BDM[,c("SiteId", "SampId")]
data.BDM <- left_join(data.BDM, temp.data.env, by=c("SampId"))   #xxx nis: changed 7.4.21
dim(data.BDM)

filename <- paste(dir.env.data,"BDM_environmental_data_",d,".dat", sep="")
write.table(data.BDM, filename,sep="\t",row.names=F,col.names=TRUE)


# # write ordered list of env fact in txt file
# cat("Information about sample/site:\n", info, "\n", 
#     "Environmental factors to prioritize: \n", prio, "\n",
#     "Environmental factors to explore: \n", explo, "\n",
#     "Environmental factors to exclude: \n", excl, "\n",
#     file = paste(dir.env.data,"list_ordered_env_factors_",d,".txt", sep=""))

