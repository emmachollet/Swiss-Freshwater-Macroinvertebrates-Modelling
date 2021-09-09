## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Convert BDM habitat data ----
##
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 29, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
graphics.off()
cat("\14")


# = libraries and global definitions ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# libraries:

if ( !require("readxl") ) { install.packages("readxl"); library(readxl) }
# source("utilities.r")

# files and directories

#inp.file.BDM <- "../Processed data/BDM_Habitat_Data_2010-11-12-13-14-15_KG_NIS_NM_Final.xlsx"
inp.file.BDM <- "../Original data/BDM_EPT_Kopfdaten_vonChristianStickelberger_20200811.csv"

dir.output              <- "../Processed data/Environmental data/"

d <- "2020-08-11" # date when the Kopfdaten dataset was exported, needed to date output files

# global definitions

substrate_types <- 10:0
names(substrate_types) <- c("mobile Blocks > 250 mm",
                            "moss",
                            "hydrophytes",
                            "coarse organic matter",
                            "coarse inorganic sediments 25-250mm",
                            "gravel",
                            "helophytes",
                            "(organic) fine sediments <0.1mm",
                            "sand and silt 0.1-0.25mm",
                            "natural or artifical surfaces",
                            "algae (or clay)")

revised_names_stypes <- c("mobile_blocks",
                          "moss",
                          "hydrophytes",
                          "coarse_organic_matter",
                          "coarse_inorganic_sediments",
                          "gravel",
                          "helophytes",
                          "fine_sediments",
                          "sand_silt",
                          "natural_artifical_surfaces",
                          "algae_clay")

cov_class <- 1:4
names(cov_class) <- c("1-5%",
                      "6-10%",
                      "11-50%",
                      "51-100%")

lo_cov_class <- c(0.01,0.06,0.11,0.51)
up_cov_class <- c(0.05,0.1, 0.5, 1)

mid_cov_class <- (lo_cov_class+up_cov_class)/2

names(lo_cov_class) <- names(cov_class)
names(up_cov_class) <- names(cov_class)
names(mid_cov_class) <- names(cov_class)

convert_covclass <- function(cl,classes,values)
{
  values <- as.numeric(values)
  classes <- as.numeric(classes)
  cl <- as.numeric(cl)
  val=rep(NA,length(cl))
  for(i in 1:length(cl))
  {
    temp <- values[which(classes==cl[i])]
    if(sum(temp)>0) val[i] <- temp
  }
  return(val)
}

translate_names <- function(x, oldnames,newnames)
{
  xnew <- x  
  indo <- which(oldnames==x)
  if(sum(indo)>0) xnew <- newnames[indo] else cat("name",x,"not known\n")
  return(xnew)
}

# convert_covclass(cl=4:2,classes=cov_class,values=lo_cov_class)

v_class <- c(2,4,5,3,1)
names(v_class) <- paste0("v_",c(">150","75-150","25-75","5-25","<5")) # cm/s
revised_names_vclass <-  c("v_150","v_75_150","v_25_75","v_5_25","v_5")
v_class <- sort(v_class)

convert_vclass <- function(v)
{
  v <- as.numeric(v)
  if(is.na(v)) c=NA else {
    if(v<5) c=1 
    if(v>=150) c=2
    if(v>=5 & v<25) c=3
    if(v>=75 & v<150) c=4
    if(v>=25 & v<75) c=5}
  return(c)
}

#convert_vclass(76)
#sapply(c(76,1,50),convert_vclass)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# = convert data BDM ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n\n start converting BDM Data...\n\n ")

# > read data BDM ====
# ~~~~~~~~~~~~~~~~~~~~~

data.BDM <- read.table(inp.file.BDM,header=T, stringsAsFactors = F, sep=";")
# colnames(data.BDM)[15:25] <- paste0("S_",colnames(data.BDM)[15:25])

# > get substrate fractions and flow velocity from the sampled area ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MP = "BDM"
if(MP=="BDM") data <- data.BDM

substrate.frac <- as.data.frame(matrix(nrow=nrow(data),ncol=length(substrate_types)+4))

# colnames(substrate.frac) <- c("SiteId","Day","Month","Year","Date","X","Y",
#                               "Altitude","Watercourse","Locality",names(substrate_types))
# 

colnames(substrate.frac) <- c("aIdstao","aIdbeob","Jahr",
                              names(substrate_types),"error")

substrate.frac[,c("aIdstao","aIdbeob","Jahr")] <- data[,c("aIdstao","aIdbeob","Jahr")]

#formatC(as.numeric(data[,"JOUR"]),digits=NULL,width=2,flag=0)

# substrate.frac$Day         <- formatC(as.numeric(data$Day),width=2,flag=0) # to get leading zeros
# substrate.frac$Month       <- formatC(as.numeric(data$Month),width=2,flag=0)
# substrate.frac$Year        <- data$Year
# substrate.frac$Date        <- as.Date(paste(substrate.frac$Year,
#                                             substrate.frac$Month,
#                                             substrate.frac$Day,sep="-"))
# substrate.frac$X           <- data$X
# substrate.frac$Y           <- data$Y
# substrate.frac$Altitude    <- data$Altitude 
# substrate.frac$Watercourse <- data$Watercourse
# substrate.frac$Locality    <- data$Locality


flow <- as.data.frame(matrix(NA,nrow=nrow(substrate.frac),ncol=6))
colnames(flow) <- c(names(sort(v_class)),"flow error")

for(i in 1:nrow(data))
{
  for(j in 1:8)
  {
    sv <- NA
    if(!is.na(data[i,paste0("Teilprobe_P",j)]) )
    {
      sv <- unlist(strsplit(as.character(data[i,paste0("Teilprobe_P",j)]),split=""))
    }
    if(length(sv)>1)
    {
      s <- which(sv=="s"); v <- which(sv=="v")
      # add substrate types
      ind1 <- s+1; ind2 <- v-1
      if(sum(ind1)>0 & sum(ind2)>0) {  # s and v exist
        if (ind1==ind2) {stype = sv[ind1]} else {stype = paste0(sv[ind1],sv[ind2])} # to account for cases where s is two-digit
      } else {  # if v is missing and only s available
        if(sum(ind1)>0)
        { 
          if(length(sv)==2) {stype = sv[ind1]}
          if(length(sv)==3 & ind1 ==2) {stype=paste0(sv[ind1],sv[ind1+1])}
        }
      }
      stype <- as.numeric(stype)
      if(!is.na(stype))
      {
        ind <- which(substrate_types==stype)
        substrate.frac[i,names(ind)] <- sum(substrate.frac[i,names(ind)],1/8,na.rm=TRUE)
      }
      # add flow velocity class
      
      if(sum(v>0)) {
        vcl <- as.numeric(sv[v+1])
        flow[i,names(v_class)[vcl]] <- sum(flow[i,names(v_class)[vcl]],1/8,na.rm=TRUE)
      } #  else {cat("error for flow velocity in line",i," data missing, comment: ",
      #   unlist(data[i,c("Bemerkung/Pruefen?","geprueft","Bemerkung Nicolas",   
      #  "Entscheidung")]),"\n") }
    }
  }
  rsum <- sum(substrate.frac[i,c(4:14)],na.rm=T)
  
  if(rsum==1) # set other substrate classes to zero
  {
    ind.na <- is.na(substrate.frac[i,c(4:14)])
    substrate.frac[i,c(4:14)][ind.na] <- 0
  }
  
  if(MP=="BDM")
  {
    if(rsum!=1 & rsum>0) {
      cat("error for substrate in line",i,", rowsum=",rsum,"\n")
      substrate.frac[i,"error"]<- paste("substrate area=",rsum)}
    if(rsum==0) cat("substrate data missing in line",i,", rowsum=",rsum,"\n")
  }
  
  rsum2 <- sum(flow[i,1:5],na.rm=T)
  if(rsum2==1) 
  {
    ind.na <- is.na(flow[i,]) 
    flow[i,][ind.na] <- 0
  }
  
  if(MP=="BDM")
  {
    if(rsum2!=1 & rsum2>0) {
      cat("error for flow in line",i,", rowsum=",rsum2,"\n")
      flow[i,"flow error"]<- paste("flow area=",rsum)
    }
    if(rsum2==0) cat("flow data missing in line",i,", rowsum=",rsum2,"\n")
  }
}

colnames(flow)[1:5] <- lapply(colnames(flow)[1:5],translate_names,oldnames=names(v_class),newnames=revised_names_vclass)
colnames(substrate.frac)[4:14] <- lapply(colnames(substrate.frac)[4:14],translate_names,
                                          oldnames=names(substrate_types),newnames=revised_names_stypes)

substrate.frac <- cbind(substrate.frac[,c(1:3,c(20,11,15,16,19,18,14,12,13,17,21)-7,ncol(substrate.frac))],
                        flow[,c(1,3,5,4,2,6)])

substrate.frac$"MonitoringProgram" <- MP

write.table(substrate.frac,paste(dir.output,MP,"_Habitat_samplesubstratefraction_",d,".dat",sep=""),
            sep="\t",col.names = TRUE, row.names = FALSE)


# > collect coverage classes in a matrix ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subs.cov.cl <- matrix(ncol=11,nrow=0)  # coverage classes
colnames(subs.cov.cl) <- revised_names_stypes

dat.sc.cl <- data[,paste0("Substrat_",10:0)]
colnames(dat.sc.cl) <- paste0("covclass_",revised_names_stypes)

dat.sc.mc <- apply(dat.sc.cl,2,convert_covclass,classes=cov_class,values=mid_cov_class) # midvalues of classes

dat.sc.rs <- rowSums(dat.sc.mc,na.rm=T) #rowsums

dat.sc.csn <- dat.sc.mc  # normalized that sum of coverages are 1
colnames(dat.sc.csn) <- paste0("normcov_",revised_names_stypes)

for(i in 1:nrow(dat.sc.mc))
{
  dat.sc.csn[i,] <- dat.sc.mc[i,]/dat.sc.rs[i]
}

for(i in 1:nrow(dat.sc.csn))  # differentiate between 0 and true NA
{
  if(dat.sc.rs[i]>0)
  {
    ind.na <- which(is.na(dat.sc.csn[i,]))
    dat.sc.csn[i,ind.na] <- 0
  }
}

cind.Bem <- match("error",colnames(substrate.frac))
data_substrate_coverage <- cbind(substrate.frac[,c(1:10)],dat.sc.cl,dat.sc.csn,substrate.frac[,cind.Bem])
ind.names <- c((ncol(data_substrate_coverage)+1-length(cind.Bem)):ncol(data_substrate_coverage))
colnames(data_substrate_coverage)[ind.names] <- colnames(substrate.frac)[cind.Bem]


#some checks 
# rowSums(dat.sc.csn)

ind <-  which(is.na(rowSums(dat.sc.csn)))
if(length(ind)>0) 
{
  cat(MP, ": normalized coverage classes not available for ",length(ind)," samples, i:",ind)
  data.problem <- data_substrate_coverage[ind,]
  if(MP=="BDM") data.problem.BDM <- data.problem
}

ind <- which(rowSums(dat.sc.csn)<0.999)
if(length(ind)>0) 
{
  cat(MP,": sum of normalized coverage classe not 1 for ",length(ind)," samples, i:",ind)
  data.problem2 <- data_substrate_coverage[ind,]
  if(MP=="BDM") data.problem2.BDM <- data.problem2
}

write.table(data_substrate_coverage,paste(dir.output,MP,"_Habitat_substrate_coverage_",d,".dat",sep=""),
            sep="\t",col.names = TRUE, row.names = FALSE)

cat(".... ",MP," done\n")

  