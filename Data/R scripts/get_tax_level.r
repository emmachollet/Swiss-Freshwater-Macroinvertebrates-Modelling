# ======================================================
#
#           Get taxonomic level
#
# ======================================================


get_tax_level <- function(data.inv=data.inv, program=c("NAWA","BDM"), data.tax=data.tax){
  # required input: 
  # data.inv = the wide-format monitoring data of the invertebrates, with taxa in columnnames starting with Occurrence_
  # program = name(s) of the program(s) that you want to analyse, as they appear in the data.inv within the coluöm named "MonitoringProgram"
  # data.tax = a taxonomic dictionary with all taxa in Switzerland and their taxonomic resolution
  
  
  # make a vector of taxa in Switzerland and their respective taxonomic level
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
  
  # get a vector with taxa names in the whole data
  occ.taxa.names <- colnames(data.inv[,grepl("Occurrence", colnames(data.inv))])
  taxa.names <- substring(occ.taxa.names,nchar("Occurrence.")+1)
  taxa.names <- sub("group.", "",taxa.names)
  
  # name the taxa.names, based on the taxonomic resolution of that taxon
  ind <- match(taxa.names, data.tax$taxon)
  names(taxa.names) <- data.tax$lowest.level[ind]
  
  # a matrix to place the taxonomic resolutions in that I find with the loop below
  tax.resolutions <- c( "phylum", "class", "order", "family", "genus", "species")
  resolution <- matrix(NA,nrow = length(tax.resolutions),ncol=length(program))
  rownames(resolution) <- tax.resolutions
  colnames(resolution) <- program
  
  for (prog in program){
    if (prog=="ALL") 
    {
      program.data <- data.inv
    } else {
      # subset the data for the specific program of interest
      program.data <- subset(data.inv, data.inv[,"MonitoringProgram"] %in% prog)
    }
    
    # which taxa are actually present in each program
    taxa.names.present <- names(which(colSums(program.data[,occ.taxa.names], na.rm=TRUE)>0))
    taxa.names.present <- substring(taxa.names.present,nchar("Occurrence.")+1)
    taxa.names.present <- sub("group.", "", taxa.names.present)
    ind.present <- match(taxa.names.present, taxa.names)
    names(taxa.names.present) <- names(taxa.names)[ind.present]
    
    # add the number of taxa in each taxonomic resolution to the `resolution` matrix that I made before
    resolution["phylum",prog] <-    length(which(names(taxa.names.present)== "phylum"))
    resolution["class",prog] <-    length(which(names(taxa.names.present)== "class"))
    resolution["order",prog] <-    length(which(names(taxa.names.present)== "order"))
    resolution["family",prog] <-    length(which(names(taxa.names.present)== "family"))
    resolution["genus",prog] <-    length(which(names(taxa.names.present)== "genus"))
    resolution["species",prog] <-    length(which(names(taxa.names.present)== "species"))
  }
  print(resolution)

}







