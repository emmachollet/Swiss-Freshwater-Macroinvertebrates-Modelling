

# Split data in training and testing datasets for extrapolation
split.data.ml <- function(data, training.ratio, variable = "random", bottom = T){
  
  # add an argument in the fct, eg splitcol, which would be here SiteId
  # first split by this argument, then select rows (aka samples) according to their SiteId
  
  if( variable == "random"){
    sample.size <- floor(ratio*nrow(data)) # size of training set is 80% of the whole dataset
    set.seed(2021) # used for reproducing the same output of simulations
    rind.train  <- sort(sample(nrow(data), size = sample.size))
    data.train  <- data[rind.train,]
    data.test   <- data[-rind.train,]
  } else {
    
    v <- data[,variable]
    q <- rank(v, na.last = F)/length(v) # put the NA in the testing set
    rind.train <- if(bottom == T){ which(q < training.ratio)
    } else { which(q > 1-training.ratio) 
    }
    data.train  <- data[rind.train,]
    data.test   <- data[-rind.train,]
  }
  
  splitted.data <- list("Split1" = list("Training data" = data.train, "Testing data" = data.test))
  
  return(splitted.data)
}

# Split the data in 3 training and testing data sets for CV
split.data <- function(data){
  
      set.seed(2021)  
    
      folds <- groupKFold(data$SiteId, 3)
      
      train1 <- data[folds$Fold1,]
      test1 <- data[-folds$Fold1,]
    
      train2 <- data[folds$Fold2,]
      test2 <- data[-folds$Fold2,]
    
      train3 <- data[folds$Fold3,]
      test3 <- data[-folds$Fold3,]
      
  return(list("Split1" = list("Training data" = train1, "Testing data" = test1), 
              "Split2" = list("Training data" = train2, "Testing data" = test2), 
              "Split3" = list("Training data" = train3, "Testing data" = test3)))
}

# Old splitting function from Jonas
split.data.manual <- function(data,sd){
    repeat{
        inv.data <- data[sample(nrow(data)),]
       
        folds <- groupKFold(inv.data$SiteId, 3) 
        
        #Group by SiteId to not avoid data leakage and then assign randomly to folds 1-3 
        inv.data.fold <- as.data.frame(inv.data %>% group_by(SiteId) %>% summarize(folds = sample(seq(3),1)))
        inv.data <- left_join(inv.data, inv.data.fold, by = "SiteId") 
        
        fold1 <- inv.data[which(inv.data$fold == 1),]
        fold2 <- inv.data[which(inv.data$fold == 2),]
        fold3 <- inv.data[which(inv.data$fold == 3),]
        
        #remove fold information
        fold1 <- subset(fold1, select=-c(folds))
        fold2 <- subset(fold2, select=-c(folds))
        fold3 <- subset(fold3, select=-c(folds))
        
        
        # Combine the folds manually for the joint model
        train1 <- bind_rows(fold1, fold2)
        test1 <- fold3
        message(dim(test1)[1])
        
        train2 <- bind_rows(fold2, fold3)
        test2 <- fold1
        message(dim(test2)[1])
        
        train3 <- bind_rows(fold1, fold3)
        test3 <- fold2
        message(dim(test3)[1])
        
        if(sd(c(dim(fold1)[1],dim(fold2)[1],dim(fold3)[1])) < sd){
            break
        }
    }
    
    
    return(list("Split1" = list("Training data" = train1, "Testing data" = test1), 
                "Split2" = list("Training data" = train2, "Testing data" = test2), 
                "Split3" = list("Training data" = train3, "Testing data" = test3)))
}

# function to center data. If there is only one split i.e. the whole occ data, it needs to be inputed as a list with the data as the only element (i.e. list(inv.occ))
# split <-- list, with one of the three splits done by split.data
# cv <-- boolean, T if we have 3 splits for CV, F if we don't do CV
center.data.old <- function(split, CV){
  
    #split <- list(inv.occ)
    #split <- splits[[1]]
    training.data <- split[[1]]
    
    #extract column names to access predictors and invertebrate data seperately
    inv.names <- colnames(select(training.data, contains("Occurrence.group."), contains("Occurrence.")))
    env.names <- colnames(select(training.data, - all_of(inv.names)))
    
    
    #Here I convert to numeric to count observations across sites
    inv.data <- as.data.frame(apply(training.data[, inv.names],2,as.numeric))
    
    # drop TAXA without observations or only presence at the selected sites:
    ind <- apply(inv.data,2,sum, na.rm = T) <= 0
    inv.data <- inv.data[, !ind]
    n.taxa <- ncol(inv.data)
    
    ind <-  apply(inv.data,2,sum, na.rm = T) == ncol(inv.data)
    inv.data <- inv.data[, !ind]
    n.taxa <- ncol(inv.data)
    
    names.selected <- colnames(inv.data)
    training.data <- cbind(training.data[,env.names],training.data[,names.selected])
 
    #if center is true substract the mean of each predictor, check if its divded by sd, I added the division by sd
    mean.env.cond <- apply(select(training.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")), 2, function(k){
        #mean.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
        mean(k, na.rm = TRUE)
    })
    sd.env.cond <- apply(select(training.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")), 2, function(k){
        sd(k, na.rm = TRUE)
    })
    for(i in 1:length(select(training.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")))){
        #i = 6
        #message(i)
        #Bit ugly good, check that the indices are right (i.e. the ones for the env data)
        training.data[i+4] <- as.matrix(training.data[i+4]) - mean.env.cond[i]
    }
    for(i in 1:length(select(training.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")))){
        #i = 6
        training.data[i+4] <- as.matrix(training.data[i+4]) / sd.env.cond[i]
    }
    
    if(CV == F){
        return(list( "Entire dataset" = training.data))
    }else{
        
        testing.data <- split[[2]]
        
        #Here I make sure that the same species are dropped from the training and testing set.
        testing.data <- cbind(testing.data[,env.names],testing.data[,names.selected])
        
        # join the environmental conditions to the occurrence data
        #test.predictors <- left_join(testing.data[, c("SiteId", "SampId")], env.cond.orig, by = c("SiteId", "SampId")
        
        for(i in 1:length(select(testing.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")))){
            #i = 6
            #message(i)
            #Bit ugly good, check that the indices are right (i.e. the ones for the env data)
            testing.data[i+4] <- as.matrix(testing.data[i+4]) - mean.env.cond[i]
        }
        for(i in 1:length(select(testing.data, all_of(env.names), - c("SiteId", "SampId","X", "Y")))){
            #i = 6
            testing.data[i+4] <- as.matrix(testing.data[i+4]) / sd.env.cond[i]
        }
        
        return(list("Training data" = training.data, "Testing data" = testing.data, "Mean" = mean.env.cond, "SD" = sd.env.cond))
    }
    
}

center.data <- function(data, split, CV, extrapol, dl, mean.dl, sd.dl, env.fact.full){
  
  training.data <- split[[1]]
  
  # extract column names to access predictors and invertebrate data separately
  inv.names <- colnames(select(training.data, contains("Occurrence.")))
  env.names <- env.fact.full
  info.names <- colnames(select(training.data, - all_of(inv.names), - all_of(env.names)))
  
  # convert to numeric to count observations across sites
  inv.data <- as.data.frame(apply(training.data[, inv.names],2,as.numeric))
  
  # WAS ADDED TO PRE-PROCESSING COULD BE REMOVED
  # drop taxa without observations or only presence at the selected sites
  ind <- apply(inv.data,2,sum, na.rm = T) <= 0
  inv.data <- inv.data[, !ind]
  n.taxa <- ncol(inv.data)
  
  ind <-  apply(inv.data,2,sum, na.rm = T) == ncol(inv.data)
  inv.data <- inv.data[, !ind]
  n.taxa <- ncol(inv.data)
  
  names.selected <- colnames(inv.data)
  training.data <- cbind(training.data[, info.names],training.data[,env.names],training.data[,names.selected])
  
  # if center is true substract the mean of each predictor, check if its divided by sd, I added the division by sd
  mean.env.cond <- apply(select(training.data, all_of(env.names)), 2, function(k){
    mean(k, na.rm = TRUE)
  })
  sd.env.cond <- apply(select(training.data, all_of(env.names)), 2, function(k){
    sd(k, na.rm = TRUE)
  })
  
  if(dl == T){
    
    for(env in env.names){
      training.data[env] <- training.data[env] -  mean.dl[env]
    }
    
    for(env in env.names){
      training.data[env] <- training.data[env] / sd.dl[env]
    }
    
  }else{
    
    for(env in env.names){
      training.data[env] <- training.data[env] -  mean.env.cond[env]
    }
    
    for(env in env.names){
      training.data[env] <- training.data[env] / sd.env.cond[env]
    }
  }

  if(CV == F & extrapol == F){
    return(list( "Entire dataset" = training.data))
  }else{
    
    testing.data <- split[[2]] # load testing data 
    
    # Here I make sure that the same species are dropped from the training and testing set.
    testing.data <- cbind(testing.data[,info.names],testing.data[,env.names],testing.data[,names.selected])

    if(dl == T){
      
      for(env in env.names){
        testing.data[env] <- testing.data[env] -  mean.dl[env]
      }
        
      for(env in env.names){
        testing.data[env] <- testing.data[env] / sd.dl[env]
      }
      
    }else{
      
      for(env in env.names){
        testing.data[env] <- testing.data[env] -  mean.env.cond[env]
      }
      
      for(env in env.names){
        testing.data[env] <- testing.data[env] / sd.env.cond[env]
      }
    }
    return(list("Training data" = training.data, "Testing data" = testing.data, "Mean" = mean.env.cond, "SD" = sd.env.cond))
  }
}


preprocess.data <- function(data.env, data.inv, prev.inv, env.fact.full, dir.workspace, BDM, dl, CV, extrapol, extrapol.info = c()){
    
    # Merge data sets
    cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))
    list.taxa <- colnames(data.inv)[cind.taxa]
    data <- data.env[, c("SiteId", "SampId", "X", "Y", env.fact.full)] %>%
        left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
    dim(data)
    
    # Print information
    cat("\nSummary information of original datasets before preprocessing:\n",
        length(unique(data$SampId)), "samples,\n",
        length(unique(data$SiteId)), "sites,\n",
        length(list.taxa), "taxa (including missing values) with taxonomic level:\n")
    
    print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
    
    # Drop taxa with too many NAs
    cind.taxa <- which(grepl("Occurrence.", colnames(data)))
    too.many.na <- c()
    threshold <- ifelse(BDM, 100, 200)
    for(i in cind.taxa){
      if(sum(is.na(data[,i])) > threshold){ too.many.na <- c(too.many.na, i)}
    }
    cat("\nThe following", length(too.many.na), "taxa are excluded because they have more than ", threshold, "NAs :\n", gsub("Occurrence.", "", colnames(data)[too.many.na]), "\n")
    data<- data[, -too.many.na]
    
    
    # Drop rows with incomplete taxa or influence factors
    rind <- !apply(is.na(data), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\n", sum(!rind),"sites/samples excluded because of incomplete taxa or influence factors.\n"))
    data <- data[rind,]
    dim(data)
    
    # Delete taxa with only presence or absence (happens because we deleted rows)
    # Delete taxa with less than 5% of only present or absent
    cind.taxa <- which(grepl("Occurrence.",colnames(data)))
    no.samples <- nrow(data)
    five.perc <- 0.05*no.samples
    cind.rem <- c()
    for (j in cind.taxa) {
        if(sum(data[,j]) < five.perc | sum(data[,j]) > no.samples - five.perc){ cind.rem <- c(cind.rem, j) }
    }
    cat("\nThe following", length(cind.rem), "taxa are excluded because they have less than 5% or more than 95% of prevalence:\n", gsub("Occurrence.", "", colnames(data)[cind.rem]), "\n")
    data <- data[,-cind.rem]
    dim(data)
    
    # Split data
    
    prefix <- ifelse(BDM, "BDM_", "All_")
    
    if(CV == T){
        
        # Split for CV
        file.name <- paste0(dir.workspace, prefix, "SplitsForCV.rds")
        
        # If the file with the three different splits already exist, just read it
        if (file.exists(file.name) == T ){
            
            if(exists("splits") == F){ splits <- readRDS(file = file.name)
            cat("\nFile with data splits already exists, we read it from", file.name, "and save it in object 'splits'.\n")}
            else{
                cat("\nList with data splits already exists as object 'splits' in this environment.\n")
            }
        } else {
            
            cat("\nNo data splits exist yet, we produce it and save it in", file.name, ".\n")
            splits <- split.data(data)
            saveRDS(splits, file = file.name)
        }
    } else if (extrapol) {
      
      training.ratio <- as.numeric(extrapol.info["training.ratio"])
      variable <- extrapol.info["variable"]
      
      # Split for extrapolation
      file.name <- paste0(dir.workspace, prefix, "SplitForExtrapolation_", training.ratio, variable, ".rds")
      
      # If the file with the three different splits already exist, just read it
      if (file.exists(file.name) == T ){
        
        if(exists("splits") == F){ splits <- readRDS(file = file.name)
        cat("\nFile with data splitted for extrapolation already exists, we read it from", file.name, "and save it in object 'splits'.\n")}
        else{
          cat("\nList with data splits already exists as object 'splits' in this environment.\n")
        }
      } else {
        
        cat("\nNo data splits exist yet, we produce it and save it in", file.name, ".\n")
        splits <- split.data.ml(data, training.ratio = training.ratio, variable = variable)
        saveRDS(splits, file = file.name)
      }
    } else {
      splits <- list(data)
    }
    
    # Normalize data
    
    # Calculate mean and sd for env data for normalization with data leakage
    mean.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        mean(k, na.rm = TRUE)
    })
    sd.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        sd(k, na.rm = TRUE)
    })
    
    # Normalize the data (each folds for CV and whole data set else)
    if(CV == T | extrapol == T){
        
        # Center the splits
        centered.splits.tmp <- lapply(splits, FUN = center.data, CV = CV, extrapol = extrapol, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full = env.fact.full)
        
        # To control if the splits have the same taxa
        # # Splits don't have same taxa columns, should be harmonized
         col1 <- colnames(centered.splits.tmp$Split1$`Training data`)
        # col2 <- colnames(centered.splits.tmp$Split2$`Training data`)
        # col3 <- colnames(centered.splits.tmp$Split3$`Training data`)
        # 
        # # Make union of column names
        # col.names <- unique(c(col1, col2, col3))
        # list.taxa <- col.names[which(grepl("Occurrence.", col.names))]
         list.taxa <- col1[which(grepl("Occurrence.", col1))]
        # 
        # # Make intersection of column names (remove taxa not common to each split)
        # rem.taxa <- unique(c(col1[-which(col1 %in% col2)],
        #                     col1[-which(col1 %in% col3)],
        #                     col2[-which(col2 %in% col1)],
        #                     col2[-which(col2 %in% col3)],
        #                     col3[-which(col3 %in% col1)],
        #                     col3[-which(col3 %in% col2)]))
        # cat("\nFollowing taxa are not in all splits, consider removing them:",length(rem.taxa), rem.taxa, "\n")

        # Extract necessary information
        centered.data <- lapply(centered.splits.tmp,"[", 1:2) # only the splits without the mean, sd info
        normalization.data <- lapply(centered.splits.tmp,"[", 3:4) # the mean and sd of the splits
        
        # Normalize the folds but replace '0' and '1' by factors
        centered.data.factors <- lapply(centered.data, function(split){
            
            return(lapply(split, function(fold){
                
                cind.taxa <- which(grepl("Occurrence.",colnames(fold)))
                #Replace "0" and "1" by "absent" and "present" and convert them to factors
                for (i in cind.taxa ) {
                    fold[which(fold[,i] == 0),i] <- "absent"
                    fold[which(fold[,i] == 1),i] <- "present"
                    fold[,i] = as.factor(fold[,i])
                }
                return(fold)
            })
            )
        })

    } else {
        
        centered.data <- center.data(split = splits, CV = CV, extrapol = extrapol, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full = env.fact.full)
        normalization.data <- list("Mean" = mean.dl, "SD" = sd.dl)
        centered.data.factors <- centered.data
        
        # Replace '0' and '1' by factors
        cind.taxa <- which(grepl("Occurrence.",colnames(centered.data.factors[[1]])))
        #Replace "0" and "1" by "absent" and "present" and convert them to factors
        for (i in cind.taxa ) {
            centered.data.factors[[1]][which(centered.data.factors[[1]][,i] == 0),i] <- "absent"
            centered.data.factors[[1]][which(centered.data.factors[[1]][,i] == 1),i] <- "present"
            centered.data.factors[[1]][,i] = as.factor(centered.data.factors[[1]][,i])
        }
        
        cind.taxa <- which(grepl("Occurrence.",colnames(centered.data$`Entire dataset`)))
        list.taxa <- colnames(centered.data$`Entire dataset`)[cind.taxa]
    }
    
    # Test centering
    
    # if(CV == T){
    #     if(mean(centered.data$Split1$`Training data`$temperature) <= 0.001){
    #         cat("\nThe data is normalized.\n")
    #     }else{
    #         cat("\nThe data isn't normalized.\n")
    #         break()
    #     }
    # }else if (exists("centered.data") == T){
    #     if(mean(centered.data$`Entire dataset`$temperature) <= 0.001){
    #         cat("\nThe data is normalized.\n")
    #     }else{
    #         cat("\nThe data isn't normalized.\n")
    #         break()
    #     }
    # }else{
    #     cat("\nThe data isn't normalized.\n")
    #     
    # }
    
    # Update prevalence dataframe with new list of taxa and number of samples
    prev.inv <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa),]
    for (j in list.taxa) { 
      prev <- sum(data[,j]) / dim(data)[1]
      prev.inv[which(prev.inv$Occurrence.taxa == j), "Prevalence"] <- prev
    }
    
    prev.inv[which(grepl("Occurrence.group", prev.inv$Occurrence.taxa)), "Taxonomic.level"] <- "group"
    
    # Print information
    cat("\nSummary information of final datasets after preprocessing:\n",
        length(unique(data$SampId)), "samples,\n",
        length(unique(data$SiteId)), "sites,\n",
        length(list.taxa), "taxa (without missing values) with prevalence between 5% and 95% with taxonomic level:\n")
        
    print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
    
    preprocessed.data <- list("data" = data, "splits" = splits, "list.taxa" = list.taxa, # "rem.taxa" = if(CV){ rem.taxa },
                              "centered.data" = centered.data, "centered.data.factors" = centered.data.factors, 
                              "normalization.data" = normalization.data, "prev.inv" = prev.inv)
    return(preprocessed.data)
}

# Apply NULL model
apply.null.model <- function(data, list.taxa, prev.inv){
  
  # Make a list with the "outputs" of null model
  list.outputs <- vector(mode = 'list', length = length(list.taxa))
  names(list.outputs) <- list.taxa
  
  for (j in list.taxa){
    
    temp.output <- c("Likelihood", "Performance")
    temp.list <- vector(mode = 'list', length = length(temp.output))
    names(temp.list) <- temp.output
    
    no.pres <- sum(data[, j] == 1, na.rm = TRUE)
    no.abs  <- sum(data[, j] == 0, na.rm = TRUE)
    no.obs  <- no.pres + no.abs
    prev    <- prev.inv[prev.inv$Occurrence.taxa == j,"Prevalence"]
    likeli <- rep(c(prev, 1-prev),c(no.pres, no.abs))
    temp.list[["Likelihood"]] <- likeli
    
    st.dev <- -2 * sum(log(likeli)) / no.obs
    temp.list[["Performance"]] <- st.dev
    
    list.outputs[[j]] <- temp.list
  }
  return(list.outputs)
}

## ---- Process output from stat models to fit structure of ml models (makes plotting easier)
#JW: THE CODE IS QUITE UGLY AND DUPLICATE ATM BUT AT LEAST IT WORKS
#ECR: I won't be the one blaming you for ugly code ':D (awkward smiley)
transfrom.stat.outputs <- function(CV, stat.outputs){
    #CV = F
    #stat.outputs = stat.outputs
    
    #stat.output.list <- vector(mode = "list", length = length(stat.outputs))
    
    if ( CV == F){ # Model performance (FIT)
        stat.fit.res <- lapply(stat.outputs, function(models){
          #models <- stat.outputs[[1]]
          temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
          
          for(j in 1:length(list.taxa)){
            #j = 1
            dev.temp <- models[[2]]$deviance
            dev.temp <- dev.temp[c("Taxon", "std.deviance")]
            dev.temp <- subset(dev.temp, Taxon == list.taxa[[j]])
         
            dev.temp$Performance.train <- as.numeric(dev.temp$std.deviance)

            prop.temp <- models[[2]]$probability
            prop.temp <- subset(prop.temp, Taxon == list.taxa[[j]])
            
            prop.temp$Likelihood.train <- ifelse(prop.temp$Obs == 1, prop.temp$Pred, 1 - prop.temp$Pred)
            
            temp.list.st.dev[[j]] <- list("Prediction factors training set" = ifelse(prop.temp$Pred >= 0.5,"present","absent"),
                                          "Prediction probabilities training set" = data.frame("present" = prop.temp$Pred, "absent" = 1 - prop.temp$Pred),
                                          "Likelihood training set" = prop.temp$Likelihood.train,
                                          "Performance training set" = dev.temp$Performance.train
            )
          }
          names(temp.list.st.dev) <-  list.taxa
          #stat.output.list[[n]] <- temp.list.st.dev
          return(temp.list.st.dev)
        })
        #stat.outputs[[1]] <- stat.outputs[[1]]$deviance[c("Taxon", "std.deviance")]
        
        # for(n in 1:length(stat.outputs)){
        #     #n = 1
        #     # ECR: Maybe later fix that environment variables used here (like list.taxa) are also declared in the function
        #     # JW you are right, I'm not quite sure yet what to do about this yet
        #     temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
        #     
        #     for(j in 1:length(list.taxa)){
        #         #j = 1
        #         temp.dat.dev <- subset(stat.outputs[[n]][[2]]$deviance)
        #         
        #         temp.dat.dev <- subset(temp.dat.dev, Taxon == list.taxa[[j]])
        #         
        #         temp.dat.dev$Performance <- as.numeric(temp.dat.dev$std.deviance)
        #         
        #         temp.dat.prop <- subset(stat.outputs[[n]][[2]]$probability)
        #         temp.dat.prop <- subset(temp.dat.prop, Taxon == list.taxa[[j]])
        #         
        #         temp.dat.prop$Likelihood  <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
        #         
        #         temp.list.st.dev[[j]] <- list("Performance training set" = temp.dat.dev$Performance, "Likelihood " = temp.dat.prop$Likelihood )
        #         
        #     }
        #     names(temp.list.st.dev) <-  list.taxa
        #     stat.output.list[[n]] <- temp.list.st.dev
        #     
        # }
        # names(stat.output.list) <- names(stat.outputs)
        # return(stat.output.list)
        return(stat.fit.res) # this return statments are not really needed
    } else { # Prediction (CV)
        stat.cv.res <- lapply(stat.outputs, function(models){
            #models <- stat.outputs[[2]]
            
            stat.cv.res.splits <- lapply(models, function(split){
              #split <- models[[3]]
              temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
              
              for(j in 1:length(list.taxa)){
                #j = 1
                dev.temp <- split[[2]]$deviance
                dev.temp <- dev.temp[c("Taxon", "Type", "std.deviance")]
                dev.temp <- subset(dev.temp, Taxon == list.taxa[[j]])
                dev.temp.train <- subset(dev.temp, Type == "Training")
                dev.temp.test <- subset(dev.temp, Type == "Testing")
                
                
                dev.temp.train$Performance.train <- as.numeric(dev.temp.train$std.deviance)
                dev.temp.test$Performance.test <- as.numeric(dev.temp.test$std.deviance)

                prop.temp <- split[[2]]$probability
                prop.temp <- subset(prop.temp, Taxon == list.taxa[[j]])
                prop.temp.train <- subset(prop.temp, Type == "Training")
                prop.temp.test <- subset(prop.temp, Type == "Testing")
                
                prop.temp.train$Likelihood.train <- ifelse(prop.temp.train$Obs == 1, prop.temp.train$Pred, 1 - prop.temp.train$Pred)
                prop.temp.test$Likelihood.test <- ifelse(prop.temp.test$Obs == 1, prop.temp.test$Pred, 1 - prop.temp.test$Pred)
                
                
                temp.list.st.dev[[j]] <- list("Prediction factors training set" = ifelse(prop.temp.train$Pred >= 0.5,"present","absent"),
                                              "Prediction probabilities training set" = data.frame("present" = prop.temp.train$Pred, "absent" = 1 - prop.temp.train$Pred),
                                              "Likelihood training set" = prop.temp.train$Likelihood.train,
                                              "Performance training set" = dev.temp.train$Performance.train,
                                              
                                              "Prediction factors testing set" = ifelse(prop.temp.test$Pred >= 0.5,"present","absent"),
                                              "Prediction probabilities testing set" = data.frame("present" = prop.temp.test$Pred,"absent" = 1 - prop.temp.test$Pred),
                                              "Likelihood testing set" = prop.temp.test$Likelihood.test,
                                              "Performance testing set" = dev.temp.test$Performance.test
                                              )
                
                
                
                }
                names(temp.list.st.dev) <-  list.taxa
                return(temp.list.st.dev)})
            return(stat.cv.res.splits)
            })
            #JW: I LEAVE THIS STUFF IN RIGHT NOW FOR IDEAS BUT CAN BE DELETED LATER
            #bind rows for all three splits
            #stat.cv.res.splits[[1]]$Performance$Split <- 1
            #stat.cv.res.splits[[1]]$Prediction$Split <- 1
            
            #stat.cv.res.splits.table <- stat.cv.res.splits[[1]]
            
            # for(i in 2:length(stat.cv.res.splits)){
            #     # i = 2
            #     stat.cv.res.splits[[i]]$Performance$Split <- i
            #     stat.cv.res.splits[[i]]$Prediction$Split <- i
            #     
            #     stat.cv.res.splits.table$Performance <- rbind(stat.cv.res.splits.table$Performance, stat.cv.res.splits[[i]]$Performance)
            #     stat.cv.res.splits.table$Prediction <- rbind(stat.cv.res.splits.table$Prediction, stat.cv.res.splits[[i]]$Prediction)
            # }
            # stat.cv.res.splits.table2 <- as.data.frame(stat.cv.res.splits.table$Performance %>% group_by(Taxon) %>% summarise(Performance = mean(std.deviance, na.rm = T)))
            # stat.cv.res.splits.table3 <- as.data.frame(stat.cv.res.splits.table$Prediction %>% group_by(Taxon, SiteId) %>% summarise(Likelihood  = mean(Likelihood , na.rm = T), Pred = mean(Pred, na.rm = T)))
            #return(list("Performance" = stat.cv.res.splits.table2, "Prediction" = stat.cv.res.splits.table3))
            #return(stat.cv.res.splits)
        #})
    #     
    #     for(n in 1:length(stat.outputs)){
    #         #n = 1
    #         temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
    #         
    #         for(j in 1:length(list.taxa)){
    #             #j = 1
    #             
    #             temp.dat.dev <- subset(stat.cv.res[[n]]$Performance, Taxon == list.taxa[[j]])
    #             temp.dat.pred <- subset(stat.cv.res[[n]]$Prediction, Taxon == list.taxa[[j]])
    #             
    #             temp.list.st.dev[[j]]$Performance <- temp.dat.dev
    #             temp.list.st.dev[[j]]$Prediction <- temp.dat.pred
    #             
    #             #names(temp.list.st.dev[[2]]) <- "Performance testing set"
    #             #temp.dat.prop <- subset(stat.outputs$Split1[[2]]$probability, Type == "Testing")
    #             #temp.dat.prop$Likelihood  <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
    #             
    #             #temp.list.st.dev[[1]] <- list("Performance" = temp.dat.dev$Performance)
    #             
    #         }
    #         names(temp.list.st.dev) <-  list.taxa
    #         stat.output.list[[n]] <- temp.list.st.dev
    #         
    #     }
    #     names(stat.output.list) <- names(stat.outputs)
    # }
    
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
    # stat.cv.res <- lapply(stat_cv_nocorr, function(split){
    #   #split <- stat_cv_nocorr[[1]]
    #   dev_temp <- split[[2]]$deviance
    #   dev_temp$std.deviance <- dev_temp$std.deviance
    #   dev_temp <- dev_temp[c("Taxon", "Type", "std.deviance")]
    #   tmp_dev_test <- subset(dev_temp, Type = "Testing")
    #   return(tmp_dev_test)
    # })
    # 
    #return(stat.output.list)
    return(stat.cv.res)
    }
}

make.final.outputs.cv <- function(outputs.cv, list.models, list.taxa){
    
    outputs <- vector(mode = "list", length = length(list.models))
    names(outputs) <- list.models
    
    out <- c("Observation", #1
             "Prediction factors", #2 
             "Prediction probabilities", #3 
             "Likelihood", #4
             "Performance", #5
             "Performance splits") #6
    output.names <- paste(out, "testing set")
    
    for (l in list.models){
        
        temp.list.taxa <- vector(mode = "list", length = length(list.taxa))
        names(temp.list.taxa) <- list.taxa

        for(j in list.taxa){
            
            temp.output <- vector(mode = "list", length = length(output.names))
            names(temp.output) <- output.names
            
            for (m in output.names[3]) {
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
                perf <- outputs.cv[[n]][[l]][[j]][["Performance testing set"]]
                temp.vect[n] <- ifelse(is.numeric(perf), perf, NA)
            }
            temp.output[[5]] <- mean(temp.vect, na.rm = T)
            temp.output[[6]] <- temp.vect
            
            temp.list.taxa[[j]] <- temp.output
        }
        outputs[[l]] <- temp.list.taxa
    }
    return(outputs)
}

make.df.outputs <- function(outputs, list.models, list.taxa, 
                            list.splits = c("Split1", "Split2", "Split3"), null.model, prev.inv, CV, extrapol){
    
    if(CV | extrapol){
        outputs.cv <- outputs
        no.splits <- length(list.splits)
        no.taxa <- length(list.taxa)
        no.models <- length(list.models)
        
        # make table with perf for each split
        df.pred.perf.cv <- data.frame(matrix(nrow = no.taxa, ncol = no.splits*no.models))
        colnames(df.pred.perf.cv) <- apply(expand.grid(list.splits,list.models), 1, paste, collapse="_")
        df.pred.perf.cv$Taxa <- list.taxa
        df.fit.perf.cv <- df.pred.perf.cv

        for (s in list.splits) {
          # s <- list.splits[1]
            for (l in list.models) {
              # l <- list.models[1]
                list.taxa.temp <- names(outputs.cv[[s]][[l]])
                for (j in list.taxa.temp) {
                  # j <- list.taxa[1]
                    rind.taxa <- which(df.pred.perf.cv$Taxa == j)
                    pred.perf <- outputs.cv[[s]][[l]][[j]][["Performance testing set"]]
                    fit.perf <- outputs.cv[[s]][[l]][[j]][["Performance training set"]]
                    df.pred.perf.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(pred.perf), pred.perf, NA)
                    df.fit.perf.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(fit.perf), fit.perf, NA)
                }
            }
        }
    }
    
    # make table with mean perf across splits
    df.pred.perf <- data.frame(matrix(nrow = no.taxa, ncol = no.models))
    colnames(df.pred.perf) <- list.models
    df.pred.perf$Taxa <- list.taxa
    df.pred.perf$Prevalence <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Prevalence"]
    df.pred.perf[, "Taxonomic level"] <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"]
    df.pred.perf[,"Null_model"] <- NA
    # Add columns for explanatory power
    expl.pow <- paste0("expl.pow_", list.models)
    df.pred.perf[, expl.pow] <- NA
    df.fit.perf <- df.pred.perf
    
    for (j in list.taxa) {
        rind.taxa <- which(df.pred.perf$Taxa == j)
        df.pred.perf[rind.taxa, "Null_model"] <- null.model[[j]][["Performance"]]
        df.fit.perf[rind.taxa, "Null_model"] <- null.model[[j]][["Performance"]]
        for(l in list.models){
            if(CV | extrapol){
                splits.model <- apply(expand.grid(list.splits,l), 1, paste, collapse="_")
                # For testing/prediction
                mean.temp <- mean(as.matrix(df.pred.perf.cv[rind.taxa, splits.model]), na.rm = T)
                df.pred.perf[rind.taxa,l] <- mean.temp
                val.expl.pow <- (null.model[[j]][["Performance"]] - mean.temp) / null.model[[j]][["Performance"]]
                df.pred.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
                # df.pred.perf[rind.taxa, paste0("expl.pow_","Null_model")] <- 0
                # For training/fitting
                mean.temp <- mean(as.matrix(df.fit.perf.cv[rind.taxa, splits.model]), na.rm = T)
                df.fit.perf[rind.taxa,l] <- mean.temp
                val.expl.pow <- (null.model[[j]][["Performance"]] - mean.temp) / null.model[[j]][["Performance"]]
                df.fit.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
                # df.fit.perf[rind.taxa, paste0("expl.pow_","Null_model")] <- 0
            } else {
                perf <- outputs[[l]][[j]][["Performance training set"]]
                df.fit.perf[rind.taxa,l] <- perf
                val.expl.pow <- (null.model[[j]][["Performance"]] - perf) / null.model[[j]][["Performance"]]
                df.fit.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
                # df.fit.perf[rind.taxa, paste0("expl.pow_","Null_model")] <- 0
            }

        }
    }
    
    # Make one df with performance during training AND testing AND expl.pow
    if(CV | extrapol){
      # Merge dataframes for comparison
      common.vect <- c("Taxa", "Prevalence", "Taxonomic level", "Null_model")
      merged.df <- left_join(df.fit.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], df.pred.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], 
                             by = c(common.vect), suffix = c(".fit", ".pred"))
      
      merged.df[, paste0(list.models, ".likelihood.ratio")] <- NA
      merged.df$Big.model.diff <- NA
      merged.df$Big.pred.expl.pow.diff <- NA
      merged.df$CF0.model.diff <- NA
      merged.df$CF0.pred.expl.pow.diff <- NA
      
      
      # temp.df$Big.pred.expl.pow.diff <- apply(temp.df[,1:no.models], 1, FUN = function(x){diff(range(x))})
      # min.max <- apply(temp.df[,1:no.models], 1, FUN = function(x){range(x)})
      
      # Compute biggest difference in expl. pow. of prediction
      temp.df <- select(merged.df, "Taxa", contains("expl.pow") & contains(".pred") & !contains("diff"))
      for(j in list.taxa){
        # j <- list.taxa[1]
        # Find min and max expl. pow
        min <- min(temp.df[which(temp.df$Taxa == j), 2:(no.models+1)])
        max <- max(temp.df[which(temp.df$Taxa == j), 2:(no.models+1)])
        model.min <- sub("expl.pow_", "", sub(".pred", "", colnames(temp.df)[which(temp.df[which(temp.df$Taxa == j), ] == min)]))
        model.max <- sub("expl.pow_", "", sub(".pred", "", colnames(temp.df)[which(temp.df[which(temp.df$Taxa == j), ] == max)]))
        merged.df[which(merged.df$Taxa == j), "Big.model.diff"] <- paste(model.max, model.min, sep = "-")
        merged.df[which(merged.df$Taxa == j), "Big.pred.expl.pow.diff"] <- max - min
        
        if(grepl("CF0", colnames(temp.df)) == T){
            # Compare with CF0
            expl.pow.CF0 <- temp.df[which(temp.df$Taxa == j), which(grepl("CF0", colnames(temp.df)))]
            merged.df[which(merged.df$Taxa == j), "CF0.model.diff"] <- paste(model.max, "CF0", sep = "-")
            merged.df[which(merged.df$Taxa == j), "CF0.pred.expl.pow.diff"] <- max - expl.pow.CF0
        }
        # Compute likelihood ratio
        for (l in list.models) {
          pred <- merged.df[which(merged.df$Taxa == j), paste0(l, ".pred")]
          fit <- merged.df[which(merged.df$Taxa == j), paste0(l, ".fit")]
          merged.df[which(merged.df$Taxa == j), paste0(l, ".likelihood.ratio")] <- exp(-(pred - fit) / 2)
        }
      }
    }
    
    
    if(CV | extrapol){
        result <- list("Table predictive performance CV" = df.pred.perf.cv, "Table predictive performance" = df.pred.perf,
                       "Table fit performance CV" = df.fit.perf.cv, "Table fit performance" = df.fit.perf, "Table merged" = merged.df)
    } else {
        result <- df.fit.perf
    }
    
    return(result)
}

make.table <- function(df.pred.perf, df.fit.perf, list.models){
  list.models <- append(list.models, "Null_model")
  names(list.models) <- c()

  # calculate mean standardized deviance for the training set (i.e. quality of fit)
  table.fit.mean <- apply(df.fit.perf[list.models],2, FUN = mean)
  table.fit.mean <-t(table.fit.mean)
  table.fit.mean <- as.data.frame(table.fit.mean)
  rownames(table.fit.mean) <- "Mean std. dev. during training"
  
  # calculate corresponding standart deviation
  table.fit.sd <- apply(df.fit.perf[list.models],2, FUN = sd)
  table.fit.sd <-t(table.fit.sd)
  table.fit.sd <- as.data.frame(table.fit.sd)
  rownames(table.fit.sd) <- "SD std. dev. during training"
  
  # calculate mean standardized deviance for the testing set (i.e. quality of fit)
  table.pred.mean <- apply(df.pred.perf[list.models],2, FUN = mean)
  table.pred.mean <-t(table.pred.mean)
  table.pred.mean <- as.data.frame(table.pred.mean)
  rownames(table.pred.mean) <- "Mean std. dev. during testing"
  
  # calculate corresponding standart deviation
  table.pred.sd <- apply(df.pred.perf[list.models],2, FUN = sd)
  table.pred.sd <-t(table.pred.sd)
  table.pred.sd <- as.data.frame(table.pred.sd)
  rownames(table.pred.sd) <- "SD std. dev. during testing"
  
  # calculate mean standardized deviance for the testing set (i.e. predictive perfomance)
  names.expl.pow <- paste("expl.pow_", list.models, sep="")
  
  table.mean.exp <- apply(df.pred.perf[names.expl.pow],2, FUN = mean)
  table.mean.exp <-t(table.mean.exp)
  table.mean.exp <- as.data.frame(table.mean.exp)
  colnames(table.mean.exp) <- list.models
  rownames(table.mean.exp) <- "Mean expl. power"
  
  # calculate corresponding standart deviation
  table.sd.exp <- apply(df.pred.perf[names.expl.pow],2, FUN = sd)
  table.sd.exp <-t(table.sd.exp)
  table.sd.exp <- as.data.frame(table.sd.exp)
  colnames(table.sd.exp) <- list.models
  rownames(table.sd.exp) <- "SD expl. power"
  
  # add performance ratio
  perf.ratio <- (table.pred.mean/table.pred.sd) * table.pred.mean
  rownames(perf.ratio) <- "performance ratio"
  
  # row bind results
  table <- rbind(table.fit.mean, table.fit.sd, table.pred.mean, table.pred.sd, table.mean.exp, table.sd.exp, perf.ratio)
  
  tab1 <- table %>% gt(rownames_to_stub = T) %>% tab_header(
    title = md("**Mean predictive performance across models**") # make bold title
  ) %>%
    fmt_number(
      columns = list.models, # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  return(tab1)
}
    
make.table.species <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  
  common.vect <- c("Taxa", "Prevalence", "Null_model")
  expl.pow.vect.pred <- paste("expl.pow_", list.models, ".pred", sep="")
  expl.pow.vect.fit <- paste("expl.pow_", list.models, ".fit", sep="")
  likeli.ratio.vect <- colnames(df.merged.perf)[which(grepl("likelihood.ratio", colnames(df.merged.perf)))]
  
  # Make tables
  tmp.table <- df.merged.perf
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  # colnames(tmp.table) <- c("Taxa", "Prevalence", list.models)
  #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
  tab2 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**") # make bold title
    ) %>%
    tab_spanner(
      label = "Training",
      columns = paste0(list.models, ".fit")
    ) # %>%
    
    fit.vect <- list.models
    names(fit.vect) <- paste0(list.models, ".fit")
      
    tab2 <- tab2 %>% cols_label(
      .list = fit.vect
    ) %>%
    tab_spanner(
      label = "Testing",
      columns = paste0(list.models, ".pred")
    ) %>%
    tab_spanner(
      label = "Explanatory power for prediction",
      columns = all_of(expl.pow.vect.pred)
    ) %>%
    tab_spanner(
      label = "Explanatory power during training",
      columns = all_of(expl.pow.vect.fit)
    ) %>%
    tab_spanner(
      label = "Likelihood ratio",
      columns = all_of(likeli.ratio.vect)
    ) %>%
    tab_spanner(
      label = "Biggest expl. pow. difference",
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, CF0.model.diff, CF0.pred.expl.pow.diff)
    ) %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "CF0.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab2)
}

make.table.species.rearranged <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  no.models <- length(list.models)
  
  # Make tables
  tmp.table <- df.merged.perf[,-which(grepl("expl.pow",colnames(df.merged.perf)) & grepl(".fit",colnames(df.merged.perf)))]
  
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  # colnames(tmp.table) <- c("Taxa", "Prevalence", list.models)
  #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
  
  tab3 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**") # make bold title
    ) %>%
    cols_label(
      Null_model = "Null model",
      Big.model.diff = "Models",
      Big.pred.expl.pow.diff = "Value"
    ) %>%
    tab_spanner(
      label = "Biggest expl. pow. difference in pred. ",
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, CF0.model.diff, CF0.pred.expl.pow.diff)
    )
  for (l in 0:(no.models-1)) {
    col.group <- colnames(tmp.table)[which(grepl(list.models[no.models-l], colnames(tmp.table)) & !grepl("diff", colnames(tmp.table)) )]
    col.names <- c("Fit", "Prediction", "Expl. pow.", "Likelihood ratio")
    names(col.names) <- col.group
    
    tab3 <- tab3 %>%
      cols_move(
        columns = all_of(col.group),
        after = "Null_model"
        ) %>%
      tab_spanner(
        label = list.models[no.models-l],
        columns = all_of(col.group)
      ) %>%
      cols_label( .list = col.names
      )
  }
  
  tab3 <- tab3 %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "CF0.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab3)
}

make.table.species.rearranged.order <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  no.models <- length(list.models)
  
  # Make tables
  tmp.table <- df.merged.perf[,-which(grepl("expl.pow",colnames(df.merged.perf)) & grepl(".fit",colnames(df.merged.perf)))]
  
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  tmp.table  <- arrange(tmp.table, desc(CF0.pred.expl.pow.diff))
  
  tab3 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**") # make bold title
    ) %>%
    cols_label(
      Null_model = "Null model",
      Big.model.diff = "Models",
      Big.pred.expl.pow.diff = "Value"
    ) %>%
    tab_spanner(
      label = "Biggest expl. pow. difference in pred.",
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, CF0.model.diff, CF0.pred.expl.pow.diff)
    )
  for (l in 0:(no.models-1)) {
    col.group <- colnames(tmp.table)[which(grepl(list.models[no.models-l], colnames(tmp.table)) & !grepl("diff", colnames(tmp.table)))]
    col.names <- c("Fit", "Prediction", "Expl. pow.", "Likelihood ratio")
    names(col.names) <- col.group
    
    tab3 <- tab3 %>%
      cols_move(
        columns = all_of(col.group),
        after = "Null_model"
      ) %>%
      tab_spanner(
        label = list.models[no.models-l],
        columns = all_of(col.group)
      ) %>%
      cols_label( .list = col.names
      )
  }
  
  tab3 <- tab3 %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "CF0.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab3)
}










