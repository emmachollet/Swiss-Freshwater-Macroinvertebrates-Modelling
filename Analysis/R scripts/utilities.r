## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##               Various functions used in all the scripts
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Data wrangling ----

# Split the data in 3 three-folds for cross-validation
split.data.cv <- function(data){
  
  set.seed(2021)  
  
  folds <- groupKFold(data$SiteId, 3) # Keep same sites in same split to avoid data leakage
  
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


# Split data in training and testing datasets for ODGation
split.data.odg <- function(data, training.ratio, variable = "random", bottom = T){
  
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

# 
standardize.data <- function(data, split, CV, ODG, dl, mean.dl, sd.dl, env.fact.full){
  
  training.data <- split[[1]]
  
  # extract column names to access predictors and invertebrate data separately
  inv.names <- colnames(select(training.data, contains("Occurrence.")))
  env.names <- env.fact.full
  info.names <- colnames(select(training.data, - all_of(inv.names), - all_of(env.names)))
  
  # convert to numeric to count observations across sites
  inv.data <- as.data.frame(apply(training.data[, inv.names],2,as.numeric))
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
    
  } else {
    
    for(env in env.names){
      training.data[env] <- training.data[env] -  mean.env.cond[env]
    }
    
    for(env in env.names){
      training.data[env] <- training.data[env] / sd.env.cond[env]
    }
  }

  # Re-calculate temp2 and velocity2 with scaled variables
  training.data$temperature2 <- training.data$temperature^2
  training.data$velocity2 <- training.data$velocity^2
  
  if(CV == F & ODG == F){
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
      
    } else {
      
      for(env in env.names){
        testing.data[env] <- testing.data[env] -  mean.env.cond[env]
      }
      
      for(env in env.names){
        testing.data[env] <- testing.data[env] / sd.env.cond[env]
      }
    }
    
    # Re-calculate temp2 and velocity2 with scaled variables
    testing.data$temperature2 <- testing.data$temperature^2
    testing.data$velocity2 <- testing.data$velocity^2
    
    return(list("Training data" = training.data, "Testing data" = testing.data, "Mean" = mean.env.cond, "SD" = sd.env.cond))
  }
}


preprocess.data <- function(data.env, data.inv, prev.inv, remove.na.env.fact = T, remove.na.taxa = T, prev.restrict, env.fact.full, dir.workspace, BDM, dl, CV, ODG, ODG.info = c(), lme.temp = lme.temp){
    
  # Print information
  cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))
  list.taxa <- colnames(data.inv)[cind.taxa]
  cat("\nSummary information of original datasets before preprocessing:\n",
      length(unique(data.inv$SampId)), "samples,\n",
      length(unique(data.inv$SiteId)), "sites,\n",
      length(list.taxa), "taxa (including missing values) with taxonomic level:\n")
  
  print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
  
  if(remove.na.env.fact){
    
    # Drop rows with incomplete influence factors
    rind <- !apply(is.na(data.env[,env.fact.full]), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\n", sum(!rind),"sites/samples excluded because of incomplete environmental factor.\n"))
    data.env <- data.env[rind,]
    dim(data.env)
    
  }
  
  data.inv <- data.inv[data.inv$SampId %in% data.env$SampId, ]
  
  
  if(remove.na.taxa){
    
    # Drop taxa with too many NAs
    too.many.na <- c()
    threshold <- ifelse(BDM, 100, 200)
    for(i in cind.taxa){
      if(sum(is.na(data.inv[,i])) > threshold){ too.many.na <- c(too.many.na, i)}
    }
    cat("\nThe following", length(too.many.na), "taxa are excluded because they have more than ", threshold, "NAs :\n", gsub("Occurrence.", "", colnames(data.inv)[too.many.na]), "\n")
    data.inv <- data.inv[, -too.many.na]
    
    # Drop rows with incomplete taxa or influence factors
    rind <- !apply(is.na(data.inv), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\n", sum(!rind),"sites/samples excluded because of incomplete taxa.\n"))
    data.inv <- data.inv[rind,]
    dim(data.inv)
    
    # Delete taxa with only presence or absence (happens because we deleted rows)
    # Delete taxa with less than 5% of only present or absent
    cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
    no.samples <- nrow(data.inv)
    prev.perc <- prev.restrict*no.samples
    cind.rem <- c()
    for (j in cind.taxa) {
        if(sum(data.inv[,j]) < prev.perc | sum(data.inv[,j]) > no.samples - prev.perc){ cind.rem <- c(cind.rem, j) }
    }
    cat("\nThe following", length(cind.rem), "taxa are excluded because they have less than", 100*prev.restrict, "% or more than", 100 - 100*prev.restrict, "% of prevalence:\n", gsub("Occurrence.", "", colnames(data.inv)[cind.rem]), "\n")
    if(length(cind.rem) != 0){ data.inv <- data.inv[,-cind.rem] }
    dim(data.inv)
    
  }
  
  # Update prevalence dataframe with new list of taxa and number of samples
  cind.taxa <- which(grepl("Occurrence.", colnames(data.inv))) # update list.taxa
  list.taxa <- colnames(data.inv)[cind.taxa]
  prev.inv <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa),]
  for (j in list.taxa) { 
    prev <- sum(data.inv[,j]) / dim(data.inv)[1]
    prev.inv[which(prev.inv$Occurrence.taxa == j), "Prevalence"] <- prev
  }
  prev.inv[which(grepl("Occurrence.group", prev.inv$Occurrence.taxa)), "Taxonomic.level"] <- "group"
  
  # Reorder taxa by prevalence
  list.taxa <- list.taxa[order(match(list.taxa, prev.inv$Occurrence.taxa))]
  
  # Merge data sets
  data <- data.env[, c("SiteId", "SampId", "X", "Y", "RiverBasin", "Region", env.fact.full)] %>%
    right_join(data.inv[, c(which(colnames(data.inv) %in% c("SiteId", "SampId", list.taxa)))], by = c("SiteId", "SampId"))
  # Reorder taxa columns by prevalence
  data <- data[,c(which(!(colnames(data) %in% list.taxa)), match(list.taxa, colnames(data)))]
  dim(data)

    # Split data
    prefix <- ifelse(BDM, "BDM_", "All_")
    info.file.name <- paste0(prefix,
                             dim(data)[1], "samples_",
                             ifelse(CV, "CV_",
                                    ifelse(ODG, paste(c("ODG", paste(c(ODG.info,"_"), collapse = "")), collapse = "_"),
                                           "FIT_")),
                             ifelse(lme.temp, "lme_temp_", ""),
                             "")
    if(CV){
        
        # Split for CV
        file.name <- paste0(dir.workspace, info.file.name, "Splits.rds")
        
        # If the file with the splits exist, we read it, otherwise we split the data
        if (file.exists(file.name)){
            
            if(exists("splits") == F){ splits <- readRDS(file = file.name)
            cat("\nFile with data splits already exists, we read it from", file.name, "and save it in object 'splits'.\n")}
            else{
                cat("\nList with data splits already exists as object 'splits' in this environment.\n")
            }
        } else {
            
            cat("\nNo data splits exist yet, we produce it and save it in", file.name, ".\n")
            splits <- split.data.cv(data)
            saveRDS(splits, file = file.name)
        }
    } else if (ODG) {
      
      training.ratio <- as.numeric(ODG.info["training.ratio"])
      variable <- ODG.info["variable"]
      
      # Split for ODG
      file.name <- paste0(dir.workspace, info.file.name, "Split_", training.ratio, variable, ".rds")
      
      # If the file with the three different splits already exist, just read it
      if (file.exists(file.name) == T ){
        
        if(exists("splits") == F){ splits <- readRDS(file = file.name)
        cat("\nFile with data splitted for ODG already exists, we read it from", file.name, "and save it in object 'splits'.\n")}
        else{
          cat("\nList with data splits already exists as object 'splits' in this environment.\n")
        }
      } else {
        
        cat("\nNo data splits exist yet, we produce it and save it in", file.name, ".\n")
        splits <- split.data.odg(data, training.ratio = training.ratio, variable = variable)
        saveRDS(splits, file = file.name)
      }
    } else {
      splits <- list(data)
    }
    
    # Standardize data
    
    # Calculate mean and sd for env data for normalization with data leakage
    mean.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        mean(k, na.rm = TRUE)
    })
    sd.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        sd(k, na.rm = TRUE)
    })
    
    # Normalize the data (each folds for CV and whole data set else)
    if(CV | ODG){
        
        # Center the splits
        centered.splits.tmp <- lapply(splits, FUN = standardize.data, CV = CV, ODG = ODG, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full = env.fact.full)
        
        # Extract necessary information
        standardized.data <- lapply(centered.splits.tmp,"[", 1:2) # only the splits without the mean, sd info
        normalization.data <- lapply(centered.splits.tmp,"[", 3:4) # the mean and sd of the splits
        
        # Normalize the folds but replace '0' and '1' by factors
        standardized.data.factors <- lapply(standardized.data, function(split){
            
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
        
        standardized.data <- standardize.data(split = splits, CV = CV, ODG = ODG, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full = env.fact.full)
        normalization.data <- list("Mean" = mean.dl, "SD" = sd.dl)
        standardized.data.factors <- standardized.data
        
        # Replace '0' and '1' by factors
        cind.taxa <- which(grepl("Occurrence.",colnames(standardized.data.factors[[1]])))
        #Replace "0" and "1" by "absent" and "present" and convert them to factors
        for (i in cind.taxa ) {
            standardized.data.factors[[1]][which(standardized.data.factors[[1]][,i] == 0),i] <- "absent"
            standardized.data.factors[[1]][which(standardized.data.factors[[1]][,i] == 1),i] <- "present"
            standardized.data.factors[[1]][,i] = as.factor(standardized.data.factors[[1]][,i])
        }
        
        cind.taxa <- which(grepl("Occurrence.",colnames(standardized.data$`Entire dataset`)))
        list.taxa <- colnames(standardized.data$`Entire dataset`)[cind.taxa]
    }
    
    # Test centering
    
    # if(CV == T){
    #     if(mean(standardized.data$Split1$`Training data`$temperature) <= 0.001){
    #         cat("\nThe data is normalized.\n")
    #     }else{
    #         cat("\nThe data isn't normalized.\n")
    #         break()
    #     }
    # }else if (exists("standardized.data") == T){
    #     if(mean(standardized.data$`Entire dataset`$temperature) <= 0.001){
    #         cat("\nThe data is normalized.\n")
    #     }else{
    #         cat("\nThe data isn't normalized.\n")
    #         break()
    #     }
    # }else{
    #     cat("\nThe data isn't normalized.\n")
    #     
    # }
    
    
    # Print information
    cat("\nSummary information of final datasets after preprocessing:\n",
        length(unique(data$SampId)), "samples,\n",
        length(unique(data$SiteId)), "sites,\n",
        length(list.taxa), "taxa (without missing values) with prevalence between", 
        100*prev.restrict, "% or more than", 
        100 - 100*prev.restrict, "% of prevalence with taxonomic level:\n")
    
    print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
    
    preprocessed.data <- list("data" = data, "splits" = splits, "list.taxa" = list.taxa, # "rem.taxa" = if(CV){ rem.taxa },
                              "standardized.data" = standardized.data, "standardized.data.factors" = standardized.data.factors, 
                              "normalization.data" = normalization.data, "prev.inv" = prev.inv)
    return(preprocessed.data)
}

## ---- Models utilities ----

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
transfrom.stat.outputs <- function(stat.outputs, list.taxa, CV, ODG){
    # CV = F
    # stat.outputs = stat.outputs
    
    # stat.output.list <- vector(mode = "list", length = length(stat.outputs))
    
    if ( CV == F & ODG == F ){ # Model performance (FIT)
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
            
            temp.list.st.dev[[j]] <- list("Trained model" = models[[1]],
              
                                          "Prediction factors training set" = ifelse(prop.temp$Pred >= 0.5,"present","absent"),
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
    } else { # Prediction (CV or ODG)
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
                
                
                temp.list.st.dev[[j]] <- list("Trained model" = split[[1]],

                                              "Prediction factors training set" = ifelse(prop.temp.train$Pred >= 0.5,"present","absent"),
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
    # # parm   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
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
                            list.splits = c("Split1", "Split2", "Split3"), null.model, prev.inv, CV, ODG){
    
    if(CV | ODG){
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
            if(CV | ODG){
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
    if(CV | ODG){
      # Merge dataframes for comparison
      common.vect <- c("Taxa", "Prevalence", "Taxonomic level", "Null_model")
      merged.df <- left_join(df.fit.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], df.pred.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], 
                             by = c(common.vect), suffix = c(".fit", ".pred"))
      
      merged.df[, paste0(list.models, ".likelihood.ratio")] <- NA
      merged.df$Big.model.diff <- NA
      merged.df$Big.pred.expl.pow.diff <- NA
      merged.df$chGLM.model.diff <- NA
      merged.df$chGLM.pred.expl.pow.diff <- NA
      
      
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
        
        # Compare with chGLM
        if("chGLM" %in% colnames(temp.df)){
          expl.pow.chGLM <- temp.df[which(temp.df$Taxa == j), which(grepl("chGLM", colnames(temp.df)))]
          merged.df[which(merged.df$Taxa == j), "chGLM.model.diff"] <- paste(model.max, "chGLM", sep = "-")
          merged.df[which(merged.df$Taxa == j), "chGLM.pred.expl.pow.diff"] <- max - expl.pow.chGLM
        }
        
        # Compute likelihood ratio
        for (l in list.models) {
          pred <- merged.df[which(merged.df$Taxa == j), paste0(l, ".pred")]
          fit <- merged.df[which(merged.df$Taxa == j), paste0(l, ".fit")]
          merged.df[which(merged.df$Taxa == j), paste0(l, ".likelihood.ratio")] <- exp(-(pred - fit) / 2)
        }
      }
    }
    
    
    if(CV | ODG){
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
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
    ) %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
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
  # tmp.table <- df.merged.perf[,-which(grepl("expl.pow",colnames(df.merged.perf)) & grepl(".fit",colnames(df.merged.perf)))]
  tmp.table <- df.merged.perf[,c(1, 2, which((grepl("expl.pow_",colnames(df.merged.perf)) & grepl(".pred",colnames(df.merged.perf))) | grepl("likelihood",colnames(df.merged.perf)))) ]
  
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  # colnames(tmp.table) <- c("Taxa", "Prevalence", list.models)
  #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
  
  tab3 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**")) # make bold title
    # ) %>%
    # cols_label(
    #   Null_model = "Null model",
    #   Big.model.diff = "Models",
    #   Big.pred.expl.pow.diff = "Value"
    # ) %>%
    # tab_spanner(
    #   label = "Biggest expl. pow. difference in pred. ",
    #   columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
    # )
  for (l in 0:(no.models-1)) {
    col.group <- colnames(tmp.table)[which(grepl(list.models[no.models-l], colnames(tmp.table)) & !grepl("diff", colnames(tmp.table)) )]
    col.names <- c("Expl. pow.", "Likelihood ratio")
    names(col.names) <- col.group
    
    tab3 <- tab3 %>%
      cols_move(
        columns = all_of(col.group),
        after = "Prevalence"
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
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
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
  tmp.table  <- arrange(tmp.table, desc(chGLM.pred.expl.pow.diff))
  
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
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
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
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
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

pred.stat.models <- function(res.extracted, matrix.predictors){
  # res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
  
  # Name dimensions of parameters within stanfit object
  # colnames(res.extracted[["alpha_taxa"]]) <- list.taxa
  # colnames(res.extracted[["mu_beta_comm"]]) <- env.fact.full
  # colnames(res.extracted[["sigma_beta_comm"]]) <- env.fact.full
  
  # Extract inputs (x), observations (y), and parameters at maximum posterior
  # x <- as.matrix(env.cond[,env.fact.full])
  # colnames(x) <- env.fact.full
  # y <- as.matrix(occur.taxa)
  ind.maxpost <- which.max(res.extracted[["lp__"]])
  mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
  sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
  mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
  sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
  alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
  beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
  
  # Name dimensions of maximum posterior community parameters
  # names(mu.beta.comm.maxpost) <- env.fact.full
  # names(sigma.beta.comm.maxpost) <- env.fact.full
  # rownames(beta.taxa.maxpost) <- env.fact.full
  # 
  ### Calibration results
  # Check if site effects AND latent variables are disabled
  z <- matrix(rep(alpha.taxa.maxpost,nrow(matrix.predictors)),nrow=nrow(matrix.predictors),byrow=TRUE) + 
    matrix.predictors%*%beta.taxa.maxpost
  
  p.maxpost <- 1/(1+exp(-z))
  
  return(p.maxpost)
  # res.extracted   <- rstan::extract(model,permuted=TRUE,inc_warmup=FALSE)
  # 
  # # extract taxon and env variable names
  # #output <- list("deviance" = tibble(), "probability" = tibble(), "parameters" = tibble()) #this is where the data is gathered in the end for the return
  # #training.data <- data.splits[[1]]
  # 
  # inv.names <- list.taxa
  # env.names <- colnames(env.fact.test)
  # #env.names <- colnames(select(training.data, colnames(training.data), -contains("Occurrence.")))
  # 
  # env.cond <- env.fact.test # useless
  # # occur.taxa <- training.data[,inv.names]
  # # comm.corr = comm.corr
  # # inf.fact <- colnames(select(env.cond, -SiteId, -SampId, -X, -Y))
  # # inf.fact <- colnames(select(env.cond, -SiteId, -SampId, -X, -Y))
  # 
  # 
  # # sites <- occur.taxa$SiteId
  # # samples <- occur.taxa$SampId
  # 
  # # n.sites <- length(unique(sites))
  # # n.samples <- length(samples)
  # 
  # # occur.taxa$SiteId <- NULL
  # # occur.taxa$SampId <- NULL
  # 
  # # occur.taxa <- occur.taxa[, ind]
  # # n.taxa <- ncol(occur.taxa)
  # # 
  # # replace NA with -1
  # # occur.taxa.na.encod <- occur.taxa
  # # for(j in 1:ncol(occur.taxa.na.encod)) occur.taxa.na.encod[,j] <- ifelse(is.na(occur.taxa.na.encod[,j]),
  # #                                                                         -1, occur.taxa.na.encod[,j])
  # # extract model parameters
  # colnames(res.extracted[["alpha_taxa"]]) <- list.taxa
  # 
  # colnames(res.extracted[["mu_beta_comm"]]) <- env.names
  # colnames(res.extracted[["sigma_beta_comm"]]) <- env.names
  # 
  # # Extract inputs (x), observations (y), and parameters at maximum posterior
  # x <- as.matrix(env.fact.test)
  # colnames(x) <- env.names
  # #y <- as.matrix(occur.taxa)
  # ind.maxpost <- which.max(res.extracted[["lp__"]])
  # mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
  # sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
  # mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
  # sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
  # alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
  # beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
  # 
  # # Name dimensions of maximum posterior community parameters
  # names(mu.beta.comm.maxpost) <- env.names
  # names(sigma.beta.comm.maxpost) <- env.names
  # rownames(beta.taxa.maxpost) <- env.names
  # 
  # ### Calibration results
  # alpha.mat <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE)
  # z <- alpha.mat + x%*%beta.taxa.maxpost
  # 
  # p.maxpost <- 1/(1+exp(-z))
  # 
  # train.p <- as_tibble(p.maxpost, stringsAsFactors = FALSE)
  # colnames(train.p) <- list.taxa
  # #train.p$SiteId <- sites
  # #train.p$SampId <- samples
  # #train.p <- gather(train.p, Taxon, Pred, -SiteId, -SampId)
  # #train.p <- left_join(train.p, train.obs, by = c("SiteId", "SampId", "Taxon"))
  # 
  # return(train.p[, taxon])
}


