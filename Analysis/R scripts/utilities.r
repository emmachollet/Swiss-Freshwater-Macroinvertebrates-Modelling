## ---- Split the data in 3 training and testing data sets for CV ----
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

center.data <- function(data, split, CV, dl, mean.dl, sd.dl, env.fact.full){
  
  training.data <- split[[1]]
  
  # extract column names to access predictors and invertebrate data separately
  inv.names <- colnames(select(training.data, contains("Occurrence.")))
  env.names <- env.fact.full
  info.names <- colnames(select(training.data, - all_of(inv.names), - all_of(env.names)))
  
  # convert to numeric to count observations across sites
  inv.data <- as.data.frame(apply(training.data[, inv.names],2,as.numeric))
  
  # WAS ADDED TO PRE-PROCESSING COULD BE REMOVED ####
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
  
    for(env in env.names){
    training.data[env] <- training.data[env] -  mean.env.cond[env]
    }
  
    for(env in env.names){
      training.data[env] <- training.data[env] / sd.env.cond[env]
    }

  if(CV == F){
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

## ---- Process output from stat models to fit structure of ml models (makes plotting easier)
transfrom.stat.outputs <- function(CV, stat.outputs){
  #CV = T
  #stat.outputs = stat.outputs
  
  stat.output.list <- vector(mode = "list", length = length(stat.outputs))
  
  if ( CV == F){# Model performance (FIT)
    
    #stat.outputs[[1]] <- stat.outputs[[1]]$deviance[c("Taxon", "std.deviance")]
    
    for(n in 1:length(stat.outputs)){
      #n = 1
      temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
      
      for(j in 1:length(list.taxa)){
        #j = 1
        temp.dat.dev <- subset(stat.outputs[[n]][[2]]$deviance)
        
        temp.dat.dev <- subset(temp.dat.dev, Taxon == list.taxa[[j]])
        
        temp.dat.dev$Performance <- as.numeric(temp.dat.dev$std.deviance)
        
        temp.dat.prop <- subset(stat.outputs[[n]][[2]]$probability)
        temp.dat.prop <- subset(temp.dat.prop, Taxon == list.taxa[[j]])
        temp.dat.prop$Likelihood  <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
        
        temp.list.st.dev[[j]] <- list("Performance training set" = temp.dat.dev$Performance, "Likelihood " = temp.dat.prop$Likelihood )
        
      }
      names(temp.list.st.dev) <-  list.taxa
      stat.output.list[[n]] <- temp.list.st.dev
      
    }
    names(stat.output.list) <- names(stat.outputs)
    return(stat.output.list)
  }else{ # Prediction (CV)
    stat.cv.res <- lapply(stat.outputs, function(models){
      #models <- stat.outputs[[1]]
      
      stat.cv.res.splits <- lapply(models, function(split){
        #split <- models[[1]]
        dev_temp <- split[[2]]$deviance
        dev_temp$std.deviance <- dev_temp$std.deviance
        dev_temp <- dev_temp[c("Taxon", "Type", "std.deviance")]
        tmp_dev_test <- subset(dev_temp, Type == "Testing")
        
        temp.dat.prop <- split[[2]]$probability
        temp.dat.prop$Likelihood  <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
        tmp_prop_test <- subset(temp.dat.prop, Type == "Testing")
        
        return(list("Performance" = tmp_dev_test, "Prediction" = tmp_prop_test))
      })
      #bind rows for all three splits
      stat.cv.res.splits[[1]]$Performance$Split <- 1
      stat.cv.res.splits[[1]]$Prediction$Split <- 1
      
      stat.cv.res.splits.table <- stat.cv.res.splits[[1]]
      
      for(i in 2:length(stat.cv.res.splits)){
        # i = 2
        stat.cv.res.splits[[i]]$Performance$Split <- i
        stat.cv.res.splits[[i]]$Prediction$Split <- i
        
        stat.cv.res.splits.table$Performance <- rbind(stat.cv.res.splits.table$Performance, stat.cv.res.splits[[i]]$Performance)
        stat.cv.res.splits.table$Prediction <- rbind(stat.cv.res.splits.table$Prediction, stat.cv.res.splits[[i]]$Prediction)
      }
      stat.cv.res.splits.table2 <- as.data.frame(stat.cv.res.splits.table$Performance %>% group_by(Taxon) %>% summarise(Performance = mean(std.deviance, na.rm = T)))
      stat.cv.res.splits.table3 <- as.data.frame(stat.cv.res.splits.table$Prediction %>% group_by(Taxon, SiteId) %>% summarise(Likelihood  = mean(Likelihood , na.rm = T), Pred = mean(Pred, na.rm = T)))
      return(list("Performance" = stat.cv.res.splits.table2, "Prediction" = stat.cv.res.splits.table3))  
    })
    
    for(n in 1:length(stat.outputs)){
      #n = 1
      temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
      
      for(j in 1:length(list.taxa)){
        #j = 1
        
        temp.dat.dev <- subset(stat.cv.res[[n]]$Performance, Taxon == list.taxa[[j]])
        temp.dat.pred <- subset(stat.cv.res[[n]]$Prediction, Taxon == list.taxa[[j]])
        
        temp.list.st.dev[[j]]$Performance <- temp.dat.dev
        temp.list.st.dev[[j]]$Prediction <- temp.dat.pred
        
        #names(temp.list.st.dev[[2]]) <- "Performance testing set"
        #temp.dat.prop <- subset(stat.outputs$Split1[[2]]$probability, Type == "Testing")
        #temp.dat.prop$Likelihood  <- ifelse(temp.dat.prop$Obs == 1, temp.dat.prop$Pred, 1 - temp.dat.prop$Pred)
        
        #temp.list.st.dev[[1]] <- list("Performance" = temp.dat.dev$Performance)
        
      }
      names(temp.list.st.dev) <-  list.taxa
      stat.output.list[[n]] <- temp.list.st.dev
      
    }
    names(stat.output.list) <- names(stat.outputs)
  }
  
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
  return(stat.output.list)

preprocess.data <- function(data.env, data.inv, env.fact.full, dir.workspace, BDM, dl, CV){
    
    # Drop columns with (taxa with) too many NAs
    too.many.na <- c()
    for(i in 1:dim(data.inv)[2]){
        if(sum(is.na(data.inv[,i])) > 200){ too.many.na <- c(too.many.na, i)}
    }
    cat("\nThe following", length(too.many.na), "taxa are excluded because too many missing information:", colnames(data.inv)[too.many.na], "\n")
    data.inv <- data.inv[, -too.many.na]
    
    # Merge data sets
    cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
    data <- data.env[, c("SiteId", "SampId", "X", "Y", env.fact.full)] %>%
        left_join(data.inv[, c(1, 2, cind.taxa)], by = c("SiteId", "SampId"))
    dim(data)
    
    # Drop rows with incomplete taxa or influence factors
    rind <- !apply(is.na(data), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\n", sum(!rind),"sites/samples excluded because of incomplete taxa or influence factors.\n"))
    data <- data[rind,]
    
    # Delete taxa with only presence or absence (happens because we deleted rows)
    cind.taxa <- which(grepl("Occurrence.",colnames(data)))
    no.samples <- nrow(data)
    cind.rem <- c()
    for (j in cind.taxa) {
        if(sum(data[,j]) == 0 | sum(data[,j]) == no.samples){ cind.rem <- c(cind.rem, j) }
    }
    cat("\nThe following", length(cind.rem), "taxa are excluded because only present or absent after preprocessing:", colnames(data)[cind.rem], "\n")
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
    } else {
        splits <- list(data)
    }
    
    # Normalize data
    
    # Calculate mean and sd for env data for normalisation with data leakage
    mean.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        mean(k, na.rm = TRUE)
    })
    sd.dl <- apply(select(data, all_of(env.fact.full)), 2, function(k){
        sd(k, na.rm = TRUE)
    })
    
    # Normalize the data (each folds for CV and whole data set else)
    if(CV == T){
        
        # Center the splits
        centered.splits.tmp <- lapply(splits, FUN = center.data, CV = CV, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full = env.fact.full)
        # centered.splits.tmp <- lapply(splits, FUN = center.data, CV = CV)
        
        # Splits don't have same taxa columns, should be harmonized
        col1 <- colnames(centered.splits.tmp$Split1$`Training data`)
        col2 <- colnames(centered.splits.tmp$Split2$`Training data`)
        col3 <- colnames(centered.splits.tmp$Split3$`Training data`)
        
        # Make union of column names
        col.names <- unique(c(col1, col2, col3))
        list.taxa <- col.names[which(grepl("Occurrence.",col.names))]
        
        
        # Make intersection of column names (remove taxa not common to each split)
        rem.taxa <- unique(c(col1[-which(col1 %in% col2)],
                            col1[-which(col1 %in% col3)],
                            col2[-which(col2 %in% col1)],
                            col2[-which(col2 %in% col3)],
                            col3[-which(col3 %in% col1)],
                            col3[-which(col3 %in% col2)]))

        cat("\nFollowing taxa are not in all splits, consider removing them:",length(rem.taxa), rem.taxa, "\n")
        # for (n in 1:length(centered.splits.tmp)) {
        #     for (m in 1:2) {
        #         centered.splits.tmp[[n]][[m]] <- centered.splits.tmp[[n]][[m]][, -which(colnames(
        #             centered.splits.tmp[[n]][[m]]) %in% rem.taxa)]
        #     }
        # }
        # cind.taxa <- which(grepl("Occurrence.",colnames(centered.data$Split1$`Training data`)))
        # list.taxa <- colnames(centered.data$Split1$`Training data`)[cind.taxa]

        
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
        
        centered.data <- center.data(split = splits, CV = CV, data = data, dl = dl, mean.dl = mean.dl, sd.dl = sd.dl, env.fact.full)
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
    if(CV == T){
        if(mean(centered.data$Split1$`Training data`$temperature) <= 0.001){
            cat("\nThe data is normalized.\n")
        }else{
            cat("\nThe data isn't normalized.\n")
            break()
        }
    }else if (exists("centered.data") == T){
        if(mean(centered.data$`Entire dataset`$temperature) <= 0.001){
            cat("\nThe data is normalized.\n")
        }else{
            cat("\nThe data isn't normalized.\n")
            break()
        }
    }else{
        cat("\nThe data isn't normalized.\n")
        
    }
    
    preprocessed.data <- list("data" = data, "splits" = splits, "list.taxa" = list.taxa, "rem.taxa" = if(CV){ rem.taxa },
                              "centered.data" = centered.data, "centered.data.factors" = centered.data.factors, 
                              "normalization.data" = normalization.data)
    return(preprocessed.data)
}
