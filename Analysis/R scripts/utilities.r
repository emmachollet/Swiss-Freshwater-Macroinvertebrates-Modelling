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
  
  #split <- splits[[1]]
  # split <- list(data)
  training.data <- split[[1]]
  
  #extract column names to access predictors and invertebrate data seperately
  
  
  #inv.names <- colnames(select(training.data, contains("Occurrence.")))
  #env.names <- colnames(select(training.data, - all_of(inv.names)))
  
  inv.names <- colnames(select(training.data, contains("Occurrence.")))
  env.names <- env.fact.full
  info.names <- colnames(select(training.data, - all_of(inv.names), - all_of(env.names)))
  
  # convert to numeric to count observations across sites
  inv.data <- as.data.frame(apply(training.data[, inv.names],2,as.numeric))
  
  # drop taxa without observations or only presence at the selected sites
  ind <- apply(inv.data,2,sum, na.rm = T) <= 0
  inv.data <- inv.data[, !ind]
  n.taxa <- ncol(inv.data)
  
  ind <-  apply(inv.data,2,sum, na.rm = T) == ncol(inv.data)
  inv.data <- inv.data[, !ind]
  
  # spot taxa columns with too many NA
  too.many.na <- c()
  for(i in 1:dim(inv.data)[2]){
      if(sum(is.na(inv.data[,i])) > 200){ too.many.na <- c(too.many.na, i)}
  }
  
  # remove col with too many NA
  inv.data <- inv.data[, -too.many.na]
  
  # remove remaining NA
  inv.data <- na.omit(inv.data)
  
  
  n.taxa <- ncol(inv.data)
  
  names.selected <- colnames(inv.data)
  training.data <- cbind(training.data[, info.names],training.data[,env.names],training.data[,names.selected])
  
  #if center is true substract the mean of each predictor, check if its divded by sd, I added the division by sd
  mean.env.cond <- apply(select(training.data, all_of(env.names)), 2, function(k){
    #mean.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
    mean(k, na.rm = TRUE)
  })
  sd.env.cond <- apply(select(training.data, all_of(env.names)), 2, function(k){
    sd(k, na.rm = TRUE)
  })
    
      
  
    for(env in env.names){
    #cat(env)
    #env <- env.names[1]
    training.data[env] <- training.data[env] -  mean.env.cond[env]
    }
  
    for(env in env.names){
      #env <- env.names[1]
      training.data[env] <- training.data[env] / sd.env.cond[env]
    }

  if(CV == F){
    return(list( "Entire dataset" = training.data))
  }
  else{
    testing.data <- split[[2]]# load testing data 
    
    #Here I make sure that the same species are dropped from the training and testing set.
    testing.data <- cbind(testing.data[,info.names],testing.data[,env.names],testing.data[,names.selected])
    
    if(dl == T){
      
      for(env in env.names){
        # env <- env.names[1]

        testing.data[env] <- testing.data[env] -  mean.dl[env]
        
      }
      for(env in env.names){
        # env <- env.names[1]
        testing.data[env] <- testing.data[env] / sd.dl[env]
        
      }
      
    }else{
      
      for(env in env.names){
        # env <- env.names[1]
        testing.data[env] <- testing.data[env] -  mean.env.cond[env]
        
      }
      
      for(env in env.names){
        # env <- env.names[1]
        testing.data[env] <- testing.data[env] / sd.env.cond[env]
      
      }
      
    }
    return(list("Training data" = training.data, "Testing data" = testing.data, "Mean" = mean.env.cond, "SD" = sd.env.cond))
  }
}
