

## ---- Split the data in 3 training and testing data sets for CV ----

# data <-- data.inv
# sd <-- standard deviation of the sizes of the three folds, put 1
split.data <- function(data,sd){
    #data <- inv.full
    #sd <- 1
    #set.seed(2021)
    repeat{
        inv.data <- data[sample(nrow(data)),]
        
        #
        #Order sites by the number of samples taken there
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
center.splits <- function(split,cv){
    #split <- list(inv.occ)
    occur.taxa <- split[[1]]
    inf.fact <- colnames(select(env.cond.orig, -SiteId, -SampId))
    
    # join the environmental conditions to the occurrence data
    env.cond <- left_join(occur.taxa[, c("SiteId", "SampId")], env.cond.orig, by = c("SiteId", "SampId"))
    
    # Pre-process data ####
    # drop rows with incomplete influence factors:
    ind <- !apply(is.na(env.cond[,inf.fact]),1,FUN=any)
    ind <- ifelse(is.na(ind),FALSE,ind)
    occur.taxa <- occur.taxa[ind, ]
    env.cond   <- env.cond[ind, c("SiteId", "SampId", inf.fact)]
    print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))
    
    sites <- occur.taxa$SiteId
    samples <- occur.taxa$SampId
    
    n.sites <- length(unique(sites))
    n.samples <- length(samples)
    
    #occur.taxa$SiteId <- NULL
    #occur.taxa$SampId <- NULL
    occur.taxa <- as.data.frame(occur.taxa)
    # drop TAXA without observations at the selected sites:
    ind <- apply(occur.taxa[, !(colnames(occur.taxa) %in% c("SiteId", "SampId"))],2,sum, na.rm = TRUE) <= 0
    occur.taxa <- occur.taxa[, !ind]
    n.taxa <- ncol(occur.taxa[, !(colnames(occur.taxa) %in% c("SiteId", "SampId"))])
    
    unique.sites <- unique(sites)
    siteIND <- match(sites, unique.sites)
    
    #if center is true substract the mean of each predictor, check if its divded by sd, I added the division by sd
    mean.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
        mean(k, na.rm = TRUE)
    })
    sd.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
        sd(k, na.rm = TRUE)
    })
    for(i in 1:length(env.cond[!(colnames(env.cond) %in% c("SiteId", "SampId"))])){
        #i = 6
        #message(i)
        env.cond[i+2] <- as.matrix(env.cond[i+2]) - mean.env.cond[i]
    }
    for(i in 1:length(env.cond[!(colnames(env.cond) %in% c("SiteId", "SampId"))])){
        #i = 6
        env.cond[i+2] <- as.matrix(env.cond[i+2]) / sd.env.cond[i]
    }
    
    if(cv == F){
        return(list(env.cond, occur.taxa))
    }else{
        
        test.y <- split[[2]]
        # join the environmental conditions to the occurrence data
        test.predictors <- left_join(test.y[, c("SiteId", "SampId")], env.cond.orig, by = c("SiteId", "SampId"))
        
        # Drop rows in y.test and predictors.test where the latter has any NAs
        ind <- !apply(is.na(test.predictors[,inf.fact]),1,FUN=any)
        ind <- ifelse(is.na(ind),FALSE,ind)
        test.predictors <- test.predictors[ind, ]
        test.y <- test.y[ind, ]
        
        for(i in 1:length(test.predictors[!(colnames(test.predictors) %in% c("SiteId", "SampId"))])){
            #i = 6
            #message(i)
            test.predictors[i+2] <- as.matrix(test.predictors[i+2]) - mean.env.cond[i]
        }
        for(i in 1:length(test.predictors[!(colnames(test.predictors) %in% c("SiteId", "SampId"))])){
            #i = 6
            test.predictors[i+2] <- as.matrix(test.predictors[i+2] / sd.env.cond[i])
        }
        
        
        
        # Keep the test sites
        test.sites <- test.predictors$SiteId
        test.samples <- test.predictors$SampId
        names(test.sites) <- test.samples
        return(list(env.cond, occur.taxa, test.predictors, test.y))
    }
}