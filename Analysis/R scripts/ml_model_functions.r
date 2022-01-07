

# Split data in training and testing datasets
split.data <- function(data, training.ratio, variable = "random", bottom = T){
    
    # add an argument in the fct, eg splitcol, which would be here SiteId
    # first split by this argument, then select rows (aka samples) according to their SiteId
    
    if( variable == "random"){
        sample.size <- floor(ratio*nrow(data)) # size of training set is 80% of the whole dataset
        # set.seed(100) # used for reproducing the same output of simulations
        set.seed(77) # try another one
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
    
    splitted.data <- list("Training data" = data.train, "Testing data" = data.test)

    return(splitted.data)
}


# Define metric function for trainControl within caret package
stand.dev <- function(data, lev = c("present","absent"), model = NULL){
    
    no.obs <- dim(data)[1]
    
    likeli <- 1:no.obs
    
    likeli[which(data$obs == lev[1])] <- data[which(data$obs == lev[1]), lev[1]]
    likeli[which(data$obs == lev[2])] <- data[which(data$obs == lev[2]), lev[2]]
    
    likeli[which(likeli < 0.01)] <- 0.01
    
    st.dev <- -2 * sum(log(likeli)) / no.obs
    names(st.dev) <- "StandardizedDeviance"
    return(st.dev)
}


# Function to select the best performance of stand dev (i.e the lowest one)
lowest <- function (x, metric, maximize = F){
    
    best <- which.min(x[, metric])
    return(best)
}

# Function to apply ML algorithms
apply.ml.model <- function(splitted.data, list.algo, list.taxa, env.fact, selec.metric = "StandardizedDeviance", CV = T, prev.inv, ...){
    
    data.train <- splitted.data[["Training data"]]
    data.test <- splitted.data[["Testing data"]]
    
    # Adapt list taxa to taxa actually present in data.train
    list.taxa <- list.taxa[list.taxa %in% colnames(data.train)]
    
    # Make a list to store the outputs of each model
    outputs <- vector(mode = 'list', length = no.algo)
    names(outputs) <- list.algo
    
    for(k in 1:no.algo){
    
        algorithm = list.algo[k]
        
        # Make a list with the outputs of the algorithm for each taxon in list.taxa
        list.outputs <- vector(mode = 'list', length = length(list.taxa))
        names(list.outputs) <- list.taxa
        
        # if testing data set exists, create outputs for results on training and testing data sets
        # else only list outputs for training data set
        if(CV == T){which.set <- c("training set", "testing set")
        } else {which.set <- c("training set")}
        out <- c("Observation", #1
                 "Prediction factors", #2 
                 "Prediction probabilities", #3 
                 "Likelihood", #4
                 "Performance")#, #5
                 # "Confusion matrix") #6
        output.names <- c("Trained model",
                          # "Variable importance",
                          c(outer(out, which.set, FUN = paste)))
        
        for (j in 1:length(list.taxa)){
            
            temp.list <- vector(mode = 'list', length = length(output.names))
            names(temp.list) <- output.names
            
            temp.train <- na.omit(data.train[, c("SiteId", "SampId",
                                                 "X", "Y", 
                                                 list.taxa[j], env.fact)]) # create a temporary training dataset with the taxon and env fact, to 
            if(CV == T){temp.test <- na.omit(data.test[, c("SiteId", "SampId",
                                                   "X", "Y", 
                                                   list.taxa[j], env.fact)])
                        temp.sets <- list(temp.train, temp.test)
            } else {    temp.sets <- list(temp.train)     } 
            
            f <- reformulate(env.fact, list.taxa[j]) # write formula (target variable ~ explanatory variables) to apply the model
            
            # Why is train() better than using the algorithm's fct directly ?
            # Because train() does other things like:
            # 1. Cross validating the model
            # 2. Tune the hyper parameters for optimal model performance
            # 3. Choose the optimal model based on a given evaluation metric
            # 4. Preprocess the predictors (what we did so far using preProcess())
            
            set.seed(2021) # for reproducibility of the folds
            
            folds <- groupKFold(temp.train$SiteId, 3) # create 3 folds, grouped by SiteId
            bal1 <- summary(temp.train[folds$Fold1, list.taxa[j]])
            bal2 <- summary(temp.train[folds$Fold2, list.taxa[j]])
            bal3 <- summary(temp.train[folds$Fold3, list.taxa[j]])
            bal.all <- bind_rows(bal1, bal2, bal3, .id = "Fold")
            
            any(bal.all == 0)
            if( any(bal.all == 0) & grepl("svm", algorithm) ){
            
                cat("During training, classes are too unbalanced for the algorithm", algorithm,
                "\nin fold 1:", summary(temp.train[folds$Fold1, list.taxa[j]]),
                "\nin fold 2:", summary(temp.train[folds$Fold2, list.taxa[j]]),
                "\nin fold 3:", summary(temp.train[folds$Fold3, list.taxa[j]]),
                "\nWe apply Null model instead.", "\n")
                
                model <- "NULL_MODEL"
                prev.taxa <- prev.inv[which(prev.inv[,"Occurrence.taxa"] == list.taxa[j]), "Prevalence"]
                
            } else {
                
                cat(paste("\nApplying", algorithm, "to", j, list.taxa[j], "\n"))
                
                # folds1 <- create_folds(temp.train[,list.taxa[j]], k =3) # create stratified folds
                # folds2 <- groupKFold(temp.train$SiteId, 2) # try to have less folds but balanced better
    
                # Define the training control
                train.control <- trainControl(
                    method = 'cv',                   # k-fold cross validation
                    number = 3,                      # number of folds
                    index = folds,                   # provide indices computed with groupKFold for the k-fold CV
                    # repeats = 1,                   # for repeated k-fold cross-validation 'repeatedcv' only: the number of complete sets of folds to compute
                    # savePredictions = 'final',     # saves predictions for optimal tuning parameter
                    classProbs = T,                  # should class probabilities be returned
                    #summaryFunction = twoClassSummary, # default metric function choosen (accuracy, ROC)
                    summaryFunction = stand.dev ,
                    selectionFunction = lowest       # we want to minimize the metric
                )
                
                # Tried another training control to avoid unbalenced classes
                # train.control2 <- trainControl(
                #     method = 'boot',                   # k-fold cross validation
                #     number = 3,                      # number of folds
                #     classProbs = T,                  # should class probabilities be returned
                #     summaryFunction = stand.dev ,
                #     selectionFunction = lowest       # we want to minimize the metric
                # ) 
                # model <- train(x = temp.train[,env.fact], y = temp.train[,list.taxa[j]], method = algorithm, trControl = train.control) # alternative to formula writing
                model <- train(f, data = temp.train, metric = selec.metric, method = algorithm, trControl = train.control)
                # model2 <- train(f, data = temp.train, metric = selec.metric, method = algorithm, trControl = train.control2) # Tried with the alternative training control
            
            temp.list[["Trained model"]] <- model
            # temp.list[["Variable importance"]] <- varImp(model)
            
            for(n in 1:length(which.set)){
                # Observation
                temp.list[[paste(out[1],which.set[n])]] <- temp.sets[[n]]
                n.obs <- dim(temp.sets[[n]])[1]
                # Prediction factors
                temp.list[[paste(out[2],which.set[n])]] <- if(model[1] == "NULL_MODEL"){
                    ifelse(prev.taxa > 0.5, rep("present", n.obs), rep("absent", n.obs))
                } else { predict(model, temp.sets[[n]]) }
                # Prediction probabilities
                temp.list[[paste(out[3],which.set[n])]] <- if(model[1] == "NULL_MODEL"){
                    data.frame("absent" = rep(1-prev.taxa, n.obs), "present" = rep(prev.taxa, n.obs))
                } else { predict(model, temp.sets[[n]], type = 'prob') }
                # Likelihood
                likeli <- 1:nrow(temp.sets[[n]])
                for(i in 1:nrow(temp.sets[[n]])){
                    if(temp.sets[[n]][i,list.taxa[j]] == "present"){
                        likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "present"]
                    } else if (temp.sets[[n]][i,list.taxa[j]] == "absent" ){
                        likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "absent"]
                    }
                }
                temp.list[[paste(out[4],which.set[n])]] <- likeli
                # Performance
                temp.list[[paste(out[5],which.set[n])]] <- -2 * sum(log(likeli)) / nrow(temp.sets[[n]])
                # Confusion matrix
                # temp.list[[paste(out[6],which.set[n])]] <- if(model[1] != "NULL_MODEL"){confusionMatrix(reference = temp.sets[[n]][,list.taxa[j]], 
                #                                                    data = temp.list[[paste(out[2],which.set[n])]], mode='everything', positive='present')}
                # 
            }
            }
            list.outputs[[j]] <- temp.list
        }
    
        outputs[[k]] <- list.outputs
    }
    
    return(outputs)
}
# 
# calc.stand.deviance <- function(likeli){
#     
#     return()
#     
# }

apply.null.model <- function(data, list.taxa, prev.inv){
    
    # Make a list with the "outputs" of null model
    list.outputs <- vector(mode = 'list', length = length(list.taxa))
    names(list.outputs) <- list.taxa
    
    for (j in 1:length(list.taxa)){
            
        temp.output <- c("Likelihood", "Performance")
        temp.list <- vector(mode = 'list', length = length(temp.output))
        names(temp.list) <- temp.output
        
        no.pres <- sum(data[, list.taxa[j]] == 1, na.rm = TRUE)
        no.abs  <- sum(data[, list.taxa[j]] == 0, na.rm = TRUE)
        no.obs  <- no.pres + no.abs
        prev    <- prev.inv[prev.inv$Occurrence.taxa == list.taxa[j],"Prevalence"]
        likeli <- rep(c(prev, 1-prev),c(no.pres,no.abs))
        temp.list[["Likelihood"]] <- likeli
        
        st.dev <- -2 * sum(log(likeli)) / no.obs
        temp.list[["Performance"]] <- st.dev
        
        list.outputs[[j]] <- temp.list
    }
    return(list.outputs)
}

# from here should be useless ####

stand.deviance.null.model <- function(list.taxa, prev.inv, data.inv){
    
    # Make a list with the deviance for each taxon in list.taxa
    list.deviance <- vector(mode = 'list', length = length(list.taxa))
    names(list.deviance) <- list.taxa
    
    for (j in 1:length(list.taxa)){
        
        no.pres <- sum(data.inv[, list.taxa[j]] == "present", na.rm = TRUE)
        no.abs  <- sum(data.inv[, list.taxa[j]] == "absent", na.rm = TRUE)
        no.obs  <- no.pres + no.abs
        prev    <- prev.inv[prev.inv$Occurrence.taxa == list.taxa[j],"Prevalence"]
        likeli <- rep(c(prev, 1-prev),c(no.pres,no.abs))
        
        
        d.taxon <- -2 * sum(log(likeli)) / no.obs
        
        list.deviance[[j]] <- d.taxon
        
        }
    
    return(list.deviance)
}

stand.deviance.ml.model <- function(list.taxa, output){
    # Make a list with the deviance for each taxon in list.taxa
    list.deviance <- vector(mode = 'list', length = length(list.taxa))
    names(list.deviance) <- list.taxa
    
    for (j in 1:length(list.taxa)){
        
        temp.df <- data.frame(output[[list.taxa[j]]][["Observation testing set"]][list.taxa[j]], stringsAsFactors = FALSE)
        temp.df$likeli <- 1:nrow(temp.df)
        
        for(i in 1:nrow(temp.df)){
            if(temp.df[i,list.taxa[j]] == "present"){
                temp.df[i,"likeli"] <- output[[list.taxa[j]]][["Prediction on testing set (probabilities)"]][i, "present"]
            } else if (temp.df[i,list.taxa[j]] == "absent" ){
                temp.df[i, "likeli"] <- output[[list.taxa[j]]][["Prediction on testing set (probabilities)"]][i, "absent"]
            }
        }
        
        likeli <- temp.df$likeli
        no.obs <- nrow(temp.df)
        
        d.taxon <- -2 * sum(log(likeli)) / no.obs
        
        list.deviance[[j]] <- d.taxon
        
    }
    
    return(list.deviance)
}

pred.perf <- function(list.taxa, dev.null, dev.ml){
    
    # Make a list with the different predictive performance for each taxon in list.taxa
    list.deviance <- vector(mode = 'list', length = length(list.taxa))
    names(list.deviance) <- list.taxa
    
    for (j in 1:length(list.taxa)){
        
        temp.pred.perf <- c("Deviance null model", "Deviance ml algorithm", "Explanatory power")
        temp.list <- vector(mode = 'list', length = length(temp.pred.perf))
        names(temp.list) <- temp.pred.perf
        
        temp.list[[1]] <- dev.null[[list.taxa[j]]]
        temp.list[[2]] <- dev.ml[[list.taxa[j]]]
        temp.list[[3]] <- (dev.null[[list.taxa[j]]] - dev.ml[[list.taxa[j]]]) / dev.null[[list.taxa[j]]]
        
        list.deviance[[j]] <- temp.list
    }
    
    return(list.deviance)
}

