

exit <- function() { invokeRestart("abort") }  



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
    
    # likeli[which(likeli < 0.01)] <- 0.01
    
    st.dev <- -2 * sum(log(likeli)) / no.obs
    names(st.dev) <- "StandardizedDeviance"
    return(st.dev)
}


# Function to select the best performance of stand dev (i.e the lowest one)
lowest <- function (x, metric, maximize = F){
    
    best <- which.min(x[, metric])
    return(best)
}

apply.ml.model <- function(splitted.data, list.taxa, env.fact, algorithm, selec.metric = "StandardizedDeviance", train.control = trainControl(), ...){
    
    data.train <- splitted.data[["Training data"]]
    data.test <- splitted.data[["Testing data"]]
    
    # Make a list with the outputs of the algorithm for each taxon in list.taxa
    list.outputs <- vector(mode = 'list', length = length(list.taxa))
    names(list.outputs) <- list.taxa
    
    if(nrow(data.test)>0){which.set <- c("training set", "testing set")
    } else {which.set <- c("training set")}
    out <- c("Observation", #1
             "Prediction factors", #2 
             "Prediction probabilities", #3 
             "Likelihood", #4
             "Performance", #5
             "Confusion matrix") #6
    output.names <- c("Trained model", "Variable importance", c(outer(out, which.set, FUN = paste)))
    
    
    for (j in 1:length(list.taxa)){
        
        print(paste("Applying", algorithm, "to", list.taxa[j]))
        
        # For each taxon we do a list with the important outputs
        # temp.output <- c("Trained model",                                 #1
        #                  "Variable importance",                           #2
        #                  "Training performance",                          #2
        #                  "Prediction on testing set (presence/absence)",  #4
        #                  "Prediction on testing set (probabilities)",     #5
        #                  "Observation testing set",                       #6
        #                  "Likelihood testing set",                        #7
        #                  "Testing performance",                           #8
        #                  "Confusion matrix")                              #9
        temp.list <- vector(mode = 'list', length = length(output.names))
        names(temp.list) <- output.names
        
        temp.train <- na.omit(data.train[, c("SiteId", "SampId", "X", "Y", list.taxa[j], env.fact)]) # create a temporary training dataset with the taxon and env fact, to 
        temp.test <- na.omit(data.test[, c("SiteId", "SampId", "X", "Y", list.taxa[j], env.fact)])
        temp.sets <- list(temp.train, temp.test)
        f <- reformulate(env.fact, list.taxa[j]) # not obligatory, just if we want to express it as a fct
        
        # Why is train() better than using the algorithm's fct directly ?
        # Because train() does other things like:
        # 1. Cross validating the model
        # 2. Tune the hyper parameters for optimal model performance
        # 3. Choose the optimal model based on a given evaluation metric
        # 4. Preprocess the predictors (what we did so far using preProcess())
        
        # model <- train(x = temp.train[,env.fact], y = temp.train[,list.taxa[j]], method = algorithm, trControl = train.control)
        model <- train(f, data = temp.train, metric = selec.metric, method = algorithm, trControl = train.control)
        
        temp.list[["Trained model"]] <- model
        temp.list[["Variable importance"]] <- varImp(model)
        
        for(n in 1:length(which.set)){
            # Observation
            temp.list[[paste(out[1],which.set[n])]] <- temp.sets[[n]]
            # Prediction factors
            temp.list[[paste(out[2],which.set[n])]] <- predict(model, temp.sets[[n]])
            # Prediction probabilities
            temp.list[[paste(out[3],which.set[n])]] <- predict(model, temp.sets[[n]], type = 'prob')
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
            temp.list[[paste(out[6],which.set[n])]] <- confusionMatrix(reference = temp.sets[[n]][,list.taxa[j]], 
                                                               data = temp.list[[paste(out[2],which.set[n])]], mode='everything', positive='present')
            
        }
            # # temp.list[[2]] <- NA
            # temp.list[["Training performance"]] <- min(model$results[,selec.metric], na.rm = T)
            #     # model$results[which(model$results[,selec.metric] == model$bestTune[,1]), selec.metric]
            # # ?? not sure its the right final standerdized deviance of the model ####
            # 
            # if (nrow(temp.test)>0){ # if there is a testing set we predict on it, else all predictions are empty
            # 
            #     temp.list[["Prediction on testing set (presence/absence)"]] <- predict(model, temp.test)
            #     temp.list[["Prediction on testing set (probabilities)"]] <- predict(model, temp.test, type = 'prob')
            #     temp.list[["Observation testing set"]] <- temp.test
            #     
            #     likeli <- 1:nrow(temp.test)
            #     for(i in 1:nrow(temp.test)){
            #         if(temp.test[i,list.taxa[j]] == "present"){
            #             likeli[i] <- temp.list[["Prediction on testing set (probabilities)"]][i, "present"]
            #         } else if (temp.test[i,list.taxa[j]] == "absent" ){
            #             likeli[i] <- temp.list[["Prediction on testing set (probabilities)"]][i, "absent"]
            #         }
            #     }
            #     
            #     temp.list[["Likelihood testing set"]] <- likeli
            #     temp.list[["Testing performance"]] <- -2 * sum(log(likeli)) / nrow(temp.test)
            #     temp.list[["Confusion matrix"]] <- confusionMatrix(reference = temp.test[,list.taxa[j]], 
            #                                       data = temp.list[["Prediction on testing set (presence/absence)"]], mode='everything', positive='present')
            # }
        list.outputs[[j]] <- temp.list
    }
    
    return(list.outputs)
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
        
        no.pres <- sum(data[, list.taxa[j]] == "present", na.rm = TRUE)
        no.abs  <- sum(data[, list.taxa[j]] == "absent", na.rm = TRUE)
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

