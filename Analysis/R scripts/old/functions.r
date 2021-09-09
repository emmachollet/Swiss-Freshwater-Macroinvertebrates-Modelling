## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Functions needed to sort factors, construct data sets, ?? ----
## 
## --- "Bridging gap in macroinvertebrates community assembly" -- Project ---
## 
## --- June 21, 2021 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ---- Function to construct environmental data set with important factors ----

# sortEnvFact <- function(data.env, rank.env){
#     
#     # check if colnames are the same
#     if ( length(setdiff(colnames(data.inv), colnames(rank.env))) != 0){
#         print(paste("Colnames of", file.env.data, "don't match with colnames of", file.rank.env), sep=" ")
#     } # should be 0
#     
#     # replace "_" by "." in colnames to be consistent (also with colnames of invertebrate data)
#     # we could also change upper or lower case letter but hard to be consistent
#     colnames(rank.env) <- gsub("_", ".", colnames(rank.env))
#     colnames(data.env) <- gsub("_", ".", colnames(data.env))
#     
#     # sort the columns in 4 categories: the sample/site information = 3,
#     # the environmental factors to, keep in priority = 2, keep to explore = 1, exclude = 0
#     info    <- colnames(rank.env)[which(rank.env[1,] == 3)]
#     prio    <- colnames(rank.env)[which(rank.env[1,] == 2)]
#     explo   <- colnames(rank.env)[which(rank.env[1,] == 1)]
#     excl    <- colnames(rank.env)[which(rank.env[1,] == 0)]
#     
#     data.info <- data.env[, info]
#     data.env <- data.env[, c("SiteId", "SampId", prio, explo)]
#     
#     return(data.env)
#     # return(data.info, data.env, info, prio, explo, excl)
# }

# sortEnvFact(data.env, rank.env)

## ---- Construct ALL and BDM datasets ----

# Replace "0" and "1" by "absent" and "present" and convert them to factors
cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))

for (i in cind.taxa ) {
    data.inv[which(data.inv[,i] == 0),i] <- "absent"
    data.inv[which(data.inv[,i] == 1),i] <- "present"
    data.inv[,i] = as.factor(data.inv[,i])
}

# Construct main data set with site info, occ taxa and env fact (after FLOZ)
cind <- which(colnames(data.env) == "FLOZ")-1
data <- data.env[ , -(3:cind)]   # remove col in data.env that are between SiteId and "FLOZ"

data <- data.inv %>%
    left_join(data, by = c("SiteId", "SampId"))
dim(data)

# Construct BDM data set
data.BDM <- data %>%
    filter(MonitoringProgram == "BDM") # filter the rows, for columns use select
    # data[which(data.env$MonitoringProgram =="BDM") , ]
dim(data.BDM)

# BDM SampId and row indices for additional checks
BDM.samples <- data[which(data.env$MonitoringProgram =="BDM") , "SampId"]
length(BDM.samples) # should have 886 samples
BDMind <- which(data$SampId %in% BDM.samples)

# save("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts/AllDataSets_06_22.RData")

## WORKSPACE SAVED AS AllDataSets_06_22.RData ####

# rm(list=ls())
# load("Q:/Abteilungsprojekte/siam/Emma Chollet/Data processing/Swiss Freshwater Macroinvertebrates Modelling/Analysis/R scripts/AllDataSets_06_22.RData")

## ---- Construct dataset for modelling ----

# Select environmental factors and taxa
env.fact <- c("IAR", "saprobic_cond", "urban_area", "velocity", "BREITENVAR", "temperature") # inspired from Bogdan's selection in 2020
occ.taxa <- c("Occurrence.Hydraenidae") # , "Occurrence.Gammaridae", "Occurrence.Erpobdellidae") # chosen depending on mixed prevalence

# Summary of prevalence of chosen taxa
for ( i in 1:length(occ.taxa)){
    cat("Summary of absence and presence for", occ.taxa[i], ":", summary(data.BDM[, occ.taxa[i]]), "\n")
}

# Data set for modelling
omit.NA <- TRUE
if (omit.NA == TRUE){
    taxa.env <- na.omit(data.BDM[, c(occ.taxa, env.fact)]) # remove NA to avoid problems
} else {
    taxa.env  <- data.BDM[, c(occ.taxa, env.fact)]
}
str(taxa.env) # structure of the dataframe
head(taxa.env[,c(occ.taxa,env.fact)]) # see top 6 rows and 9 columns


## ---- Create the training and test datasets ----

# Set a seed in R is used for:
# 1. Reproducing the same output of simulation studies.
# 2. Help to debug the code when dealing with pseudorandom numbers.
set.seed(100) # !! Didn't understand well what it does, we do it each time we train a model so it means for each training we generate random numbers ? ####


# Step 1: Get row numbers for the training data
trainRowNumbers <- createDataPartition(taxa.env[,occ.taxa], p=0.8, list=FALSE) # split into training 80% and test 20% datasets

# Step 2: Create the training  dataset
trainData <- taxa.env[trainRowNumbers,]

# Step 3: Create the test dataset
testData <- taxa.env[-trainRowNumbers,]

# Store X and Y for later use.
x = trainData[, env.fact]
y = trainData[, occ.taxa]


## ---- Explore descriptive statistics ----

skimmed <- skim(trainData) # produce dataframe containing descriptive statistics for each column
# skimmed[, c(1:5, 9:11, 13:15)]


## ---- Impute missing values ----

# The dataset has few missing values across all columns, we may to do well to impute it.
# Impute, means to fill it up with some meaningful values.
# If the feature is a continuous variable, it is a common practice to replace the missing values with the mean of the column.
# And if it's a categorical variable, replace the missings with the most frequently occurring value, aka, the mode.
# This is quite a basic and a rather rudimentary approach.
# Instead we can actually predict the missing values by considering the rest of the available variables as predictors.
# A popular algorithm to do imputation is the k-Nearest Neighbors.

# Create the knn imputation model on the training data
preProcess_missingdata_model <- preProcess(trainData, method='knnImpute')
preProcess_missingdata_model

# It has centered (subtract by mean) 6 variables, ignored 3, used k=5 (considered 5 nearest neighbors) to predict the missing values
# and finally scaled (divide by standard deviation) 6 variables.

# Use the imputation model to predict the values of missing data points
trainData <- predict(preProcess_missingdata_model, newdata = trainData)
anyNA(trainData) # all missing values are successfully imputed


## ---- One-Hot Encoding ----

# If we have a categorical column as one of the features, it needs to be converted to numeric
# in order for it to be used by the machine learning algorithms.
# Just replacing the categories with a number may not be meaningful especially if there is no intrinsic ordering amongst the categories.
# So we convert the categorical variable with as many binary (1 or 0) variables as there are categories.

## NOT NEEDED AS ANY FEATURE (env predictor) IS CATEGORICAL/AS FACTORS

# # Creating dummy variables is converting a categorical variable to as many binary variables as here are categories.
# dummies_model <- dummyVars(Occurrence.Hydraenidae ~ ., data=trainData)
#
# # Create the dummy variables using predict. The Y variable (Occurrence.Hydraenidae) will not be present in trainData_mat.
# trainData_mat <- predict(dummies_model, newdata = trainData)
#
# # # Convert to dataframe
# trainData <- data.frame(trainData_mat)
#
# # # See the structure of the new dataset
# str(trainData)


## ---- Preprocess the data ----

# A lot of preprocessing types are available in caret (range, center, scale, pca, ...)


## ---- Visualize the importance of variables ----

featurePlot(x = trainData[, env.fact],
            y = trainData[, occ.taxa],
            plot = "box",   # to visualize two boxplots, each for absent and present
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"),
                          y = list(relation="free")))

featurePlot(x = trainData[, env.fact],
            y = trainData[, occ.taxa],
            plot = "density", # to visualize density of absence and presence in relation to env predictors
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"),
                          y = list(relation="free")))


## ---- Feature selection using recursive feature elimination (rfe) ----

# We might need a rigorous way to determine the important variables first before feeding them to the ML algorithm.
# RFE works in 3 broad steps:
# Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
# Step 2: Keeping priority to the most important variables, iterate through by building models of given subset sizes, that is,
#         subgroups of most important predictors determined from step 1. Ranking of the predictors is recalculated in each iteration.
# Step 3: The model performances are compared across different subset sizes to arrive at the optimal number and list of final predictors.


set.seed(100)
options(warn=-1)


ctrl <- rfeControl(functions = rfFuncs, # based on random forest algorithm
                   method = "repeatedcv", # method of cross-validation
                   repeats = 5,
                   verbose = FALSE)

lmProfile <- rfe(x=trainData[, env.fact], y=trainData[, occ.taxa], # receives the output of the rfeControl() as values
                 sizes = c(1:3, 5, 6), # what all model sizes (the number of most important features) the rfe should consider.
                 rfeControl = ctrl)

lmProfile


## ---- Train the model ----

# Here we will train a Multivariate Adaptive Regression Splines (MARS), a form of regression analysis that models nonlinearities.

modelLookup('earth') # describe the MARS algorithm (called 'earth' in R to avoid copyright conflicts)

# Set the seed for reproducibility
set.seed(100)

# ..
f <- reformulate(".", occ.taxa)

# Train the model using MARS and predict on the training data itself
model_mars = train(f, data=trainData, method='earth')
fitted <- predict(model_mars)

# Why is train() better than using the algorithm's fct directly ?
# Because train() does other things like: (see later in the tutorial)
# 1. Cross validating the model
# 2. Tune the hyper parameters for optimal model performance
# 3. Choose the optimal model based on a given evaluation metric
# 4. Preprocess the predictors (what we did so far using preProcess())

model_mars # look at information of output (method, nb of predictors, classes, pre-processing, hyperparameter tuning, ...)
plot(model_mars, main="Model Accuracies with MARS") # show how the various iterations of hyperparameter search performed


## ---- Compute variable importance ----

varimp_mars <- varImp(model_mars)
plot(varimp_mars, main="Variable Importance with MARS")


## ---- Prepare the test dataset ----

# We preprocess and transform the data as we did for training dataset.
# In the tutorial they followed a sequence:
# Missing Value imputation -> One-Hot Encoding -> Range Normalization
# But here we just impute missing value.

# Step 1: Impute missing values
testData2 <- predict(preProcess_missingdata_model, testData)

# # Step 2: Create one-hot encodings (dummy variables)
# testData3 <- predict(dummies_model, testData2)
#
# # Step 3: Transform the features to range between 0 and 1
# testData4 <- predict(preProcess_range_model, testData3)

# View
head(testData2[, 1:7])


## ---- Predict on test datset ----

# Predict on testData
predicted <- predict(model_mars, testData2)
head(predicted)


## ---- Confusion matrix ----

# Compute the confusion matrix
confusionMatrix(reference = testData[,occ.taxa], data = predicted, mode='everything', positive='present')


## ---- Control cv and summaries of results before training the model ----

# The train() function takes a trControl argument that accepts the output of trainControl().
# The latest control:
# 1. The cross validation method used (bootstrap sampling, k-fold cv, leave one/group out, ...)
# 2. How the results should be summarised (twoClassSummary if Y is binary class or multiClassSummary if the Y has more than 2 categories)

# Define the training control
fitControl <- trainControl(
    method = 'cv',                   # k-fold cross validation
    number = 5,                      # number of folds
    savePredictions = 'final',       # saves predictions for optimal tuning parameter
    classProbs = T,                  # should class probabilities be returned
    summaryFunction=twoClassSummary  # results summary function
)

## ---- Hyperparameter tuning ----

# There are two main ways to do hyper parameter tuning using the train():
# 1. Set the tuneLength (number of unique values for the tuning parameters caret will consider
#   while forming the hyper parameter combinations, Caret will automatically determine the values each parameter should take)
# 2. Define and set the tuneGrid (explicitly control what values should be considered for each parameter)

## Using tuneLength

# Step 1: Tune hyper parameters by setting tuneLength
set.seed(100)
model_mars2 = train(Occurrence.Hydraenidae ~ ., data=trainData, method='earth', tuneLength = 5, metric='ROC', trControl = fitControl)
model_mars2

# Step 2: Predict on testData and Compute the confusion matrix
predicted2 <- predict(model_mars2, testData2)
confusionMatrix(reference = testData[, occ.taxa], data = predicted2, mode='everything', positive='present')

## Using tuneGrid

# Step 1: Define the tuneGrid
marsGrid <-  expand.grid(nprune = c(2, 4, 6, 8, 10),
                         degree = c(1, 2, 3))

# Step 2: Tune hyper parameters by setting tuneGrid
set.seed(100)
model_mars3 = train(Occurrence.Hydraenidae ~ ., data=trainData, method='earth', metric='ROC', tuneGrid = marsGrid, trControl = fitControl)
model_mars3

# Step 3: Predict on testData and Compute the confusion matrix
predicted3 <- predict(model_mars3, testData2)
confusionMatrix(reference = testData[, occ.taxa], data = predicted3, mode='everything', positive='present')


## ---- Evaluate performance of multiple ML algorithms -----

# The resamples() function can be provided by multiple machine learning models and will collectively evaluate them.
# We first need to train multiple algorithms.

set.seed(100)

# Train the model using Adaboost
model_adaboost = train(Occurrence.Hydraenidae ~ ., data=trainData, method='adaboost', tuneLength=2, trControl = fitControl)
model_adaboost

set.seed(100)

# Train the model using Random Forest
model_rf = train(Occurrence.Hydraenidae ~ ., data=trainData, method='rf', tuneLength=5, trControl = fitControl)
model_rf

# set.seed(100)
#
# # Train the model using xgBoost Dart # takes a lot of time
# model_xgbDART = train(Occurrence.Hydraenidae ~ ., data=trainData, method='xgbDART', tuneLength=5, trControl = fitControl, verbose=F)
# model_xgbDART

set.seed(100)

# Train the model using SVM
model_svmRadial = train(Occurrence.Hydraenidae ~ ., data=trainData, method='svmRadial', tuneLength=15, trControl = fitControl)
model_svmRadial

# Compare model performances using resample()
models_compare <- resamples(list(ADABOOST=model_adaboost, RF=model_rf, # XGBDART=model_xgbDART,
                                 MARS=model_mars3, SVM=model_svmRadial))

# Summary of the models performances
summary(models_compare)

# Draw box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales) # see how algorithms perform in terms of ROC, Sepcificity and Sensitivity


## ---- Compare prediction performance ----

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv",
                             number=10,
                             repeats=3,
                             savePredictions=TRUE,
                             classProbs=TRUE)

algorithmList <- c('rf', 'adaboost', 'earth', # 'xgbDART',
                   'svmRadial')

set.seed(100)
models <- caretList(Occurrence.Hydraenidae ~ ., data=trainData, trControl=trainControl, methodList=algorithmList)
results <- resamples(models)
summary(results)

# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)


# to print in a pdf try
print(bwplot())

## ---- Ensemble prediction ----

# We can do it with caret, but I won't for now :)



