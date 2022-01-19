# youpi

## ---- Produce pdf of a list of plots ----

print.pdf.plots <- function(list.plots, width = 12, height = width*3/4, dir.output, info.file.name = "", file.name){
  
  pdf(paste0(dir.output, info.file.name, file.name), paper = 'special', width = width, height = height, onefile = TRUE)
  
  for (n in 1:length(list.plots)) {
    
    plot(list.plots[[n]])
    
  }
  
  dev.off()
  
}

## ---- Inputs for Swiss map plot ----

map.inputs <- function(dir.env.data, data.env){
    
    # try to understand what Bogdan did
    inputs = list()
    
    # Obtain simple feature for borders of Switzerland, and major lakes and rivers
    inputs$ch <- st_read(paste0(dir.env.data,"workspace.gdb"), layer="switzerland", stringsAsFactors = F)
    inputs$ch <- filter(inputs$ch, NAME=="Schweiz")
    
    inputs$rivers.major <- st_read(paste0(dir.env.data,"workspace.gdb"), layer = "major_rivers", stringsAsFactors = F)
    inputs$lakes.major <- st_read(paste0(dir.env.data,"workspace.gdb"), layer = "major_lakes", stringsAsFactors = F)
    
    # Get coordinates of all sites 
    inputs$xy <- select(data.env, SiteId, X, Y)
    inputs$xy <- distinct(inputs$xy)
    
    return(inputs)
}

## ---- Explorative plots -----

map.env.fact <- function(inputs, env.fact, data.env, env.explan, dir.output, file.prefix, file.name=NA){
    
    plot.title <- paste(gsub("_", "", file.prefix), "environmental dataset")
    
    if(!is.na(file.name)){
        pdf(paste0(dir.output, file.prefix, file.name), paper = 'special', 
            width = 11, 
            onefile = TRUE)
    }
    
    for (k in 1:length(env.fact)){
        
        variable <- env.fact[k]
        temp.data.env <- data.env[, c("X","Y", env.fact[k])]
        no.rows <- nrow(temp.data.env)
        no.na <- sum(is.na(temp.data.env[,variable]))
        explanation <- env.explan[which(env.explan$column.name == variable), "explanation"]
        
        g <- ggplot()
        g <- g + geom_sf(data = inputs$ch, fill="#E8E8E8", color="black")
        g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
        g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
        g <- g + geom_point(data = temp.data.env, aes(x=X, y=Y, color= temp.data.env[, variable]), size= 3, alpha= 0.6)
        
        if( is.numeric(temp.data.env[,variable]) == TRUE ){
            # Set up scales
            k.min <- round(min(temp.data.env[,variable], na.rm = T), 1)
            k.max <- round(max(temp.data.env[,variable], na.rm = T), 1)
            k.int <- (k.max - k.min)/5 # ; if.int <- round(if.int)
            
            g <- g + scale_colour_gradient2(name = variable, 
                                            #low = "coral", # useless
                                            high = "firebrick3",
                                            space = "Lab", na.value = "grey20", guide = "colourbar")
            # g <- g + guides(color = guide_bins(axis = T, show.limits = T, override.aes = list(alpha = 1)))
        }
        
        g <- g + theme_void(base_size = 18)
        g <- g + theme(# plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(colour="transparent"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
        g <- g + labs(title = paste("Geographic distribution of",variable),
                      subtitle = paste0(no.na, " NAs out of ", no.rows, " samples \n", explanation), colour = variable)
        
        print(g)
        
        p <- ggplot(data = temp.data.env, aes(temp.data.env[, variable]))
        
        if( is.numeric(temp.data.env[,variable]) == TRUE ){
            
            q <- p + geom_histogram(col="grey20", fill="firebrick3", alpha = 0.2,
                                    binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)))
            q <- q + labs(title=paste("Histogram for", variable), x=variable, y="Frequency")
            
            r <- p + geom_density(fill="firebrick3", alpha = 0.2)
            r <- r + labs(x=variable, col="grey20", y="Density")
            
            p <- ggarrange(q,r, ncol = 2)
            
        } else {
            p <- p + geom_bar(col="grey20", fill="firebrick3", alpha = 0.2)
            p <- p + labs(title=paste("Histogram for", variable))
        }
        
        print(p)
        # print(ggarrange(g,p, nrow =2))
        
    }
    if(!is.na(file.name)){
        dev.off()
    }
    
}

plot.data.envvstax <- function(data, env.fact, list.taxa){

  no.taxa <- length(list.taxa)
  list.plots <- list()
  n = 1
  
  for (j in 1:no.taxa){
    
    plot.data <- data[,c(list.taxa[j],env.fact)]
    plot.data <- na.omit(plot.data)
    
    list.plots[[n]] <- featurePlot(x = plot.data[,env.fact],
                                    y = plot.data[,list.taxa[j]],
                                    plot = "box",   # to visualize two boxplots, each for absent and present
                                    strip=strip.custom(par.strip.text=list(cex=.7)),
                                    scales = list(x = list(relation="free"),
                                                  y = list(relation="free")),
                                    main = list.taxa[j])
    # n <- n + 1
    # list.plots[[n]] <- featurePlot(x = plot.data[,env.fact],
    #                                 y = plot.data[,list.taxa[j]],
    #                                 plot = "density", # to visualize density of absence and presence in relation to env predictors
    #                                 strip=strip.custom(par.strip.text=list(cex=.7)),
    #                                 scales = list(x = list(relation="free"),
    #                                               y = list(relation="free")),
    #                                 main = list.taxa[j])
    n <- n + 1
    
    # Plot pairs of env. fact. vs taxa, too heavy to print
    # list.plots[[n]] <- featurePlot(x = data[,env.fact],
    #                                 y = data[,list.taxa[j]],
    #                                 plot = "pairs",
    #                                 ## Add a key at the top
    #                                 auto.key = list(columns = 2))
    # n <- n + 1
    # to plot the last featurePlots but by pairs, really heavy

  }

  return(list.plots)
  
}


## ---- Models analysis plots ----

plot.df.perf <- function(df.perf, list.models, list.taxa, CV){
    
    size.val <- ifelse(length(list.taxa) < 25, 5, 1)
    title <- ifelse(CV, "Table of models predictive performance \n during Cross Validation",
                    "Table of models quality of fit \n during Calibration")
    
    if(CV){ cind <- 1:which(colnames(df.perf) == "Taxa")
    } else { cind <- 1:which(colnames(df.perf) == "Taxa")
      # cind <- 1:which(colnames(df.perf) %in% list.models | colnames(df.perf) == "Taxa") # ECR: was there for FIT
      }
    temp.df <- df.perf[,cind]
    # temp.df$models <- row.names(df.perf)
    melted.df <- melt(temp.df, id = "Taxa")
    p <- ggplot(data = melted.df, aes(x = Taxa, y = variable, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(midpoint = 1, low = "#007139", mid ="grey70", high = "#c2141b", 
                             limits = c(0, 2)) +
        scale_x_discrete(limits = list.taxa) +
        labs(title = title, 
             x = "", y = "", fill = "Standardized \n deviance") +
        theme(plot.title = element_text(hjust = 0.5, colour = "black"), 
              axis.title.x = element_text(face="bold", colour="darkgreen", size = 2),
              axis.text.x = element_text(angle=90),
              axis.title.y = element_text(face="bold", colour="darkgreen", size = 2),
              legend.title = element_text(face="bold", colour="black", size = 10))  +
        geom_text(aes(x = Taxa, y = variable, label = round(value, 2)),
                  color = "black", fontface = "bold", size = size.val)
    
    return(list(p))
    
}



# Compare models
model.comparison <- function(df.perf, list.models, CV){
    
    list.models.temp <- c("#000000" = "Null_model")
    list.models <- c(list.models.temp, list.models)
    
    # Make a vector of colors
    col.vect <- names(list.models)
    names(col.vect) <- list.models
    
    no.taxa <- nrow(df.perf)
    
    if(CV){
        title <- paste("Models comparison in predictive performance")
    } else {
        title <- paste("Models comparison in quality of fit")
    }
    title <- paste(title, "for", no.taxa, "taxa")
    
    plot.data <- df.perf
    plot.data$null.perf <- plot.data[,"Null_model"]
    
    plot.data <- gather(plot.data, key = model, value = performance, all_of(c("Null_model",list.models)))
    plot.data <- gather(plot.data, key = expl.pow_model, value = expl.pow, all_of(paste0("expl.pow_", list.models)))
    
    # Prevalence vs stand dev
    p1 <- ggplot()
    p1 <- p1  + geom_point(data = plot.data, aes(x = Prevalence, y = performance, 
                                                 colour = model, 
                                                 shape = plot.data[,"Taxonomic level"]), 
                           # alpha = 0.4,
                           size = 3)
    #p1 <- p1 + geom_line(data = plot.data, aes(x = Prevalence, y = null.perf), linetype = "dashed", alpha=0.4, show.legend = FALSE) # to plot null model as dash line between data points
    p1 <- p1 + stat_function(fun=function(x) -2*(x*log(x) + (1-x)*log(1-x))) # to plot null model as function line
    p1 <- p1  + labs(y = "Standardized deviance",
                     x = "Prevalence (%)",
                     shape = "Taxonomic level",
                     color = "Model",
                     title = title)
    p1 <- p1 + scale_colour_manual(values=col.vect)
    p1 <- p1 + theme_bw(base_size = 20)
    p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=6)))
    
    # Boxplots
    p2 <- ggplot(plot.data, aes(x=model, y = performance, fill = model), alpha = 0.4) 
    p2 <- p2 + geom_boxplot()
    p2 <- p2 + scale_fill_manual(values=col.vect)
    p2 <- p2 + labs(title = title)
    p2 <- p2 + scale_x_discrete(limits = rev(list.models))
    p2 <- p2 + coord_flip()
    p2 <- p2 + theme_bw(base_size = 20)
    p2 <- p2 + theme(legend.position = "none")
    p2 <- p2 + labs(x="Models",
                    y="Standardized deviance",
                    # fill = "Models",
                    title = title)
    # Boxplots
    p3 <- ggplot(plot.data, aes(x=reorder(model, -performance, na.rm = T), y = performance, fill = model), alpha = 0.4) 
    p3 <- p3 + geom_boxplot()
    p3 <- p3 + scale_fill_manual(values=col.vect)
    p3 <- p3 + labs(title = title)
    # p3 <- p3 + scale_x_discrete(limits = rev(list.models))
    p3 <- p3 + coord_flip()
    p3 <- p3 + theme_bw(base_size = 20)
    p3 <- p3 + theme(legend.position = "none")
    p3 <- p3 + labs(x="Models",
                    y="Standardized deviance",
                    # fill = "Models",
                    title = title,
                    subtitle = "Ordered by decreasing mean")

    return(list(p1,p2,p3))
    
}

plot.dl.perf <- function(df.perf.dl.comb, list.models){
  
  # make a vector of colors
  col.vect <- names(list.models)
  names(col.vect) <- list.models
  col.vect <- c("Null_model" = "#000000", col.vect)
  
  title <- paste("Models comparison in predictive performance with and without data leakage")
  
  plot.data <- df.perf.dl.comb
  plot.data$null.perf <- plot.data[,"Null_model"]
  
  plot.data <- gather(plot.data, key = model, value = performance, all_of(c("Null_model",list.models)))
  plot.data <- gather(plot.data, key = expl.pow_model, value = expl.pow, all_of(paste0("expl.pow_", list.models)))
  
  # Boxplots
  p3 <- ggplot(plot.data, aes(x=model, y = performance, fill = as.factor(DL)), alpha = 0.4) 
  p3 <- p3 + geom_boxplot() #facet_wrap(~ DL) 
  #p3 <- p3 + scale_fill_manual(values=col.vect)
  p3 <- p3 + labs(title = title)
  p3 <- p3 + scale_x_discrete(limits = rev(list.models))
  p3 <- p3 + theme_bw(base_size = 12)
  p3 <- p3 + labs(x="Models",
                  y="Standardized deviance",
                  title = title,
                  fill = "Data leakage")
  
  return(list(p3))
  
}



# Plot model performance against hyperparameters

plot.perf.hyperparam <- function(outputs, list.algo, list.taxa){
    
    list.algo <- list.algo[! list.algo %in% c("glm")] # glm doesn't have hyperparameters
    no.algo <- length(list.algo)
    no.taxa <- length(list.taxa)
    
    list.plots <- list()
    
    for(j in list.taxa){
        
        temp.list.plots <- vector(mode = 'list', length = no.algo) 
        names(temp.list.plots) <- list.algo 
        
        for (l in list.algo){ 
            
            # cat(paste("Computing performance of hyperparameters of", list.algo[l], "for", which(list.taxa == j), j))
            
            trained.mod <- outputs[[l]][[j]][["Trained model"]]
            # Show how the various iterations of hyperparameter search performed
            if(trained.mod != "NULL_MODEL"){temp.list.plots[[l]] <- plot(trained.mod, main = l)}
        }
        
        title <- paste("Performance for each hyperparameter for taxa", j)
        p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
        
        list.plots[[j]] <- p
    }
    return(list.plots)
}

#  Plot variable importance
plot.varimp <- function(outputs, list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  
  list.plots <- list()
  n = 1
  
  # ECR: Will try to fix the table later
  # # print directly table with variable importance for each algo and taxa (too complicated to put in a fct)
  # file.name <- "TableVarImp.pdf"
  # pdf(paste0(dir.plots.output, info.file.name, file.name), paper = 'special', width = 12, height = 9, onefile = TRUE)
  # temp.df <- data.frame(matrix(ncol = no.algo*no.taxa, nrow = no.env.fact))
  # colnames(temp.df) <- c(outer(list.algo, list.taxa, FUN = paste))
  # rownames(temp.df) <- env.fact
  # for (j in list.taxa) {
  #     for (l in list.algo) {
  #         for (k in 1:no.env.fact) {
  #             temp.df[env.fact[k],paste(l, j)] <- outputs[[l]][[j]][["Variable importance"]][["importance"]][env.fact[k],1]
  #         }
  #     }
  # }
  # temp.df$mean.imp <- rowMeans(temp.df)
  # temp.df <- as.matrix(temp.df)
  # par(mar=c(1,5,15,3)+ 0.2, xaxt = "n")
  # plot(temp.df, 
  #      #key = NULL,
  #      digits = 2, text.cell=list(cex=0.5),
  #      # axis.col=list(side=3, las=2), 
  #      axis.row = list(side=2, las=1),
  #      col = viridis,
  #      xlab = "",
  #      ylab = "",
  #      cex.axis = 0.5,
  #      srt = 45,
  #      main = "Variable importance for ML algorithm applied to taxa"
  # )
  # axis(1, at=seq(1:ncol(temp.df)+1), labels = FALSE)
  # text(seq(1:ncol(temp.df)+1), par("usr")[4] + 0.15, srt = 50, 
  #      labels = colnames(temp.df), adj= 0, cex = 0.5, xpd = T)
  # 
  # dev.off()
  
  for(j in list.taxa){
    
    temp.list.plots <- vector(mode = 'list', length = no.algo)
    names(temp.list.plots) <- list.algo  
    
    for (l in list.algo){

      # Plot variable importance
      perf <- round(outputs[[l]][[j]][["Performance training set"]], digits = 3)
      plot.title <- paste("Algo:", l, 
                          "\nPerformance on training set:", 
                          ifelse(length(perf) > 1, perf, perf[1] ))
      
      temp.list.plots[[l]] <- ggplot(outputs[[l]][[j]][["Variable importance"]]) +
        ggtitle(plot.title)
    }
    
    title <- j
    p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
      
    list.plots[[n]] <- p
    n <- n + 1
  }
  
  return(list.plots)
}

# plot.table.varimp <- function(outputs, list.algo, list.taxa, env.fact){
#   
#   temp.list <- list()
#   no.algo <- length(list.algo)
#   no.taxa <- length(list.taxa)
#   no.env.fact <- length(env.fact)
#   
#   temp.df <- data.frame(matrix(ncol = no.algo*no.taxa, nrow = no.env.fact))
#   colnames(temp.df) <- c(outer(list.algo, list.taxa, FUN = paste))
#   rownames(temp.df) <- env.fact
#   
#   for (j in 1:no.taxa) {
#     for (l in 1:no.algo) {
#       for (k in 1:no.env.fact) {
#         temp.df[env.fact[k],paste(list.algo[l], list.taxa[j])] <- outputs[[l]][[j]][["Variable importance"]][["importance"]][env.fact[k],1]
#       }
#     }
#   }
#   temp.df$mean.imp <- rowMeans(temp.df)
#   temp.df <- as.matrix(temp.df)
#   temp.list[[1]] <- plot(temp.df, 
#                            #key = NULL,
#                            digits = 2, text.cell=list(cex=0.5),
#                            axis.col=list(side=3, las=2), axis.row = list(side=2, las=1),
#                            col = viridis,
#                            xlab = "",
#                            ylab = "Environmental factor",
#                            main = "Variable importance for ML algorithm applied to taxa")
#  return(temp.list)
# }

# Plot ICE manually 

plot.ice.per.taxa <- function(taxa, outputs, list.algo, env.fact, normalization.data){
    
    no.algo <- length(list.algo)
    
    # Initialize final list with all arranged plots
    list.plots <- list()
    
    for(k in env.fact){
        
        # Make temporary list of ICE plots for env.fact k for each algorithm
        temp.list.plots <- vector(mode = 'list', length = no.algo)
        names(temp.list.plots) <- list.algo
        
        for(l in list.algo){
            
            # Extract trained model
            trained.mod <- outputs[[l]][[taxa]][["Trained model"]]
            if(trained.mod != "NULL_MODEL"){ # need this because some models are "NULL" (svmRadial didn't work for all taxa)
                
                # Extract environmental dataframe of each sample
                env.df <- outputs[[l]][[j]][["Observation training set"]][, env.fact]
                
                # Make range of values to test for env. fact. k
                no.steps <- 10 # 30 was way too long to compute ....
                m <- min(env.df[,k])
                M <- max(env.df[,k])
                range.test <- seq(m, M, length.out = no.steps)
                
                # Make range of backward normalized values for labelling x axis
                # !! Might mathematically false 
                m2 <- (m * normalization.data$SD[k]) + normalization.data$Mean[k]
                M2 <- (M * normalization.data$SD[k]) + normalization.data$Mean[k]
                range.orig.fact <- round(seq(m2, M2, length.out = no.steps), digits = 1)
                
                # I think that the 2'600 samples are too much, we could maybe randomly select 500 of these
                set.seed(2021)
                env.df <- env.df[sample(nrow(env.df), size = 300),]
                no.samples <- nrow(env.df)
                
                # Make a dataframe for predicted values for each sample
                pred.df <- data.frame(matrix(nrow = no.samples, ncol = no.steps))
                
                for(n in 1:no.samples){
                    for (s in 1:no.steps) {
                        
                        # Make test vector for sample n with with each value in range s of env. fact k 
                        env.fact.test <- env.df[n,]
                        env.fact.test[,k] <- range.test[s]
                        
                        # Use function predict with type probability for each test env. vector
                        pred.df[n,s] <- predict(trained.mod, env.fact.test, type = 'prob')[,"present"]
                    }
                }
                
                plot.data <- melt(pred.df)
                plot.data$rind <- 1:no.samples
                
                p <- ggplot(plot.data, aes(x = variable, y = value, group=factor(rind))) 
                p <- p + geom_line(aes(color=factor(rind)), alpha = 0.3, show.legend = FALSE) # remove legend that would be the number of samples
                p <- p + labs(title = paste("Model:", l),
                              x = k,
                              y = "Predicted probability")
                p <- p + scale_x_discrete(labels = factor(range.orig.fact))
                p <- p + theme_bw(base_size = 20)
                
                temp.list.plots[[l]] <- p
                
            }
        }
        
        title <- paste("ICE of", taxa, "for", k)
        q <- grid.arrange(grobs = temp.list.plots, ncol = 2)
        list.plots[[k]] <- q
    }
}


# Plot PDP

plot.pdp <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.taxa <- length(env.fact)
  
  # make a vector of colors
  col.vect <- names(list.algo)
  names(col.vect) <- list.algo
  col.vect <- c("Null_model" = "#000000", col.vect)
  
  list.plots <- list()
  
  if ( algo == "all"){
    
    for(j in list.taxa){
      plot.title <- paste("PDP for", j)
      plot.data <- data.frame(matrix(ncol = 4, nrow = 0))
      colnames(plot.data) <- c("value","factor","model","fct")
      plot.data$factor <- as.character(plot.data$factor)
      plot.data$model <- as.character(plot.data$model)
        
      for (k in 1:no.env.fact) {
        print(paste("Producing PDP of", k, env.fact[k], "for", j,  j))
        temp.df <- data.frame(pdp::partial(outputs[[1]][[j]][["Trained model"]], pred.var = env.fact[k])[,env.fact[k]])
        colnames(temp.df) <- "value"
        temp.df$factor <- env.fact[k]
        for (l in list.algo){ 
        temp.df[,l] <- pdp::partial(outputs[[l]][[j]][["Trained model"]], pred.var = env.fact[k])[,"yhat"]
        }
        temp.df <- gather(temp.df, key = model, value = fct, -value, - factor)
        plot.data <- union(plot.data,temp.df)
      }
      
      p <- ggplot(plot.data, aes(x = value, y = fct, color = model))
      p <- p + geom_line()
      p <- p + facet_wrap( ~ factor, scales = "free_x", 
                           #labeller=label_parsed, 
                           strip.position="bottom")
      p <- p  + labs(# title = plot.title,
                     x = "",
                     y = plot.title, # paste("f( environmental factor )"),
                     color = "Model")
      p <- p + scale_colour_manual(values=col.vect)
      p <- p + theme_bw(base_size = 20)
      p <- p + theme(axis.text=element_text(size=14),
                       plot.title = element_blank())
      p <- p + guides(colour = guide_legend(override.aes = list(size=6)))
      p <- p + coord_cartesian(ylim = c(-2,2))
      # p <- p + ggtitle(plot.title)
      # print(p)
      list.plots[[j]] <- p
    }
  } else {
    for(j in 1:no.taxa){
      temp.list.plots <- vector(mode = 'list', length = no.env.fact)
      names(temp.list.plots) <- env.fact  
      
      for (k in 1:no.env.fact) {
          print(paste("Producing PDP of", k, env.fact[k], "for", j,  j))
          plot.title <- paste("PDP of", env.fact[k])
        
        p <- partial(outputs[[algo]][[j]][["Trained model"]], pred.var = env.fact[k])
        p <- autoplot(p, ylab = paste("f(", env.fact[k], ")")) +
          theme_light() +
          ggtitle(plot.title)
        p <- p + coord_cartesian(ylim = c(-2,2))
        temp.list.plots[[k]] <- p
        
      }
      title <- paste(algo, "applied to", j)
      p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
        
      list.plots[[j]] <- p
    }
    
  }
  return(list.plots)
}

# Plot multiple predictors PDP

plot.mult.pred.pdp <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.taxa <- length(env.fact)
  
  list.plots <- list()
  n = 1
  
  # For now set number of env fact, algo and taxa to minimum
  # Later can be looped but computationally heavy
  k1 = 5
  k2 = 6
  l = 2
  j = 1
  
  p <- outputs[[l]][[j]][["Trained model"]] %>%
    partial(pred.var = c(env.fact[k1], env.fact[k2]))
  
  # p1 <- plotPartial(p,# default 2 fact. PDP
  #                   main = "Basic 2 predictors PDP")
  
  rwb <- colorRampPalette(c("blue", "white", "red"))
  p2 <- plotPartial(p, contour = TRUE, # add contour lines and colors
                    col.regions = rwb,
                    main = paste(l, "applied to", j, "\nPDP with two predictors"))
  
  # p3 <- plotPartial(p, levelplot = FALSE, # 3D surface
  #                   zlab = "bruh", colorkey = TRUE,
  #                   screen = list(z = -20, x = -60),
  #                   main = "3D surface 2 predict. PDP")
  
  # title <- paste(l, "applied to", j)
  # p <- as_ggplot(grid.arrange(p1, p2, p3, ncol = 3, top = title))
    
  list.plots[[n]] <- p2
  n <- n + 1
  
  return(list.plots)
}

# Plot ICE

plot.ice <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.fact <- length(env.fact)
  
  list.plots <- list()
  n = 1
  
  if ( algo == "all"){
    for(j in 1:no.taxa){
      for (l in 1:no.algo){ 
        
        temp.list.plots <- vector(mode = 'list', length = no.env.fact)
        names(temp.list.plots) <- env.fact  
        
        for (k in 1:no.env.fact) {
          plot.title <- paste("ICE of", env.fact[k])
          
          p <- partial(outputs[[l]][[j]][["Trained model"]], pred.var = env.fact[k], ice = TRUE, alpha = 0.1)
          p <- autoplot(p, smooth = TRUE, ylab = paste("f(", env.fact[k], ")")) +
            theme_light() +
            ggtitle(plot.title)
          p <- p + coord_cartesian(ylim = c(-2,2))
          temp.list.plots[[k]] <- p
          
        }
        title <- paste(l, "applied to", j)
        p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
          
        list.plots[[n]] <- p
        n <- n + 1
      }
    }  
  } else {
    for(j in 1:no.taxa){
        
        temp.list.plots <- vector(mode = 'list', length = no.env.fact)
        names(temp.list.plots) <- env.fact  
        
        for (k in 1:no.env.fact) {
            
            print(paste("Producing ICE of", k, env.fact[k], "for", j,  j))
            plot.title <- paste("ICE of", env.fact[k])
          
          p <- partial(outputs[[algo]][[j]][["Trained model"]], pred.var = env.fact[k], ice = TRUE, alpha = 0.1)
          p <- autoplot(p, smooth = TRUE, ylab = paste("f(", env.fact[k], ")")) +
            theme_light() +
            ggtitle(plot.title) +
            coord_cartesian(ylim = c(-2,2))
          temp.list.plots[[k]] <- p
          
        }
        title <- paste(algo, "applied to", j)
        p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
          
        list.plots[[n]] <- p
        n <- n + 1
    }
  }
  return(list.plots)
}


map.ml.pred.taxa <- function(taxa, inputs, outputs, list.algo, CV){
        
        m <- ifelse(CV, "testing set", "training set")
        taxon <- sub("Occurrence.", "", taxa)
        cat("Constructing ggplot for:", taxon, "\n")
        
        df.st.dev <- data.frame("model" = list.algo)
            
        temp.df <- data.frame(outputs[[l]][[taxa]][[paste("Observation", m)]][, c("X","Y", taxa)])
        for (l in 1:no.algo) {
            temp.df[,list.algo[l]] <- outputs[[l]][[taxa]][[paste("Prediction probabilities",m)]][,"present"]
            df.st.dev[l, "st.dev"] <-  round(outputs[[l]][[taxa]][[paste("Performance",m)]], digits = 3)
            
        }
        plot.data <- gather(temp.df, key = model, value = pred, -X, -Y, -taxa)
        subtitle <- paste(paste(df.st.dev$model, df.st.dev$st.dev), collapse = " - ")
        
        # Map geometries
        g <- ggplot()
        g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
        g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
        g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
        g <- g + geom_point(data = plot.data, aes(plot.data[,"X"], plot.data[,"Y"], size = plot.data[,"pred"], 
                                                  color = plot.data[,taxa]), 
                            alpha =0.7) #alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
        g <- g + facet_wrap(~ model,
                            #labeller=label_parsed, 
                            strip.position="bottom")
        
        # Configure themes and labels
        g <- g + theme_void()
        g <- g + theme(plot.title = element_text(size = 16),#, hjust = 0.5),
                       panel.grid.major = element_line(colour="transparent"),
                       plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                       legend.title = element_text(size=14))
        
        g <- g + labs(title = paste("Geographic distribution to compare observation\nand model prediction on", m,":", taxa),
                      subtitle = paste("Stand. dev.:", subtitle),
                      x = "",
                      y = "",
                      size = "Probability of\noccurrence",
                      color = "Observation")
        
        # Configure legends and scales
        g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
                        # alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
                        color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
        g <- g + scale_size(range = c(1,3.5))
        g <- g + scale_color_manual(values=c(absent = "#c2141b", present = "#007139"), labels=c("Absence", "Presence"))
        g <- g + scale_shape_identity() # Plot the shape according to the data
        
        
    return(g)
    
}

response.ml.pred.taxa <- function(taxa, outputs, list.algo, env.fact, algo = list.algo[1], CV){
    
        m <- ifelse(CV, "testing set", "training set")
    
        taxon <- sub("Occurrence.", "", taxa)
        cat("Constructing ggplot for:", taxon, "\n")
        
        plot.data <- outputs[[algo]][[taxa]][[paste("Observation",m)]]
        plot.data$pred <- outputs[[algo]][[taxa]][[paste("Prediction probabilities",m)]][,"present"]
        plot.data <- gather(plot.data, key = factors, value = value, -SiteId, -SampId, -X, -Y, -taxa, -pred)
        
        g <- ggplot(data = plot.data, aes(x = value, y = pred, color = plot.data[,taxa]))
        
        g <- g + geom_point(alpha = 0.35)
        g <- g + theme_bw(base_size=15)
        g <- g + facet_wrap( ~ factors, scales = "free_x", 
                             #labeller=label_parsed, 
                             strip.position="bottom")

        g <- g + scale_color_manual(name = "Observation", values=c(absent = "#c2141b", present = "#007139"), labels = c("Absence", "Presence"))
        g <- g + labs(title = paste("Probability of occurrence vs explanatory variables:",algo, "applied to",paste(taxon)),
                      subtitle = paste("Predicted on", m),
                      x = "Explanatory variable",
                      y = "Predicted probability of occurrence",
                      color = "Observation")
        g <- g + theme(strip.background = element_blank(),
                       strip.placement = "outside",
                       plot.title = element_text(size=10))
        g <- g + ylim(0,1)
        g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
       
    return(g)
}
