# youpi

## ---- Produce pdf of a list of plots ----

print.pdf.plots <- function(list.plots, width = 12, dir.output, info.file.name = "", file.name){
  
  pdf(paste0(dir.output, info.file.name, file.name), paper = 'special', width = width, height = width*3/4, onefile = TRUE)
  
  for (n in 1:length(list.plots)) {
    
    print(list.plots[[n]])
    
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

# Compare models
model.comparison <- function(outputs, null.model, list.algo, list.taxa, prev.inv){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  
  col.vect <- c("Null model" = "#000000", 
                # "adaboost" = "#948B8B",
                # "earth" = "#048504",
                "svmRadial" = "#DB1111",
                "glm" = "#030AE8",
                # "CT2" = "#A84E05",
                "rf" = "#790FBF")
  
  list.plots <- list()
  n = 1
  # Plot performance on training and testing set
  if( is.null(outputs[[1]][[1]][["Testing performance"]]) == FALSE){ # check if testing set isn't empty
    c <- c("Training performance", "Testing performance")
  } else {
    c <- c("Training performance") # if empty, then only plot for training set
  }
  for (m in 1:length(c)){
    if(c[m] == "Training performance"){
      title <- paste("Comparison of models performance on training set")
    } else {
      title <- paste("Comparison of models performance on testing set")
    }
    plot.data <- data.frame(matrix(ncol = no.algo, nrow = no.taxa))
    colnames(plot.data) <- list.algo
    plot.data$Taxa <- list.taxa
    plot.data$Prevalence <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Prevalence"]
    plot.data[, "Taxonomic level"] <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"]
    plot.data[,"Null model"] <- NA
    expl.pow <- paste0("expl.pow_", list.algo)
    plot.data[,expl.pow] <- NA
    for (j in 1:no.taxa) {
      plot.data[j,"Null model"] <- null.model[[j]][["Performance"]]
      for(l in 1:no.algo){
        plot.data[j,list.algo[l]] <- outputs[[l]][[j]][[c[m]]]
        val.expl.pow <- null.model[[j]][["Performance"]] - outputs[[l]][[j]][[c[m]]] / null.model[[j]][["Performance"]]
        plot.data[j, paste0("expl.pow_",list.algo[l])] <- val.expl.pow
        # plot.data[j, "expl.pow_null.model"] <- NA
      }
    }
    plot.data$null.perf <- plot.data[,"Null model"]
    plot.data <- gather(plot.data, key = model, value = performance, all_of(c("Null model",list.algo)))
    plot.data <- gather(plot.data, key = expl.pow_model, value = expl.pow, all_of(paste0("expl.pow_", list.algo)))
    
    # Prevalence vs stand dev
    p1 <- ggplot()
    p1 <- p1  + geom_point(data = plot.data, aes(x = Prevalence, y = performance, 
                                                 colour = model, 
                                                 shape = plot.data[,"Taxonomic level"]), 
                           alpha = 0.4,
                           size = 3)
    p1 <- p1 + geom_line(data = plot.data, aes(x = Prevalence, y = null.perf), linetype = "dashed", alpha=0.4, show.legend = FALSE)
    p1 <- p1  + labs(y = "Standardized deviance",
                      x = "Prevalence (%)",
                      shape = "Taxonomic level",
                      color = "Model")
    p1 <- p1 + scale_colour_manual(values=col.vect)
    p1 <- p1 + theme_bw(base_size = 20)
    p1 <- p1 + theme(axis.text=element_text(size=14),
                     plot.title = element_blank())
    p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=6)))
    
    # # Prevalence vs stand dev with LOESS curves
    # p2 <- ggplot()
    # p2 <- p2 + geom_point(data = plot.data, aes(x = Prevalence, y = performance, 
    #                                             colour = model, 
    #                                             shape = plot.data[,"Taxonomic level"]), alpha = 0.4)
    # p2 <- p2 + stat_smooth(data = plot.data, aes(x=Prevalence, y = performance, 
    #                            colour = model), size = 2, geom='line', alpha=0.50, 
    #                        # se=TRUE, 
    #                        method = "loess")
    # p2 <- p2 + geom_line(data = plot.data, aes(x = Prevalence, y = null.perf), linetype = "dashed", alpha=0.4, show.legend = FALSE)
    # p2 <- p2  + labs(y = "Standardized deviance",
    #                  x = "Prevalence (%)",
    #                  shape = "Taxonomic level",
    #                  color = "Model")
    # p2 <- p2 + scale_colour_manual(values=col.vect) #,
    # p2 <- p2 + theme_bw(base_size = 20)
    # p2 <- p2 + theme(axis.text=element_text(size=14),
    #                  plot.title = element_blank())
    # p2 <- p2 + guides(colour = guide_legend(override.aes = list(size=6)))

    # Boxplots
    p3 <- ggplot(plot.data, aes(x=model, y =performance, fill = model), alpha = 0.4) 
    p3 <- p3 + geom_boxplot()
    p3 <- p3 + scale_fill_manual(values=col.vect)
    p3 <- p3 + coord_flip()
    p3 <- p3 + theme_bw(base_size = 20)
    p3 <- p3 + labs(x="Model",
                    y="Standardized deviance")
    
    p <- as_ggplot(grid.arrange(p3, 
                                p1,
                                # arrangeGrob(p1, p2, ncol = 2), 
                                top = title))
    list.plots[[n]] <- p
    n <- n + 1
  }
  
  # Plot performance for each taxa after resampling
  n = 3
  for (j in 1:no.taxa){
    
    temp.list.trmod <- vector(mode = 'list', length = no.algo)
    names(temp.list.trmod) <- list.algo
    for (l in 1:no.algo){
      temp.list.trmod[[l]] <- outputs[[l]][[j]][["Trained model"]]
    }
    
    # Use caret function resamples 
    models.compare <- resamples(temp.list.trmod, metric = "StandardizedDeviance")
    summary(models.compare)
    
    # Draw box plots to compare models from this resampling
    scales <- list(x=list(relation="free"), y=list(relation="free"))
    title <- paste("Comparison of model performance after resampling for",list.taxa[[j]])
    list.plots[[n]] <- bwplot(models.compare, scales=scales, main = title) # see how algorithms perform in terms of metric
    n = n + 1  
    }
  return(list.plots)
}

#  Plot variable importance
plot.varimp <- function(outputs, list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  
  list.plots <- list()
  n = 1
  
  for(j in 1:no.taxa){
    
    temp.list.plots <- vector(mode = 'list', length = no.algo)
    names(temp.list.plots) <- list.algo  
    
    for (l in 1:no.algo){
      
      # Show how the various iterations of hyperparameter search performed
      # print(plot(outputs[[l]][[j]][["Trained model"]], main = title))
      
      # Plot variable importance
      perf <- round(outputs[[l]][[j]][["Training performance"]], digits = 3)
      plot.title <- paste("Algo:", list.algo[l], 
                          "\nStandardized deviance:", 
                          ifelse(length(perf) > 1, perf, perf[1] ))
      
      temp.list.plots[[l]] <- ggplot(outputs[[l]][[j]][["Variable importance"]]) +
        ggtitle(plot.title)
      print(c(list.taxa[j], list.algo[l]))
      
      # # Plot confusion matrix
      # plot(output[[j]][["Confusion matrix"]][["table"]], main = list.taxa[j])
      # 
      
      }
    title <- paste(list.taxa[j])
    p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
      
    list.plots[[n]] <- p
    n <- n + 1
    }
  return(list.plots)
}

plot.table.varimp <- function(outputs, list.algo, list.taxa, env.fact){
  
  temp.list <- list()
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.fact <- length(env.fact)
  
  temp.df <- data.frame(matrix(ncol = no.algo*no.taxa, nrow = no.env.fact))
  colnames(temp.df) <- c(outer(list.algo, list.taxa, FUN = paste))
  rownames(temp.df) <- env.fact
  
  for (j in 1:no.taxa) {
    for (l in 1:no.algo) {
      for (k in 1:no.env.fact) {
        temp.df[env.fact[k],paste(list.algo[l], list.taxa[j])] <- outputs[[l]][[j]][["Variable importance"]][["importance"]][env.fact[k],1]
      }
    }
  }
  temp.df$mean.imp <- rowMeans(temp.df)
  temp.df <- as.matrix(temp.df)
  temp.list[[1]] <- plot(temp.df, 
                           #key = NULL,
                           digits = 2, text.cell=list(cex=0.5),
                           axis.col=list(side=3, las=2), axis.row = list(side=2, las=1),
                           col = viridis,
                           xlab = "",
                           ylab = "Environmental factor",
                           main = "Variable importance for ML algorithm applied to taxa")
 return(temp.list)
}

# Plot PDP

plot.pdp <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.taxa <- length(env.fact)
  
  col.vect <- c(# "Null model" = "#000000", 
                # "adaboost" = "#948B8B",
                # "earth" = "#048504",
                "svmRadial" = "#DB1111",
                "glm" = "#030AE8",
                # "CT2" = "#A84E05",
                "rf" = "#790FBF")
  
  list.plots <- list()
  
  if ( algo == "all"){
    
    for(j in 1:no.taxa){
      plot.title <- paste("PDP for", list.taxa[j])
      plot.data <- data.frame(matrix(ncol = 4, nrow = 0))
      colnames(plot.data) <- c("value","factor","model","fct")
      plot.data$factor <- as.character(plot.data$factor)
      plot.data$model <- as.character(plot.data$model)
        
      for (k in 1:no.env.fact) {
        print(paste("Producing PDP of", k, env.fact[k], "for", j,  list.taxa[j]))
        temp.df <- data.frame(partial(outputs[[1]][[j]][["Trained model"]], pred.var = env.fact[k])[,env.fact[k]])
        colnames(temp.df) <- "value"
        temp.df$factor <- env.fact[k]
        for (l in 1:no.algo){ 
        temp.df[,list.algo[l]] <- partial(outputs[[l]][[j]][["Trained model"]], pred.var = env.fact[k])[,"yhat"]
        }
        temp.df <- gather(temp.df, key = model, value = fct, -value, - factor)
        plot.data <- union(plot.data,temp.df)
      }
      
      p <- ggplot(plot.data, aes(x = value, y = fct, color = model))
      p <- p + geom_line()
      p <- p + geom_smooth(method = "loess")
      p <- p + facet_wrap( ~ factor, scales = "free_x", 
                           #labeller=label_parsed, 
                           strip.position="bottom")
      p <- p  + labs(x = "",
                     y = paste("f( environmental factor )"),
                     color = "Model")
      p <- p + scale_colour_manual(values=col.vect)
      p <- p + theme_bw(base_size = 20)
      p <- p + theme(axis.text=element_text(size=14),
                       plot.title = element_blank())
      p <- p + guides(colour = guide_legend(override.aes = list(size=6)))
      p <- p + ggtitle(plot.title)
      # print(p)
      list.plots[[j]] <- p
    }
  } else {
    for(j in 1:no.taxa){
      temp.list.plots <- vector(mode = 'list', length = no.env.fact)
      names(temp.list.plots) <- env.fact  
      
      for (k in 1:no.env.fact) {
        plot.title <- paste("PDP of", env.fact[k])
        
        p <- partial(outputs[[algo]][[j]][["Trained model"]], pred.var = env.fact[k])
        p <- autoplot(p, smooth = TRUE, ylab = paste("f(", env.fact[k], ")")) +
          theme_light() +
          ggtitle(plot.title)
        temp.list.plots[[k]] <- p
        
      }
      title <- paste(algo, "applied to", list.taxa[j])
      p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
        
      list.plots[[n]] <- p
      n <- n + 1
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
  k1 = 1
  k2 = 2
  l = 1
  j = 1
  
  p <- outputs[[l]][[j]][["Trained model"]] %>%
    partial(pred.var = c(env.fact[k1], env.fact[k2]))
  
  p1 <- plotPartial(p,# default 2 fact. PDP
                    main = "Basic 2 predictors PDP")
  
  rwb <- colorRampPalette(c("red", "white", "blue"))
  p2 <- plotPartial(p, contour = TRUE, # add contour lines and colors
                    col.regions = rwb,
                    main = "Tuned 2 predictors PDP")
  
  p3 <- plotPartial(p, levelplot = FALSE, # 3D surface
                    zlab = "bruh", colorkey = TRUE,
                    screen = list(z = -20, x = -60),
                    main = "3D surface 2 predict. PDP")
  
  title <- paste(list.algo[l], "applied to", list.taxa[j])
  p <- as_ggplot(grid.arrange(p1, p2, p3, ncol = 3, top = title))
    
  list.plots[[n]] <- p
  n <- n + 1
  
  return(list.plots)
}

# Plot ICE

plot.ice <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
  no.algo <- length(list.algo)
  no.taxa <- length(list.taxa)
  no.env.taxa <- length(env.fact)
  
  list.plots <- list()
  n = 1
  
  if ( algo == "all"){
    for(j in 1:no.taxa){
      for (l in 1:no.algo){ 
        
        temp.list.plots <- vector(mode = 'list', length = no.env.fact)
        names(temp.list.plots) <- env.fact  
        
        for (k in 1:no.env.fact) {
          plot.title <- paste("ICE of", env.fact[k])
          
          p <- partial(outputs[[l]][[j]][["Trained model"]], pred.var = env.fact[k], ice = TRUE, alpha = 0.5)
          p <- autoplot(p, smooth = TRUE, ylab = paste("f(", env.fact[k], ")")) +
            theme_light() +
            ggtitle(plot.title)
          temp.list.plots[[k]] <- p
          
        }
        title <- paste(list.algo[l], "applied to", list.taxa[j])
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
          plot.title <- paste("ICE of", env.fact[k])
          
          p <- partial(outputs[[algo]][[j]][["Trained model"]], pred.var = env.fact[k], ice = TRUE, alpha = 0.5)
          p <- autoplot(p, smooth = TRUE, ylab = paste("f(", env.fact[k], ")")) +
            theme_light() +
            ggtitle(plot.title)
          temp.list.plots[[k]] <- p
          
        }
        title <- paste(algo, "applied to", list.taxa[j])
        p <- as_ggplot(grid.arrange(grobs = temp.list.plots, top = title))
          
        list.plots[[n]] <- p
        n <- n + 1
    }
  }
  return(list.plots)
}

map.ml.pred.taxa <- function(inputs, outputs, list.taxa, list.algo, algo = list.algo[1]){
    
    list.plots <- vector(mode = 'list', length = length(list.taxa))
    names(list.plots) <- list.taxa
    
    for(j in 1:length(list.taxa)){
        
        taxon <- sub("Occurrence.", "", list.taxa[j])
        cat("Constructing ggplot for:", taxon, "\n")
        
        plot.data <- outputs[[algo]][[list.taxa[j]]][["Observation testing set"]]
        plot.data$pred <- outputs[[algo]][[list.taxa[j]]][["Prediction on testing set (probabilities)"]][,"present"]
        
        # Map geometries
        g <- ggplot()
        g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
        g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
        g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
        g <- g + geom_point(data = plot.data, aes(X, Y, size = pred, 
                                                  color = plot.data[,list.taxa[j]]), 
                            alpha =0.7) #alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))

        # Configure themes and labels
        g <- g + theme_void()
        g <- g + theme(plot.title = element_text(size = 16),#, hjust = 0.5),
                       panel.grid.major = element_line(colour="transparent"),
                       plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                       legend.title = element_text(size=14))
        
        g <- g + labs(title = paste("Geographic distribution to compare observation\nand model prediction:", 
                                    algo, "applied to", taxon),
                      # subtitle = paste("Random Forest:", paste(env.fact, collapse = " ", sep = " ")), #, "- page", j),
                      x = "",
                      y = "",
                      size = "Probability of\noccurrence",
                      color = "Observation")
        
        # Configure legends and scales
        g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
                        # alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
                        color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
        g <- g + scale_color_manual(values=c(absent = "#c2141b", present = "#007139"), labels=c("Absence", "Presence"))
        g <- g + scale_shape_identity() # Plot the shape according to the data
        
        list.plots[[j]] <- g
    }
    return(list.plots)
}

response.ml.pred.taxa <- function(outputs, list.algo, list.taxa, env.fact, algo = list.algo[1]){
    
    list.plots <- vector(mode = 'list', length = length(list.taxa))
    names(list.plots) <- list.taxa
    
    for(j in 1:length(list.taxa)){
        
        taxon <- sub("Occurrence.", "", list.taxa[j])
        cat("Constructing ggplot for:", taxon, "\n")
        
        plot.data <- outputs[[algo]][[list.taxa[j]]][["Observation testing set"]]
        plot.data$pred <- outputs[[algo]][[list.taxa[j]]][["Prediction on testing set (probabilities)"]][,"present"]
        plot.data <- gather(plot.data, key = factors, value = value, -SiteId, -SampId, -X, -Y, -list.taxa[j], -pred)
        
        g <- ggplot(data = plot.data, aes(x = value, y = pred, color = plot.data[,list.taxa[j]]))
        
        g <- g + geom_point(alpha = 0.35)
        g <- g + theme_bw(base_size=15)
        g <- g + facet_wrap( ~ factors, scales = "free_x", 
                             #labeller=label_parsed, 
                             strip.position="bottom")

        g <- g + scale_color_manual(name = "Observation", values=c(absent = "#c2141b", present = "#007139"), labels = c("Absence", "Presence"))
        g <- g + labs(title = paste("Probability of occurrence vs explanatory variables:",algo, "applied to",paste(taxon)),
                      x = "Explanatory variable",
                      y = "Probability of occurrence",
                      color = "Observation")
        g <- g + theme(strip.background = element_blank(),
                       strip.placement = "outside",
                       plot.title = element_text(size=10))
        g <- g + ylim(0,1)
        g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
        list.plots[[j]] <- g
    }
    return(list.plots)
}
