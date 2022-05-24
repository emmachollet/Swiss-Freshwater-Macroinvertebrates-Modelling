# youpi

## ---- Produce pdf of a list of plots ----

print.pdf.plots <- function(list.plots, width = 12, height = width*3/4, dir.output, info.file.name = "", file.name){
  
  pdf(paste0(dir.output, info.file.name, file.name), paper = 'special', width = width, height = height, onefile = TRUE)
  
  for (n in 1:length(list.plots)) {
    
    plot(list.plots[[n]]) # sometimes this works better
    # print(list.plots[[n]])
    
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

plot.df.perf <- function(df.perf, list.models, list.taxa, CV, title = c()){
    
    size.val <- ifelse(length(list.taxa) < 25, 5, 1)
    if(length(title) != 0){
      title <- title
    } else{
      title <- ifelse(CV, "Table of models predictive performance \n during Cross Validation",
                    "Table of models quality of fit \n during Calibration")
    }
    
    if(grepl("likelihood", title)){
        limits <- c(0,1)
        midpoint <- 0.7
        low <- "#c2141b"
        mid <- "#f9c74f"
        high <- "#007139"
        legend <- "Likelihood \nratio"
    } else {
        limits <- c(0,2)
        midpoint <- 0.7
        low <- "#007139"
        mid <- "grey70"
        high <- "#c2141b"
        legend <- "Standardized \ndeviance"
    }
    
    if(CV){ cind <- 1:which(colnames(df.perf) == "Taxa")
    } else { cind <- 1:which(colnames(df.perf) == "Taxa")
      # cind <- 1:which(colnames(df.perf) %in% list.models | colnames(df.perf) == "Taxa") # ECR: was there for FIT
      }
    temp.df <- df.perf[,cind]
    # temp.df$models <- row.names(df.perf)
    melted.df <- melt(temp.df, id = "Taxa")
    p <- ggplot(data = melted.df, aes(x = Taxa, y = variable, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(midpoint = midpoint, low = low, mid = mid, high = high, 
                             limits = limits) +
        scale_x_discrete(limits = list.taxa) +
        labs(title = title, 
             x = "", y = "", fill = legend) +
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
model.comparison <- function(df.merged.perf, list.models, CV, extrapol, select.taxa){
    
    list.models.temp <- list.models
    list.models <- c("#000000" = "Null_model", list.models)
    
    # Make a vector of colors
    col.vect <- names(list.models)
    names(col.vect) <- list.models
    
    no.taxa <- nrow(df.merged.perf)
    
    title <- c("Models performance during calibration", "Models performance during prediction")
    
    if(!CV & !extrapol){
        title <- title[1]
    }
    
    subtitle <- paste("For", no.taxa, "taxa")
    
    # ECR: To be completed if CV = F #### 
    
    if(CV | extrapol){  
      
      plot.data <- df.merged.perf

      col.fit <- c(paste0(list.models.temp , ".fit"), "Null_model")
      col.pred <- c(paste0(list.models.temp , ".pred"), "Null_model")

      plot.data1 <- gather(plot.data, key = model, value = performance.fit, all_of(col.fit))
      plot.data2 <- gather(plot.data, key = model, value = performance.pred, all_of(col.pred))
      plot.data1$model <- sub(".fit", "", plot.data1$model)
      plot.data1 <- plot.data1[,c("Taxa", "Prevalence", "Taxonomic level", "model", "performance.fit")]
      plot.data2$model <- sub(".pred", "", plot.data2$model)
      plot.data2 <- plot.data2[,c("Taxa", "Prevalence", "Taxonomic level", "model", "performance.pred")]
      
      plot.data <- left_join(plot.data1, plot.data2, by = c("Taxa", "Prevalence", "Taxonomic level", "model"))

      plot.data3 <- gather(plot.data, key = dataset, value = performance, -c("Taxa", "Prevalence", "Taxonomic level", "model") )
      plot.data3$dataset <- sub("performance.fit", "Calibration", plot.data3$dataset)
      plot.data3$dataset <- sub("performance.pred", "Prediction", plot.data3$dataset)
      
      }
    
    list.plots.temp <- vector(mode = "list", length = 2)
    
    for(n in 1:length(title)){
      list.plots.temp[[n]] <- vector(mode = "list", length = 4)
      perf <- paste0("performance", ifelse(n == 1, ".fit", ".pred" ))
      
      # Prevalence vs stand dev
      p1 <- ggplot()
      p1 <- p1  + geom_point(data = plot.data, aes_string(x = "Prevalence", y = perf, 
                                                   colour = "model"# , 
                                                   # shape = "Taxonomic level"
                                                   ), 
                             # alpha = 0.4,
                             size = 3)
      p1 <- p1 + xlim(0.045, 0.955)
      #p1 <- p1 + geom_line(data = plot.data, aes(x = Prevalence, y = null.perf), linetype = "dashed", alpha=0.4, show.legend = FALSE) # to plot null model as dash line between data points
      p1 <- p1 + stat_function(fun=function(x) -2*(x*log(x) + (1-x)*log(1-x))) # to plot null model as function line
      p1 <- p1  + labs(y = "Standardized deviance",
                       x = "Prevalence (%)",
                       shape = "Taxonomic level",
                       color = "Model",
                       title = title[n],
                       subtitle = subtitle)
      p1 <- p1 + scale_colour_manual(values=col.vect)
      p1 <- p1 + theme_bw(base_size = 20)
      p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=6)))
      
      list.plots.temp[[n]][[1]] <- p1
      
      # Prevalence vs stand dev
      temp.plot.data <- plot.data[which(plot.data$Taxa %in% select.taxa),]
      p1.2 <- ggplot()
      p1.2 <- p1.2  + geom_point(data = temp.plot.data, aes_string(x = "Prevalence", y = perf, 
                                                          colour = "model"# , 
                                                          # shape = "Taxonomic level"
      ), 
      size = 3)
      p1.2 <- p1.2 + stat_function(fun=function(x) -2*(x*log(x) + (1-x)*log(1-x))) # to plot null model as function line
      p1.2 <- p1.2  + labs(y = "Standardized deviance",
                       x = "Prevalence (%)",
                       shape = "Taxonomic level",
                       color = "Model",
                       title = title[n],
                       subtitle = paste("For", length(select.taxa), "taxa"))
      p1.2 <- p1.2 + scale_colour_manual(values=col.vect)
      p1.2 <- p1.2 + theme_bw(base_size = 20)
      p1.2 <- p1.2 + guides(colour = guide_legend(override.aes = list(size=6)))
      
      list.plots.temp[[n]][[2]] <- p1.2

      # Boxplots
      p2 <- ggplot(plot.data, aes_string(x="model", y = perf, fill = "model"), alpha = 0.4) 
      p2 <- p2 + geom_boxplot()
      p2 <- p2 + scale_fill_manual(values=col.vect)
      p2 <- p2 + labs(title = title)
      # p2 <- p2 + ylim(0,2) # ECR: only because perf problems
      p2 <- p2 + scale_x_discrete(limits = rev(list.models))
      p2 <- p2 + coord_flip()
      p2 <- p2 + theme_bw(base_size = 20)
      p2 <- p2 + theme(legend.position = "none")
      p2 <- p2 + labs(x="Models",
                      y="Standardized deviance",
                      # fill = "Models",
                      title = title[n],
                      subtitle = subtitle)
      
      list.plots.temp[[n]][[3]] <- p2

      # Boxplots
      p3 <- ggplot(plot.data, aes(y = reorder(model,-!!ensym(perf)), x = !!ensym(perf), fill = model), alpha = 0.4)
      p3 <- p3 + geom_boxplot()
      p3 <- p3 + scale_fill_manual(values=col.vect)
      p3 <- p3 + labs(title = title)
      # p3 <- p3 + ylim(0.2,1.5) # ECR: only because perf problems
      # p3 <- p3 + scale_x_discrete(limits = rev(list.models))
      #p3 <- p3 + coord_flip()
      p3 <- p3 + theme_bw(base_size = 20)
      p3 <- p3 + theme(legend.position = "none")
      p3 <- p3 + labs(x="Models",
                      y="Standardized deviance",
                      # fill = "Models",
                      title = title[n],
                      subtitle = subtitle)
                        # paste0(subtitle,
                        #                "\nOrdered by decreasing mean"))
      
      list.plots.temp[[n]][[4]] <- p3
    }
    
    if(CV | extrapol){
      list.plots <- vector(mode = "list", length = 4)
      
      for (n in 1:length(list.plots)) {
        list.plots[[n]] <- grid.arrange(grobs = list(list.plots.temp[[1]][[n]], list.plots.temp[[2]][[n]]), ncol = 2)
      }
      
      # All boxplots on one plot
      # To order them by prediction performance
      lev <- levels(with(plot.data2, reorder(model,performance.pred)))
      p4 <- ggplot(plot.data3, aes(x = model, y = performance, fill = dataset))
      p4 <- p4 + geom_boxplot()
      p4 <- p4 + scale_x_discrete(limits = lev)
      p4 <- p4 + ylim(0,ifelse(extrapol, 2.5, 1.5)) # ECR: only because perf problems
      p4 <- p4 + scale_fill_manual(values=c(Calibration = "#f1a340", Prediction = "#998ec3"))
      p4 <- p4 + theme_bw(base_size = 20)
      p4 <- p4 + labs(x="Models",
                      y="Standardized deviance",
                      fill = "",
                      title = "Models performance comparison",
                      subtitle = paste0(subtitle, "\nOrdered by increasing mean"))
      
      
      list.plots <- append(list.plots, list(p4))
      
      select.model <- c("RF", "iGLM")
      temp.plot.data <- plot.data3[which(plot.data3$Taxa %in% select.taxa),]
      temp.plot.data <- temp.plot.data[which(temp.plot.data$model %in% select.model),]
      
      p5 <- ggplot()
      p5 <- p5  + geom_point(data = temp.plot.data, aes(x = Prevalence, y = performance, 
                                                                   colour = model, 
                                                                   shape = dataset), size = 5)
      p5 <- p5 + stat_function(fun=function(x) -2*(x*log(x) + (1-x)*log(1-x))) # to plot null model as function line
      p5 <- p5  + labs(y = "Standardized deviance",
                           x = "Prevalence (%)",
                           shape = "",
                           color = "Model",
                           subtitle = paste("For", length(select.taxa), "taxa"))
      p5 <- p5 + scale_colour_manual(values=col.vect[select.model])
      p5 <- p5 + theme_bw(base_size = 25)
      p5 <- p5 + guides(colour = guide_legend(override.aes = list(size=6)))
      
      list.plots <- append(list.plots, list(p5))
      
    } else {
      list.plots <- list.plots.temp
    }
    
    return(list.plots)
    
}

plot.perf.fitvspred <- function(df.fit.perf, df.pred.perf, list.models, select.taxa = list.taxa){
  
  list.models.temp <- c("grey30" = "Null_model")
  list.models <- c(list.models.temp, list.models)
  
  # Make a vector of colors
  col.vect <- names(list.models)
  names(col.vect) <- list.models
  
  # Make dataframe to plot
  plot.data.fit <- df.fit.perf[,-which(grepl("expl.pow",colnames(df.fit.perf)))] %>%
    gather(key = model, value = performance.fit, -c("Taxa", "Prevalence", "Taxonomic level"))
  plot.data.pred <- df.pred.perf[,-which(grepl("expl.pow",colnames(df.pred.perf)))] %>%
    gather(key = model, value = performance.pred, -c("Taxa", "Prevalence", "Taxonomic level"))
  plot.data <- left_join(plot.data.fit, plot.data.pred, by = c("Taxa", "Prevalence", "Taxonomic level", "model"))
  plot.data$Balance <- abs(plot.data$Prevalence - 0.5) * 200
  
  title <- "Comparison of performance"
  p1 <- ggplot(plot.data) +
    geom_point(aes(x = performance.fit, y = performance.pred, 
                   colour = model, size = Balance, shape = plot.data[,"Taxonomic level"]), 
               alpha = 0.8
    ) +
    # geom_encircle(aes(x = performance.fit, y = performance.pred), alpha = 0.2, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="dimgrey", linetype="dashed") +
    xlim(0,1.5) + ylim(0,1.5) + 
    labs(y = "Performance during prediction",
         x = "Performance during training",
         shape = "Taxonomic level",
         color = "Model",
         size = "Taxon balance (%)",
         title = title) + 
    scale_colour_manual(values=col.vect) + 
    theme_bw()
  
  # Close up view
  title <- "Comparison of performance"
  p2 <- ggplot(plot.data) +
    geom_point(aes(x = performance.fit, y = performance.pred, 
                   colour = model, size = Balance, shape = plot.data[,"Taxonomic level"]), 
               alpha = 0.8
    ) +
    # geom_encircle(aes(x = performance.fit, y = performance.pred), alpha = 0.2, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="dimgrey", linetype="dashed") +
    xlim(0.25,1.25) + ylim(0.25,1.25) + 
    labs(y = "Performance during prediction",
         x = "Performance during training",
         shape = "Taxonomic level",
         color = "Model",
         size = "Taxon balance (%)",
         title = title,
         subtitle = "Close up view") + 
    scale_colour_manual(values=col.vect) + 
    theme_bw()
  
  # Taxa selection
  plot.data2 <- plot.data[which(plot.data$Taxa %in% select.taxa),]
  plot.data2$Taxa <- as.factor(plot.data2$Taxa)
  no.select.taxa <- length(select.taxa)
  title <- "Comparison of performance"
  subtitle <- paste("Selection of", no.select.taxa, "taxa with biggest expl. pow. diff. across models")
  p3 <- ggplot(plot.data2) +
    geom_point(aes(x = performance.fit, y = performance.pred, 
                   colour = model, size = Balance, shape = Taxa), 
               alpha = 0.8
    ) +
    # geom_encircle(aes(x = performance.fit, y = performance.pred), alpha = 0.2, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, color="dimgrey", linetype="dashed") +
    # xlim(0.25,1.25) + ylim(0.25,1.25) + 
    labs(y = "Performance during prediction",
         x = "Performance during training",
         shape = "Taxon",
         color = "Model",
         size = "Taxon balance (%)",
         title = title,
         subtitle = subtitle) + 
    scale_colour_manual(values=col.vect) + 
    theme_bw()

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

# Plot ICE 

plot.ice.per.taxa <- function(taxa, outputs, list.models, env.fact, select.env.fact, normalization.data, extrapol, no.samples, no.steps, subselect){
    
    cat("\nProducing ICE plot for taxa", taxa)
    
    no.models <- length(list.models)
    no.subselect <- length(subselect)
    env.fact.orig <- env.fact
    
    # Initialize final list with all arranged plots
    list.plots <- list()
    
    for(k in select.env.fact){
      # k <- env.fact[1]
      cat("\nFor env. fact.", k)
    
      # # Make temporary list of ICE plots for env.fact k for each algorithm
      # temp.list.plots <- vector(mode = 'list', length = no.models)
      # names(temp.list.plots) <- list.models
      # 
      plot.data <- data.frame()
      plot.data.means <- data.frame()
      plot.data.rug <- data.frame()
      
      for(l in list.models){
          # l <- list.models[1]
          
        cat("\nfor model", l)
        
        # Extract trained model
        #trained.mod <- ann.outputs.cv$Split1$ANN$Occurrence.Gammaridae$`Trained model`
        trained.mod <- outputs[[l]][[taxa]][["Trained model"]]

        if(grepl("GLM", l)){
          env.fact <- c(env.fact, "temperature2", "velocity2")
        } else {
          env.fact <- env.fact.orig
        }
        
        # Extract environmental dataframe of each sample
        if(extrapol){
          env.df <- outputs[[l]][[taxa]][["Observation training set"]][, env.fact]
          all.env.df <- bind_rows(outputs[[l]][[taxa]][["Observation training set"]][, env.fact],
                              outputs[[l]][[taxa]][["Observation testing set"]][, env.fact])
        } else {
          env.df <- outputs[[l]][[taxa]][["Observation training set"]][, env.fact]
        }
        
        # Make range of values to test for env. fact. k
        # no.steps <- 200
        if(extrapol){
          m <- min(all.env.df[,k])
          M <- max(all.env.df[,k])
        } else {
          m <- min(env.df[,k])
          M <- max(env.df[,k])
        }
        range.test <- seq(m, M, length.out = no.steps)
        
        # Make range of backward normalized values for labelling x axis
        # !! Might mathematically false 
        m2 <- (m * normalization.data$SD[k]) + normalization.data$Mean[k]
        M2 <- (M * normalization.data$SD[k]) + normalization.data$Mean[k]
        range.orig.fact <- round(seq(m2, M2, length.out = no.steps), digits = 1)
        
        # Subselect a number of samples/observations to plot
        # no.samples <- 3
        set.seed(2021)
        env.df <- env.df[sample(nrow(env.df), size = no.samples),]
        
        # Make a dataframe for predicted values for each sample
        pred.df <- data.frame(matrix(nrow = no.samples, ncol = no.steps))
        
        for(n in 1:no.samples){
          # n <- 1
          env.fact.test <- env.df[n,]
          env.fact.test <- env.fact.test[rep(seq_len(1), each = no.steps), ]
          env.fact.test[,k] <- range.test

          if(l == "ANN"){
            env.fact.test <- as.matrix(env.fact.test)
            pred.df[n,] <- predict(trained.mod, env.fact.test)[ , which(names(outputs[[l]]) == taxa)]
          } else if (l == "hGLM" | l == "chGLM"){
              pred.df[n,] <- t(pred.stat.models(model = res.extracted, taxa = taxa , env.fact.test = env.fact.test, list.taxa = list.taxa))
          } else {
          pred.df[n,] <- predict(trained.mod, env.fact.test, type = 'prob')[,"present"]
          }
        } 
        
        # Make dataframe for ICE
        colnames(pred.df) <- range.orig.fact

        # subselect predictions for plot resolution
        n <- 1
        temp.range.fact <- range.orig.fact[seq(1,200, by = subselect[n])]
        temp.pred.df <- pred.df[, which(colnames(pred.df) %in% temp.range.fact)]
        temp.pred.df$observation <- rownames(pred.df)
        
        temp.plot.data <- gather(temp.pred.df, key = variable, value = value, -observation)
        temp.plot.data$observation <- as.factor(temp.plot.data$observation)
        temp.plot.data$Model <- as.factor(l)
        temp.plot.data[, c("variable")] <- as.numeric(temp.plot.data[, c("variable")])
        
        plot.data <- bind_rows(plot.data, temp.plot.data)
        
        # Make dataframe for PDP (mean of ICE)
        means <- data.frame(colMeans(temp.pred.df[,1:(length(temp.pred.df)-1)]))
        colnames(means) <- "mean"
        means$variable <- temp.range.fact
        means$Model <- as.factor(l)
        plot.data.means <- bind_rows(plot.data.means, means)
        
        # Make dataframe for rug (datapoints of observations)
        observations <- data.frame(env.df[,k])
        colnames(observations) <- "variable"
        observations$variable <- (observations$variable * normalization.data$SD[k]) + normalization.data$Mean[k]
        observations$Model <- as.factor(l)
        plot.data.rug <- bind_rows(plot.data.rug, observations)
      }
        
            
          # temp.list.plots.subselect <- vector(mode = 'list', length = no.subselect)
          # 
          # for (n in 1:no.subselect) {
          #     
          #   cat("\nProducing plot for resolution of", subselect[n], "steps")
          #   temp.range.fact <- range.orig.fact[seq(1,200, by = subselect[n])]
          #   temp.pred.df <- pred.df[, which(colnames(pred.df) %in% temp.range.fact)]
          #   plot.data <- temp.pred.df
          #   plot.data <- melt(plot.data)
          #   plot.data$variable <- as.numeric(as.character(plot.data$variable))
          #   plot.data$rind <- 1:no.samples
          #   
          #   # Make dataframe for rug (datapoints of observations)
          #   observations <- data.frame(env.df[,k])
          #   colnames(observations) <- "variable"
          #   observations$variable <- (observations$variable * normalization.data$SD[k]) + normalization.data$Mean[k]
          #   
          #   # Make dataframe for PDP (mean of ICE)
          #   means <- data.frame(colMeans(temp.pred.df))
          #   colnames(means) <- "mean"
          #   means$variable <- temp.range.fact
          #   
            p <- ggplot(plot.data, aes(x = variable, y = value, group=factor(observation))) 
            p <- p + geom_line(aes(color=factor(observation)), alpha = 0.3, show.legend = FALSE) # remove legend that would be the number of samples
            p <- p + geom_line(data = plot.data.means,
                               aes(x = variable, y = mean), 
                               color = "grey30", size = 1.5, # alpha = 0.7, 
                               inherit.aes = F)
            p <- p + geom_rug(data = plot.data.rug,
                              aes(x = variable), 
                              color = "grey30", alpha = 0.7, inherit.aes = F)
            p <- p + facet_wrap( ~ Model,# scales = "free_x", 
                                 #labeller=label_parsed, 
                                 strip.position="top")
            p <- p + labs(title = "Individual Conditional Expectation and Partial Differential Plot", # paste(l),
                          # subtitle = paste("Resolution:", no.steps/no.subselect, "steps"),
                          x = k,
                          y = "Predicted probability")
            # p <- p + scale_x_discrete(labels = factor(range.orig.fact))
            p <- p + theme_bw(base_size = 10)
            
            # temp.list.plots.subselect[[n]] <- p
          # }
          
          # temp.list.plots[[l]] <- temp.list.plots.subselect
      # }
      
      # for (n in 1:no.subselect) {
      #   
      #   temp.temp.list.plots <- vector(mode = 'list', length = no.models)
      #   names(temp.temp.list.plots) <- list.models
      #   
      #   for (l in list.models) {
      #     temp.temp.list.plots[[l]] <- temp.list.plots[[l]][[n]]
      #   }
      # 
      #   for (l in list.models) {
      #     title <- paste("ICE of", sub("Occurrence.", "", taxa), "for", k, "with resolution of", no.steps/subselect[n], "steps")
      #     q <- grid.arrange(grobs = temp.temp.list.plots, ncol = 2, top = title)
      #     list.plots[[paste(k,subselect[n], sep = "_")]] <- q
      #   }
      # }
          list.plots[[paste(k,subselect[n], sep = "_")]] <- p
      }
    return(list.plots)
}

# Response shape plot

plot.rs.taxa <- function(taxa, outputs, list.models, env.fact, CV, extrapol){
  
  taxon <- sub("Occurrence.", "", taxa)
  cat("Constructing ggplot for:", taxon, "\n")
  
  temp.list.models <- c("Observations", list.models)
  
  # Make list of plots per model to be returned at the end
  list.plots <- vector(mode = 'list', length = length(temp.list.models))
  names(list.plots) <- temp.list.models
  
  for (l in temp.list.models) {
    # l <- temp.list.models[2]
    cat("For model", l, "\n")
    m <- "training set"
    
    if(l == "Observations"){
      plot.data1 <- outputs[["iGLM"]][[taxa]][[paste("Observation", m)]]
      plot.data1$pred <- ifelse(plot.data1[,taxa] == "present", 1, 0)
    } else { 
      plot.data1 <- outputs[[l]][[taxa]][[paste("Observation",m)]]
      plot.data1$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
    }
    plot.data1[,env.fact] <- as.data.frame(sweep(sweep(plot.data1[,env.fact], 2, normalization.data$SD[env.fact], FUN="*"), 2, normalization.data$Mean[env.fact], FUN = "+"))
    plot.data1 <- gather(plot.data1, key = factors, value = value, -SiteId, -SampId, -X, -Y, -taxa, -pred)
    
    if(CV | extrapol){
      plot.data1$set <- "Training"
      
      m <- "testing set"
      if(l == "Observations"){
        plot.data2 <- outputs[["iGLM"]][[taxa]][[paste("Observation", m)]]
        plot.data2$pred <- ifelse(plot.data2[,taxa] == "present", 1, 0)
      } else { 
        plot.data2 <- outputs[[l]][[taxa]][[paste("Observation", m)]]
        plot.data2$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
      }
      plot.data2[,env.fact] <- as.data.frame(sweep(sweep(plot.data2[,env.fact], 2, normalization.data$SD[env.fact], FUN="*"), 2, normalization.data$Mean[env.fact], FUN = "+"))
      plot.data2 <- gather(plot.data2, key = factors, value = value, -SiteId, -SampId, -X, -Y, -taxa, -pred)
      plot.data2$set <- "Testing"
      
      plot.data <- bind_rows(plot.data1, plot.data2)
    }
    
    g <- ggplot(data = plot.data, aes_string(x = "value", y = "pred", color = taxa))
    if(CV | extrapol){
      g <- g + geom_point(aes_string(shape = "set"), alpha = 0.35)
      g <- g + labs(shape = "")
      g <- g + guides(shape = guide_legend(override.aes = list(size=3)))
    } else {
      g <- g + geom_point(alpha = 0.35)
    }
    g <- g + theme_bw(base_size=15)
    g <- g + facet_wrap( ~ factors, scales = "free_x", 
                         #labeller=label_parsed, 
                         strip.position="bottom")
    
    g <- g + scale_color_manual(name = "Observation", values=c(absent = "#c2141b", present = "#007139"), labels = c("Absence", "Presence"))
    g <- g + labs(title = paste("Probability of occurrence vs explanatory variables"),
                  subtitle = paste(l, "-",paste(taxon)),
                  x = "Explanatory variable",
                  y = "Predicted probability of occurrence",
                  color = "Observation")
    g <- g + theme(strip.background = element_blank(),
                   strip.placement = "outside"# ,
                   # plot.title = element_text(size=10)
                   )
    g <- g + ylim(0,1)
    g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  
    list.plots[[l]] <- g
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
  k1 = env.fact[1]
  k2 = env.fact[2]
  l = list.algo[1]
  j = list.taxa[1]
  
  p <- outputs[[l]][[j]][["Trained model"]] %>%
    pdp::partial(pred.var = c(k1, k2))
  
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


map.env.fact.2 <- function(inputs, env.fact, data.env, dir.output, file.prefix){
  
  for (k in 1:length(env.fact)){
    k = 1
    variable <- env.fact[k]
    temp.data.env <- data.env[, c("X","Y", env.fact[k])]
    no.rows <- nrow(temp.data.env)
    no.na <- sum(is.na(temp.data.env[,variable]))
    #explanation <- env.explan[which(env.explan$column.name == variable), "explanation"]
    
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill="#E8E8E8", color="black")
    g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
    g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
    g <- g + geom_point(data = temp.data.env, aes(x=X, y=Y, color= temp.data.env[, variable]), size= 3, alpha= 0.6)
    
    
      # Set up scales
      k.min <- round(min(temp.data.env[,variable], na.rm = T), 1)
      k.max <- round(max(temp.data.env[,variable], na.rm = T), 1)
      k.int <- (k.max - k.min)/5 # ; if.int <- round(if.int)
      
      g <- g + scale_colour_gradient2(name = variable, 
                                      #low = "coral", # useless
                                      high = "firebrick3",
                                      space = "Lab", na.value = "grey20", guide = "colourbar")
      # g <- g + guides(color = guide_bins(axis = T, show.limits = T, override.aes = list(alpha = 1)))
    
    
    g <- g + theme_void(base_size = 18)
    g <- g + theme(# plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(colour="transparent"),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
    g <- g + labs(title = paste("Geographic distribution of",variable), 
                  subtitle = paste0(no.na, " NAs out of ", no.rows, " samples \n"), colour = variable) + scale_fill_gradient2(midpoint = "value of midpoint", low = "blue", mid = "white", high = "blue", 
                                                                                                                              limits = c(min,max))
    
    
    print(g)
    
    p <- ggplot(data = temp.data.env, aes(temp.data.env[, variable]))
    
    
      
      q <- p + geom_histogram(col="grey20", fill="firebrick3", alpha = 0.2,
                              binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3)))
      q <- q + labs(title=paste("Histogram for", variable), x=variable, y="Frequency")
      
      r <- p + geom_density(fill="firebrick3", alpha = 0.2)
      r <- r + labs(x=variable, col="grey20", y="Density")
      
      p <- ggarrange(q,r, ncol = 2)

    
    print(p)
    # print(ggarrange(g,p, nrow =2))

  
}
}
