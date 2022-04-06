# Code that was in the main script for analysis
# ECR 18.03.2022

# Look into the effect of data leakage ####

# save intermediate output to compare dl no dl
# saveRDS(df.pred.perf, file = paste0(dir.models.output, info.file.name, "df_perf_.rds"), version = 2)
#df.pred.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_7models_59taxa_CV_no_DL_df_perf_.rds"))
#df.pred.perf.dl <- readRDS(file = paste0(dir.models.output, "All_7models_59taxa_CV_DL_df_perf_.rds"))

# df.pred.perf.no.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_no_DL_df_perf_.rds"))
# df.pred.perf.dl <- readRDS(file = paste0(dir.models.output, "All_9models_126taxa_CV_DL_df_perf_.rds"))
#
# df.pred.perf.no.dl$DL <- F
# df.pred.perf.dl$DL <- T
# df.pred.perf.dl.comb <- rbind(df.pred.perf.no.dl,df.pred.perf.dl)
#
# #remove infite vlaues to take mean later
# df.pred.perf.dl.comb$rf[!is.finite(df.pred.perf.dl.comb$rf)] <- NA
# df.pred.perf.dl.comb$gamSpline[!is.finite(df.pred.perf.dl.comb$gamSpline)] <- NA
# df.pred.perf.dl.comb$glm[!is.finite(df.pred.perf.dl.comb$glm)] <- NA
#
# colnames(df.pred.perf.dl.comb)
#
# df.pred.perf.dl.comb.table <- as.data.frame(df.pred.perf.dl.comb %>% group_by(DL) %>%
#                                       summarise(glm.mean = mean(glm, na.rm = T),
#                                       glm.sd = sd(glm, na.rm = T),
# 
#                                       gamLoess.mean = mean(gamLoess,na.rm = T),
#                                       gamLoess.sd = sd(gamLoess,na.rm = T),
# 
#                                       svmRadial.mean = mean(svmRadial, na.rm = T),
#                                       svmRadial.sd = sd(svmRadial, na.rm = T),
# 
#                                       rf.mean = mean(rf, na.rm = T),
#                                       rf.sd = sd(rf, na.rm = T),
# 
#                                       CF0.mean = mean(CF0, na.rm = T),
#                                       CF0.sd = sd(CF0, na.rm = T),
# 
#                                       UF0.mean = mean(UF0, na.rm = T),
#                                       UF0.sd = sd(UF0, na.rm = T),
# 
# 
#                                       ANN_3L32UleakyreluFCT50epo.mean = mean(ANN_3L32UleakyreluFCT50epo, na.rm = T),
#                                       ANN_3L32UleakyreluFCT50epo.sd = sd(ANN_3L32UleakyreluFCT50epo, na.rm = T)
# 
#                                       #ANN3L64U.mean = mean(ANN3L64U, na.rm = T),
#                                       #ANN3L64U.sd = sd(ANN3L64U, na.rm = T),
# 
#                                       #ANN5L32U.mean = mean(ANN5L32U, na.rm = T),
#                                       #ANN5L32U.sd = sd(ANN5L32U, na.rm = T)
# 
#                                        ))

#saveRDS(df.pred.perf.dl.comb.table, file = paste0(dir.plots.output, "Table_means_dl.rds"), version = 2)
#df.pred.perf.dl.comb.table <- readRDS(file = paste0(dir.plots.output, "Table_means_dl.rds"))
# diff.means.dl <- df.pred.perf.dl.comb.table[1,c("glm.mean", "gamLoess.mean","svmRadial.mean",
#                          "rf.mean", "CF0.mean", "UF0.mean", "ANN_3L32UleakyreluFCT50epo.mean")] -
#   df.pred.perf.dl.comb.table[2,c("glm.mean", "gamLoess.mean","svmRadial.mean",
#                             "rf.mean", "CF0.mean", "UF0.mean", "ANN_3L32UleakyreluFCT50epo.mean")]
# 
# diff.means.dl <- df.pred.perf.dl.comb.table[1,c("glm.mean", "CF0.mean", "UF0.mean")] -
#   df.pred.perf.dl.comb.table[2,c("glm.mean", "CF0.mean", "UF0.mean")]
# 
# 
# 
# dloui <- subset(df.pred.perf.dl.comb, DL == T)
# dloui <- dloui[,1:7]
# dlno <- subset(df.pred.perf.dl.comb, DL == F)
# dlno <- dlno[,1:7]
# comb <- cbind(dlno,dloui)
# comb2 <- dloui[,1:7] - dlno[,1:7]
# comb2 <- cbind(comb2,dlno[,8])
# colnames(comb2)[8] <- "Taxa"
# 
# comb3 <- comb2
# mean_difference <- colMeans(comb2[1:7])
# sd_difference <- colSds(as.matrix(comb2[1:7]))
# names(sd_difference) <- c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN")
# 
# table2 <- comb3 %>% gt()
# 
# tmp.table <- comb2[c("Taxa", "glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN_3L32UleakyreluFCT50epo")]
# colnames(tmp.table) <- c("Taxa", "glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN")
# #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
# table1 <- tmp.table %>% gt() %>%
#   tab_header(
#     title = md("**Difference in predictive performance between dl and no dl**") # make bold title
#   ) %>%
#   fmt_number(
#     columns = c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN"), # round numbers
#     decimals = 3
#   ) %>% # remove uneccessary black lines
#   tab_options(
#     table.border.top.color = "white",
#     heading.border.bottom.color = "black",
#     row_group.border.top.color = "black",
#     row_group.border.bottom.color = "white",
#     #stub.border.color = "transparent",
#     table.border.bottom.color = "white",
#     column_labels.border.top.color = "black",
#     column_labels.border.bottom.color = "black",
#     table_body.border.bottom.color = "black",
#     table_body.hlines.color = "white")
# 
# table2 <- comb3 %>% gt() %>%
#   tab_header(
#     title = md("**Difference in predictive performance between dl and no dl**") # make bold title
#   ) %>%
#   fmt_number(
#     columns = c("glm", "UF0","CF0","gamLoess", "svmRadial", "rf", "ANN_3L32UleakyreluFCT50epo"), # round numbers
#     decimals = 3
#   ) %>% # remove uneccessary black lines
#   tab_options(
#     table.border.top.color = "white",
#     heading.border.bottom.color = "black",
#     row_group.border.top.color = "black",
#     row_group.border.bottom.color = "white",
#     #stub.border.color = "transparent",
#     table.border.bottom.color = "white",
#     column_labels.border.top.color = "black",
#     column_labels.border.bottom.color = "black",
#     table_body.border.bottom.color = "black",
#     table_body.hlines.color = "white")
# 

# # Stat model traceplots ####
# res <- stat.outputs[[1]][[1]][[1]]
# res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
# 
# print(traceplot(res,pars=c(names(res)[134:162],"lp__")))
# print(traceplot(res))

# Compare models # changed by ecr 23.03.2022
model.comparison <- function(df.merged.perf, list.models, CV, extrapol, select.taxa){
  
  list.models.temp <- list.models
  list.models <- c("#000000" = "Null_model", list.models)
  
  # Make a vector of colors
  col.vect <- names(list.models)
  names(col.vect) <- list.models
  
  no.taxa <- nrow(df.merged.perf)
  
  title <- c("Models comparison in quality of fit", "Models comparison in predictive performance")
  
  if(!CV & !extrapol){
    title <- title[1]
  }
  
  subtitle <- paste("For", no.taxa, "taxa")
  
  # ECR: To be completed if CV = F
  
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
  } else {
    list.plots <- list.plots.temp
  }
  
  return(list.plots)
  
}




# Performance vs hyperparameters ####

# Compute plots
list.plots <- plot.perf.hyperparam(outputs = outputs, 
                                   list.algo = list.algo[4], # GLM algo doesn't have hyperparameters
                                   list.taxa = list.taxa)
# Print a pdf file
name <- "PerfvsHyperparam"
file.name <- paste0(name, ".pdf")
if( file.exists(paste0(dir.plots.output, info.file.name, file.name)) == F ){ # produce pdf only if it doesn't exist yet (takes too much time)
  print.pdf.plots(list.plots = list.plots, width = 9, height = 9, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
}
# # Print a jpeg file
# file.name <- paste0(name, ".jpg")
# jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
# print(list.plots[[1]])
# dev.off()

# Variables importance ####

list.plots <- plot.varimp(outputs = outputs, list.algo = list.algo, list.taxa = list.taxa)

name <- "VarImp"
file.name <- paste0(name, ".pdf")
if( file.exists(paste0(dir.plots.output, info.file.name, file.name)) == F ){ # produce pdf only if it doesn't exist yet (takes too much time)
  print.pdf.plots(list.plots = list.plots, width = 10, height = 10, dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name)
}
# # Print a jpeg file
# file.name <- paste0(name, ".jpg")
# jpeg(paste0(dir.plots.output,"JPEG/",info.file.name,file.name))
# print(list.plots[[1]])
# dev.off()


# Functions for plotting ####

# Plot PDP

old.plot.pdp <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
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

# Plot ICE

old.plot.ice <- function(outputs, algo = "all", list.algo, list.taxa, env.fact){
  
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


