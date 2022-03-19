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


