stat_mod_cv <- function (data.splits, CV, ODG, comm.corr, sampsize, n.chain){
    
    ##To test
    #comm.corr <- T
    #CV <- F
    #sampsize        <- 11 #10000
    #n.chain         <- 1 #2
  
    # global parameters ####
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    options(mc.cores = parallel::detectCores()) #this is to run chains in parallel, comment out if there is an error in joining the chains
    rstan_options(auto_write = TRUE)
    
    model.name <- ifelse(comm.corr, "CF0", "UF0")
    site.effects    <- F
    n.latent        <- 0
    lat.site        <- F
    
    frac.occ.lat    <- 0   # taxa only considered for latent variable if relative frequency of occurrence > frac.occ.lat
    generate.res    <- F
    thin            <- 5 #
    prob.defpri     <- 0.02
    thresh.sig      <- 1
    fact.sd         <- 1  
    max.warmup      <- 1000 # max warm up iterations for runs with a lot of iterations
    
    # data.splits <- centered.splits[[1]] # to test
    # data.splits <- splits[[1]] # to test
    # data.splits <- standardized.data[[1]]
    output <- list("deviance" = tibble(), "probability" = tibble(), "parameters" = tibble()) #this is where the data is gathered in the end for the return
    training.data <- data.splits[[1]]
    
    inv.names <- colnames(select(training.data, SiteId, SampId, contains("Occurrence.")))
    env.names <- colnames(select(training.data, colnames(training.data), -contains("Occurrence.")))
    
    env.cond <- training.data[,env.names]
    occur.taxa <- training.data[,inv.names]
    # comm.corr = comm.corr
    inf.fact <- colnames(select(env.cond, -SiteId, -SampId, -X, -Y))
    # 
    # # join the environmental conditions to the occurrence data
    # env.cond <- left_join(occur.taxa[, c("SiteId", "SampId")], env.cond.orig, by = c("SiteId", "SampId"))
    # 
    # # Pre-process data ####
    # # drop rows with incomplete influence factors:
    # ind <- !apply(is.na(env.cond[,inf.fact]),1,FUN=any)
    # ind <- ifelse(is.na(ind),FALSE,ind)
    # occur.taxa <- occur.taxa[ind, ]
    # env.cond   <- env.cond[ind, c("SiteId", "SampId", inf.fact)]
    # print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))
    # 
    sites <- occur.taxa$SiteId
    samples <- occur.taxa$SampId
    
    n.sites <- length(unique(sites))
    n.samples <- length(samples)
    
    occur.taxa$SiteId <- NULL
    occur.taxa$SampId <- NULL
    # occur.taxa <- as.data.frame(occur.taxa)
    # # drop TAXA without observations at the selected sites:
    # ind <- apply(occur.taxa,2,sum, na.rm = TRUE) > 0
    # occur.taxa <- occur.taxa[, ind]
    n.taxa <- ncol(occur.taxa)
    # 
    # unique.sites <- unique(sites)
    # siteIND <- match(sites, unique.sites)
    # env.cond2 <- env.cond
    # 
    # #if center is true substract the mean of each predictor, check if its divded by sd, I added the division by sd
    # if(center){
    #   mean.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
    #     mean(k, na.rm = TRUE)
    #   })
    #   sd.env.cond <- apply(env.cond[, !(colnames(env.cond) %in% c("SiteId", "SampId"))], 2, function(k){
    #     sd(k, na.rm = TRUE)
    #   })
    #   for(i in 1:length(env.cond[!(colnames(env.cond) %in% c("SiteId", "SampId"))])){
    #     #i = 6
    #     #message(i)
    #     env.cond[i+2] <- as.matrix(env.cond[i+2]) - mean.env.cond[i]
    #   }
    #   for(i in 1:length(env.cond[!(colnames(env.cond) %in% c("SiteId", "SampId"))])){
    #     #i = 6
    #     env.cond[i+2] <- as.matrix(env.cond[i+2]) / sd.env.cond[i]
    #   }
    # }
    # 
    # 
    # replace NA with -1
    occur.taxa.na.encod <- occur.taxa
    for(j in 1:ncol(occur.taxa.na.encod)) occur.taxa.na.encod[,j] <- ifelse(is.na(occur.taxa.na.encod[,j]),
                                                                            -1, occur.taxa.na.encod[,j])
    
    #comment this better
    ind <- 1:nrow(occur.taxa)
    set.ident <- rep(TRUE,length(ind))
    set.valid <- set.ident
    
    # Run GLMs ####
    # GLM analysis of occurrence data (needed for initial conditions of Markov chains):
    
    res.glm <- data.frame(matrix(NA,nrow=ncol(occur.taxa),ncol=2+1+length(inf.fact),
                                 dimnames=list(colnames(occur.taxa),c("present","absent","alpha",inf.fact))))
    for ( j in 1:ncol(occur.taxa) )
    {
        taxon <- colnames(occur.taxa)[j]
        n <- sum(occur.taxa[set.ident,j], na.rm = TRUE) # loop through rows and taxa
        res.glm[j,1] <- n
        res.glm[j,2] <- sum(set.ident)-n
        if ( n > 5 )
        {
            occur <- occur.taxa[set.ident,j]
            res.glm.j <- glm(formula(paste("occur ~",paste(inf.fact,collapse="+"))),
                             family="binomial",
                             data=env.cond[set.ident,])
            est <- summary(res.glm.j)$coefficients[,"Estimate"]
            sd  <- summary(res.glm.j)$coefficients[,"Std. Error"]
            est <- ifelse(abs(est)>thresh.sig*sd,est,0)
            res.glm[j,3:ncol(res.glm)] <- as.numeric(est)
            res.glm[j,3] <- as.numeric(summary(res.glm.j)$coefficients[,"Estimate"][1])
        }
        else
        {
            if ( n > 0 )
            {
                occur <- occur.taxa[set.ident,j]
                res.glm.j <- glm(formula(paste("occur ~","1")),
                                 family="binomial",
                                 data=env.cond[set.ident,])
                res.glm[j,3] <- as.numeric(summary(res.glm.j)$coefficients[,"Estimate"][1])
                res.glm[j,4:ncol(res.glm)] <- 0
            }
            else
            {
                res.glm[j,3] <- -Inf
                res.glm[j,4:ncol(res.glm)] <- 0
            }
        }
    }
    
    # Model ####
    
    submodel <- ifelse(comm.corr,"corr","nocorr")
    submodel <- paste(submodel,ifelse(site.effects,"site","nosite"),sep="_")
    submodel <- paste(submodel,ifelse(n.latent>0,paste(n.latent,"latvar",sep=""),"nolatvar"),sep="_")
    if ( n.latent > 0 ) submodel <- paste(submodel,ifelse(lat.site,"site","samp"),sep="")
    file.model <- paste("Inv_JSDM_",submodel,".stan",sep="")
    
    ### Define priors ####
    
    logit <- function(p) { return(log(p/(1-p))) }
    z.range.max <- logit(0.95)-logit(0.05)
    
    mu.alpha.comm.pripar    <- 0
    sigma.alpha.comm.pripar <- length(inf.fact)*z.range.max/4
    sigma.beta.comm.pripar <- z.range.max/(0.5*apply(env.cond[inf.fact],2,function(k){
        max(k)-min(k)
    }))
    mu.beta.comm.pripar <- rep(0,length(sigma.beta.comm.pripar))
    names(mu.beta.comm.pripar) <- names(sigma.beta.comm.pripar)
    
    # compilation of input data:
    data <- list(n_sites                 = n.sites,
                 n_samples               = n.samples,
                 n_taxa                  = n.taxa,
                 n_pred                  = length(inf.fact),
                 mu_alpha_comm_pripar    = mu.alpha.comm.pripar,
                 sigma_alpha_comm_pripar = sigma.alpha.comm.pripar,
                 mu_beta_comm_pripar     = mu.beta.comm.pripar,
                 sigma_beta_comm_pripar  = sigma.beta.comm.pripar,
                 fact_sd                 = fact.sd,
                 x                       = as.matrix(env.cond[set.ident,inf.fact]),
                 y                       = as.matrix(occur.taxa.na.encod[set.ident,]))
    
    if ( site.effects ) {
        data$sigma_gamma_pripar <- 2*logit(1-prob.defpri)/10
        data$siteIND <- siteIND
    }
    if ( n.latent > 0 )
    {
        data$sigma_beta_lat <- z.range.max/4
        if ( n.latent > 1 )
        {
            data$n_latent     <- n.latent
        }
    }
    
    ### Initialize chains ####
    # definition of (critical) starting points of Markov chains:
    init <- list()
    for ( i in 1:n.chain )
    { 
        init[[i]] <- list()
        init.alpha <- res.glm[,"alpha"]
        init.beta  <- t(res.glm[,inf.fact])
        # cut extreme values:
        min.max <- quantile(init.alpha[init.alpha!=-Inf],probs=c(0.2,0.8))
        init.alpha <- ifelse(init.alpha<min.max[1],min.max[1],ifelse(init.alpha>min.max[2],min.max[2],init.alpha))
        init.alpha <- init.alpha*(1+0.1*rnorm(length(init.alpha))) 
        init[[i]][["alpha_taxa"]] <- init.alpha
        for ( k in 1:length(inf.fact) )
        {
            min.max <- quantile(init.beta[k,],probs=c(0.2,0.8))
            init.beta[k,] <- ifelse(init.beta[k,]<min.max[1],min.max[1],ifelse(init.beta[k,]>min.max[2],min.max[2],init.beta[k,]))
            init.beta[k,] <- init.beta[k,]*(1+0.1*rnorm(ncol(init.beta)))
        }
        init[[i]][["beta_taxa"]] <- init.beta
        init[[i]][["mu_alpha_comm"]] <- mean(init[[i]]$alpha_taxa)
        init[[i]][["mu_beta_comm"]]  <- apply(init[[i]]$beta_taxa,1,mean)
        s.pripar = sqrt(log(1+fact.sd^2))
        init[[i]][["sigma_alpha_comm"]] <- sigma.alpha.comm.pripar
        init[[i]][["sigma_beta_comm"]]  <- sigma.beta.comm.pripar
        
        if ( site.effects )
        {
            init[[i]][["sigma_gamma"]] <- data$sigma_gamma_pripar
            init[[i]][["gamma_site"]] <- rep(0,n.sites)
        }
    }
    
    
    ### Model definition ####
    
    # data
    # ----
    
    model.code <- paste("// Model:    Logistic, hierarchical JSDM.\n",
                        "// Submodel: ",submodel,"\n",
                        "// -----------------------------------------------------------------------\n\n",
                        "// data:\n\n",
                        "data {\n",
                        "  int                              n_taxa;\n",
                        "  int                              n_sites;\n",
                        "  int                              n_samples;\n",
                        "  int                              n_pred;\n",
                        "  real                             mu_alpha_comm_pripar;\n",
                        "  real                             sigma_alpha_comm_pripar;\n",
                        "  vector[n_pred]                   mu_beta_comm_pripar;\n",
                        "  vector[n_pred]                   sigma_beta_comm_pripar;\n",
                        "  real                             fact_sd;\n",
                        sep = "")
    if ( site.effects )
    {
        model.code <- paste(model.code,
                            "  real                             sigma_gamma_pripar;\n",
                            "  int                              siteIND[n_samples];\n",
                            sep="")
    }
    if ( n.latent > 0 )
    {
        model.code <- paste(model.code,
                            "  real                             sigma_beta_lat;\n",
                            sep="")
        if ( n.latent > 1 )
        {
            model.code <- paste(model.code,
                                "  int                            n_latent;\n",
                                sep="")
        }
    }                      
    model.code <- paste(model.code,
                        "  matrix[n_samples,n_pred]         x;\n",
                        "  int<lower=-1,upper=1>            y[n_samples,n_taxa];\n",
                        "}\n\n", 
                        sep="")
    
    # parameters
    # ----------
    
    model.code <- paste(model.code,
                        "// parameters:\n\n",
                        "parameters {\n",
                        "  real                             mu_alpha_comm;\n",
                        "  real<lower=0>                    sigma_alpha_comm;\n",
                        "  vector[n_pred]                   mu_beta_comm;\n",
                        "  vector<lower=0>[n_pred]          sigma_beta_comm;\n",
                        "  vector[n_taxa]                   alpha_taxa;\n",
                        "  matrix[n_pred,n_taxa]            beta_taxa;\n",
                        sep="")
    if (comm.corr)
    {
        model.code <- paste(model.code,
                            "//corr_matrix[n_pred]              corr_comm;\n",
                            "  cholesky_factor_corr[n_pred]     corrfact_comm;\n",
                            sep="")
    }
    if ( site.effects )
    {
        model.code <- paste(model.code,
                            "  real<lower=0>                    sigma_gamma;\n",
                            "  vector[n_sites]                  gamma_site;\n",
                            sep="")
    }
    if ( n.latent == 1 )
    {
        if ( lat.site )
        {
            model.code <- paste(model.code,
                                "  vector[n_sites]                  x_lat;\n",
                                sep="")
        }
        else
        {
            model.code <- paste(model.code,
                                "  vector[n_samples]                x_lat;\n",
                                sep="")
        }
        model.code <- paste(model.code,
                            "  vector[n_taxa]                   beta_lat;\n",
                            sep="")
    }  
    if ( n.latent > 1 )
    {
        if ( lat.site )
        {
            model.code <- paste(model.code,
                                "  matrix[n_sites,n_latent]         x_lat;\n",
                                sep="")
        }
        else
        {
            model.code <- paste(model.code,
                                "  matrix[n_samples,n_latent]       x_lat;\n",
                                sep="")
        }
        model.code <- paste(model.code,
                            "  matrix[n_latent,n_taxa]          beta_lat;\n",
                            sep="")
    }
    model.code <- paste(model.code,
                        "}\n\n", 
                        sep="")
    
    # local variables
    # ---------------
    
    model.code <- paste(model.code,
                        "// model definition:\n\n",
                        "model {\n\n",
                        "  // local variables:\n\n",
                        "  real z[n_samples,n_taxa];\n",
                        "  real p[n_samples,n_taxa];\n",
                        "  real s_pripar;\n\n",
                        sep="")
    
    # root nodes
    # ----------
    
    model.code <- paste(model.code,
                        "  // root nodes:\n\n",
                        "  mu_alpha_comm ~ normal(mu_alpha_comm_pripar,sigma_alpha_comm_pripar);\n",
                        "  s_pripar = sqrt(log(1+fact_sd^2));\n",
                        "  sigma_alpha_comm   ~ lognormal(log(sigma_alpha_comm_pripar)-0.5*s_pripar^2,s_pripar);\n",
                        "  for ( k in 1:n_pred ) {\n",
                        "    mu_beta_comm[k] ~ normal(mu_beta_comm_pripar[k],sigma_beta_comm_pripar[k]);\n",
                        "    s_pripar = sqrt(log(1+fact_sd^2));\n",
                        "    sigma_beta_comm[k]   ~ lognormal(log(sigma_beta_comm_pripar[k])-0.5*s_pripar^2,s_pripar);\n",
                        "  }\n",
                        sep="")
    if ( site.effects )
    {
        model.code <- paste(model.code,
                            "  s_pripar = sqrt(log(1+fact_sd^2));\n",
                            "  sigma_gamma ~ lognormal(log(sigma_gamma_pripar)-0.5*s_pripar^2,s_pripar);\n",
                            sep="")
    }
    if ( comm.corr )
    {
        model.code <- paste(model.code,
                            "//corr_comm ~ lkj_corr(2);\n",
                            "  corrfact_comm ~ lkj_corr_cholesky(2);\n",
                            sep="")
    }
    if ( n.latent == 1 )
    {
        if ( lat.site )
        {
            model.code <- paste(model.code,
                                "  for ( i in 1:n_sites ) {\n",
                                sep="")
        }
        else
        {
            model.code <- paste(model.code,
                                "  for ( i in 1:n_samples ) {\n",
                                sep="")
        }
        model.code <- paste(model.code,
                            "    x_lat[i] ~ normal(0,1);\n",
                            "  }\n",
                            "  for ( j in 1:n_taxa ) {\n",
                            "    beta_lat[j] ~ normal(0,sigma_beta_lat);\n",
                            "  }\n",
                            sep="")
    }
    if ( n.latent > 1 )
    {
        model.code <- paste(model.code,
                            "  for ( k in 1:n_latent ) {\n",
                            sep="")
        if ( lat.site )
        {
            model.code <- paste(model.code,
                                "    for ( i in 1:n_sites ) {\n",
                                sep="")
        }
        else
        {
            model.code <- paste(model.code,
                                "    for ( i in 1:n_samples ) {\n",
                                sep="")
        }
        model.code <- paste(model.code,
                            "      x_lat[i,k] ~ normal(0,1);\n",
                            "    }\n",
                            "    for ( j in 1:n_taxa ) {\n",
                            "      beta_lat[k,j] ~ normal(0,sigma_beta_lat);\n",
                            "    }\n",
                            "  }\n",
                            sep="")
    }
    model.code <- paste(model.code,"\n",sep="")
    
    # intermediate nodes
    # ------------------
    
    model.code <- paste(model.code,
                        "  // intermediate nodes:\n\n",
                        "  for ( j in 1:n_taxa ) {\n",
                        "    alpha_taxa[j] ~ normal(mu_alpha_comm,sigma_alpha_comm);\n",
                        sep="")
    if ( !comm.corr )
    {
        model.code <- paste(model.code,
                            "    for ( k in 1:n_pred ) {\n",
                            "      beta_taxa[k,j] ~ normal(mu_beta_comm[k],sigma_beta_comm[k]);\n",
                            "    }\n",
                            sep="")
    } else {
        model.code <- paste(model.code,
                            "  //beta_taxa[,j] ~ multi_normal(mu_beta_comm,quad_form_diag(corr_comm,sigma_beta_comm));\n",
                            "    beta_taxa[,j] ~ multi_normal_cholesky(mu_beta_comm,diag_pre_multiply(sigma_beta_comm,corrfact_comm));\n",
                            sep="")
    }
    model.code <- paste(model.code,
                        "  }\n",
                        sep="")
    if ( site.effects )
    {
        model.code <- paste(model.code,
                            "  for ( i in 1:n_sites ) {\n",
                            "    gamma_site[i] ~ normal(0,sigma_gamma);\n",
                            "  }\n",
                            sep="")
    }
    model.code <- paste(model.code,
                        "  \n",
                        sep="")
    
    # end nodes
    # ---------
    
    model.code <- paste(model.code,
                        "  // end nodes:\n\n",
                        "  for ( i in 1:n_samples ) {\n",
                        "    for ( j in 1:n_taxa ) {\n",
                        "      z[i,j] = alpha_taxa[j];\n",
                        "      for ( k in 1:n_pred ) {\n",
                        "        z[i,j] = z[i,j] + beta_taxa[k,j]*x[i,k];\n",
                        "      }\n",
                        sep="")
    if ( site.effects )
    {
        model.code <- paste(model.code,
                            "      z[i,j] = z[i,j] + gamma_site[siteIND[i]];\n",        
                            sep="")
    }
    if ( n.latent == 1 )
    {
        if ( lat.site )
        {
            model.code <- paste(model.code,
                                "      z[i,j] = z[i,j] + x_lat[siteIND[i]]*beta_lat[j];\n",        
                                sep="")
        }
        else
        {
            model.code <- paste(model.code,
                                "      z[i,j] = z[i,j] + x_lat[i]*beta_lat[j];\n",        
                                sep="")
        }
    }
    if ( n.latent > 1 )
    {
        model.code <- paste(model.code,
                            "      for ( k in 1:n_latent ) {\n",
                            "        z[i,j] = z[i,j] + x_lat[siteIND[i],k]*beta_lat[k,j];\n",
                            "      }\n",
                            sep="")
    }
    model.code <- paste(model.code,
                        "      p[i,j] = 1/(1+exp(-z[i,j]));\n",
                        "      if(y[i,j]>=0) {\n",
                        "        y[i,j] ~ bernoulli(p[i,j]);\n",
                        "      }\n",
                        "    }\n",
                        "  }\n",
                        "}\n\n",
                        sep="")
    
    # generated quantities:
    # ---------------------
    
    if ( generate.res )
    {
        
        model.code <- paste(model.code,
                            "// generated quantities:\n\n",
                            "generated quantities {\n\n",
                            "  // output variables:\n\n",
                            "  real z_mod[n_samples,n_taxa];\n",
                            "  real p_mod[n_samples,n_taxa];\n",
                            "  real deviance[n_taxa];\n",
                            "\n",
                            sep="")
        
        model.code <- paste(model.code,
                            "  // end nodes:\n\n",
                            "  for ( i in 1:n_samples ) {\n",
                            "    for ( j in 1:n_taxa ) {\n",
                            "      z_mod[i,j] = alpha_taxa[j];\n",
                            "      for ( k in 1:n_pred ) {\n",
                            "        z_mod[i,j] = z_mod[i,j] + beta_taxa[k,j]*x[i,k];\n",
                            "      }\n",
                            sep="")
        if ( site.effects )
        {
            model.code <- paste(model.code,
                                "      z_mod[i,j] = z_mod[i,j] + gamma_site[siteIND[i]];\n",        
                                sep="")
        }
        if ( n.latent == 1 )
        {
            model.code <- paste(model.code,
                                "      z_mod[i,j] = z_mod[i,j] + x_lat[siteIND[i]]*beta_lat[j];\n",        
                                sep="")
        }
        if ( n.latent > 1 )
        {
            model.code <- paste(model.code,
                                "      for ( k in 1:n_latent ) {\n",
                                "        z_mod[i,j] = z_mod[i,j] + x_lat[siteIND[i],k]*beta_lat[k,j];\n",
                                "      }\n",
                                sep="")
        }
        model.code <- paste(model.code,
                            "      p_mod[i,j] = 1/(1+exp(-z_mod[i,j]));\n",
                            "    }\n",
                            "  }\n",
                            "  for ( j in 1:n_taxa ) {\n",
                            "    deviance[j] = 0;\n",
                            "    for ( i in 1:n_samples ) {\n",
                            "      if(y[i,j]>=0) {\n",
                            "        deviance[j] = deviance[j] - 2*(y[i,j]*log(p_mod[i,j])+(1-y[i,j])*log(1-p_mod[i,j]));\n",
                            "      }\n",
                            "    }\n",
                            "  }\n",
                            "}\n\n",
                            sep="")
    }
    
    cat(model.code, file=file.model)
    
    # Run model ####
    # perform Bayesian inference:
    res <- stan(file.model,data=data,init=init,iter=sampsize,chains=n.chain,warmup=min(0.5*sampsize,max.warmup),thin=thin)
    #res2 <- stat.outputs[[1]][[1]]
    #res <- res2[[3]]
    
    res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
    #saveRDS(res.extracted, file = paste0(dir.output,"test2411.rds"))
    
    #res.extracted <- readRDS(paste0(dir.output,"test2411.rds"))
    # Name dimensions of parameters within stanfit object
    #dimnames(res.extracted[["beta_taxa"]]) <- list(1:dim(res.extracted[["beta_taxa"]][,,])[1], inf.fact, colnames(occur.taxa))
    colnames(res.extracted[["alpha_taxa"]]) <- colnames(occur.taxa)
    
    colnames(res.extracted[["mu_beta_comm"]]) <- inf.fact
    colnames(res.extracted[["sigma_beta_comm"]]) <- inf.fact
    
    # Extract inputs (x), observations (y), and parameters at maximum posterior
    x <- as.matrix(env.cond[,inf.fact])
    colnames(x) <- inf.fact
    y <- as.matrix(occur.taxa)
    ind.maxpost <- which.max(res.extracted[["lp__"]])
    mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
    sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
    mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
    sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
    alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
    beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
    
    # Name dimensions of maximum posterior community parameters
    names(mu.beta.comm.maxpost) <- inf.fact
    names(sigma.beta.comm.maxpost) <- inf.fact
    rownames(beta.taxa.maxpost) <- inf.fact
    
    ### Calibration results
    # Check if site effects AND latent variables are disabled
    z <- matrix(rep(alpha.taxa.maxpost,nrow(x)),nrow=nrow(x),byrow=TRUE) + 
        x%*%beta.taxa.maxpost
    
    p.maxpost <- 1/(1+exp(-z))
    if (CV == T | ODG == T){
        train.y <- occur.taxa
        train.n.present <- apply(train.y, 2, sum, na.rm = TRUE)
        train.n.samples <- apply(train.y, 2, function(j){sum(!is.na(j))})
        
        train.p.null <- apply(train.y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
        train.p.null <- ifelse(train.p.null==0,1e-4,train.p.null)
        train.p.null <- matrix(rep(train.p.null,nrow(train.y)),nrow=nrow(train.y),byrow=TRUE)
        
        # Calculate deviance over the samples
        train.deviance.resid <- sign(train.y-0.5)*sqrt(-2*(train.y*log(p.maxpost)+(1-train.y)*log(1-p.maxpost)))
        
        # Calculate the null ("primitive") deviance
        train.deviance.null <- -2*sum(train.y*log(train.p.null)+(1-train.y)*log(1-train.p.null), na.rm = T)
        
        # Calculate the sum of squared residual deviance for all taxa
        train.deviance.resid.taxa <- apply(train.deviance.resid^2,2,sum, na.rm=T) # NAs removed
        train.deviance.null.taxa <- -2*apply(train.y*log(train.p.null)+(1-train.y)*log(1-train.p.null),2,sum,na.rm=T)
        
        train.d <- 1-train.deviance.resid.taxa/train.deviance.null.taxa ###
        
        # Standardize deviance by number of observations in training data
        train.std.deviance <- train.deviance.resid.taxa/train.n.samples
        
        # tidy deviance data
        train.deviance <- tibble(null.deviance = train.deviance.null.taxa, 
                                 residual.deviance = train.deviance.resid.taxa, 
                                 std.deviance = train.std.deviance, 
                                 D2 = train.d,
                                 Type = "Training",
                                 Model = model.name,
                                 stringsAsFactors = F)
        
        train.deviance$Taxon <- names(train.std.deviance)
        train.deviance$n.samples <- train.n.samples[train.deviance$Taxon]
        train.deviance$n.present <- train.n.present[train.deviance$Taxon]
        
        train.deviance <- train.deviance[, c("Taxon", "Type", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
        
        # tidy probability data
        train.obs <- as_tibble(occur.taxa)
        train.obs$SiteId <- sites
        train.obs$SampId <- samples
        train.obs <- gather(train.obs, Taxon, Obs, -SiteId, -SampId)
        
        train.p <- as_tibble(p.maxpost, stringsAsFactors = FALSE)
        colnames(train.p) <- colnames(occur.taxa)
        train.p$SiteId <- sites
        train.p$SampId <- samples
        train.p <- gather(train.p, Taxon, Pred, -SiteId, -SampId)
        train.p <- left_join(train.p, train.obs, by = c("SiteId", "SampId", "Taxon"))
        train.p$Type <- "Training"
        train.p$Model <- model.name

        ### TEST results

        testing.data <- data.splits[[2]]
        test.predictors <- testing.data[,env.names]
        test.y <- testing.data[,inv.names]
        # join the environmental conditions to the occurrence data
        
        
        # Drop rows in y.test and predictors.test where the latter has any NAs
        # ind <- !apply(is.na(test.predictors[,inf.fact]),1,FUN=any)
        # ind <- ifelse(is.na(ind),FALSE,ind)
        # test.predictors <- test.predictors[ind, ]
        # test.y <- test.y[ind, ]
        # 
        # if(center){
        #   
        #   for(i in 1:length(test.predictors[!(colnames(test.predictors) %in% c("SiteId", "SampId"))])){
        #     #i = 6
        #     #message(i)
        #     test.predictors[i+2] <- as.matrix(test.predictors[i+2]) - mean.env.cond[i]
        #   }
        #   for(i in 1:length(test.predictors[!(colnames(test.predictors) %in% c("SiteId", "SampId"))])){
        #     #i = 6
        #     test.predictors[i+2] <- as.matrix(test.predictors[i+2] / sd.env.cond[i])
        #   }
        # }
        
        
        # Keep the test sites
        test.sites <- test.predictors$SiteId
        test.samples <- test.predictors$SampId
        names(test.sites) <- test.samples
        
        # Drop the test sites and keep the inputs
        test.inputs <- as.matrix(test.predictors[, inf.fact])
        
        z <- matrix(rep(alpha.taxa.maxpost,nrow(test.inputs)),nrow=nrow(test.inputs),byrow=TRUE) + 
            test.inputs%*%beta.taxa.maxpost 
        #
        p.maxpost <- 1/(1+exp(-z))
        
        # Keep only taxa present in training data ## JW:IS this really correct?
        test.y <- test.y[, c("SiteId", "SampId", names(train.n.present))]
        test.n.present <- apply(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))], 2, sum, na.rm = TRUE)
        test.n.samples <- apply(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))], 2, function(j){sum(!is.na(j))})
        
        # Calculate null predicted probabilities
        test.p.null <- apply(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))], 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
        test.p.null <- ifelse(test.p.null==0,1e-4,test.p.null)
        test.p.null <- matrix(rep(test.p.null,nrow(test.y)),nrow=nrow(test.y),byrow=TRUE)
        
        # Calculate null and residual deviances
        test.deviance.resid <- sign(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))]-0.5)*sqrt(-2*(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))]*log(p.maxpost)+(1-test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))])*log(1-p.maxpost)))
        # deviance.maxpost.test   <- sum(test.deviance.resid^2, na.rm = T)
        test.deviance.null <- -2*sum(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))]*log(test.p.null)+(1-test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))])*log(1-test.p.null), na.rm = T)
        
        test.deviance.resid.taxa <- apply(test.deviance.resid^2,2,sum, na.rm=T) # NAs removed
        test.deviance.null.taxa <- -2*apply(test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))]*log(test.p.null)+(1-test.y[!(colnames(test.y) %in% c("SiteId", "SampId"))])*log(1-test.p.null),2,sum,na.rm=T)
        
        # Calculate explanatory power: D2
        test.d <- 1-test.deviance.resid.taxa/test.deviance.null.taxa
        
        # Relative deviance
        test.std.deviance <- test.deviance.resid.taxa/test.n.samples
        
        # Tidy deviance data
        test.deviance <- tibble(null.deviance = test.deviance.null.taxa, 
                                residual.deviance = test.deviance.resid.taxa, 
                                std.deviance = test.std.deviance, 
                                D2 = test.d, 
                                Type = "Testing",
                                Model = model.name,
                                stringsAsFactors = F)
        
        test.deviance$Taxon <- names(test.std.deviance)
        test.deviance$n.samples <- test.n.samples[test.deviance$Taxon]
        test.deviance$n.present <- test.n.present[test.deviance$Taxon]
        
        test.deviance <- test.deviance[, c("Taxon", "Type", "Model", "null.deviance", "residual.deviance", "std.deviance", "D2", "n.samples", "n.present")]
        
        # tidy the probability data
        test.y <- as_tibble(test.y)
        test.y <- gather(test.y, Taxon, Obs, -SiteId, -SampId)
        
        test.p <- as_tibble(p.maxpost, stringsAsFactors = FALSE)
        colnames(test.p) <- colnames(occur.taxa)
        test.p$SiteId <- test.sites
        test.p$SampId <- test.samples
        test.p <- gather(test.p, Taxon, Pred, -SiteId, -SampId)
        test.p <- left_join(test.p, test.y, by = c("SiteId", "SampId", "Taxon"))
        
        test.p$Type <- "Testing"
        test.p$Model <- model.name

        ### Bind the k-fold results
        output$deviance <- bind_rows(output$deviance, train.deviance)
        output$deviance <- bind_rows(output$deviance, test.deviance)
        output$probability <- bind_rows(output$probability, train.p, test.p)
        return(list(res,output))
        
    }
    else{
        # Prepare tidy data from maximum posterior beta_taxa
        beta <- t(beta.taxa.maxpost)
        beta <- as_tibble(beta, stringsAsFactors = F)
        beta$Taxon <- colnames(occur.taxa)
        beta <- gather(beta, Variable, Parameter, -Taxon)
        
        p.primitive <- apply(y, 2, function(j){sum(j, na.rm=TRUE)/sum(!is.na(j))})
        p.primitive <- ifelse(p.primitive==0,1e-4,p.primitive)
        p.primitive <- matrix(rep(p.primitive,nrow(y)),nrow=nrow(y),byrow=TRUE)
        dev.resid <- sign(y-0.5)*sqrt(-2*(y*log(p.maxpost)+(1-y)*log(1-p.maxpost)))
        deviance.maxpost   <- sum(dev.resid^2, na.rm = T)
        deviance.primitive <- -2*sum(y*log(p.primitive)+(1-y)*log(1-p.primitive), na.rm = T)
        
        ### Deviance
        # Residual deviance
        deviance.maxpost.taxa <- apply(dev.resid^2,2,sum, na.rm=T) # NAs removed
        # Null deviance
        deviance.primitive.taxa <- -2*apply(y*log(p.primitive)+(1-y)*log(1-p.primitive),2,sum,na.rm=T)
        # D^2 statistic
        deviance.fit <- 1-deviance.maxpost.taxa/deviance.primitive.taxa ### D2
        # Number of samples
        n.samples <- apply(y, 2, function(j){ # Number of presence-absence observations
            sum(!is.na(j))
        })
        n.present <- apply(y, 2, sum, na.rm = TRUE)
        
        deviance <- tibble(Taxon = names(deviance.fit), Model = model.name, null.deviance = deviance.primitive.taxa, residual.deviance = deviance.maxpost.taxa, std.deviance = deviance.maxpost.taxa/n.samples, D2 = deviance.fit, n.samples = n.samples[names(deviance.fit)], n.present = n.present[names(deviance.fit)])
        
        ### Prepare tidy data for output
        # Melt the occurrence data into columns of Taxon and Obs
        obs <- as_tibble(occur.taxa)
        obs$SiteId <- sites
        obs$SampId <- samples
        obs <- gather(obs, Taxon, Obs, -SiteId, -SampId)
        
        # Melt the probabilities into columns of Taxon and Pred
        prob <- as_tibble(p.maxpost)
        colnames(prob) <- colnames(occur.taxa)
        prob$SiteId <- sites
        prob$SampId <- samples
        prob <- gather(prob, Taxon, Pred, -SiteId, -SampId)
        
        # Combine occurrence/probability data and join the site coordinates
        probability <- left_join(obs, prob, by = c("SiteId", "SampId", "Taxon"))
        probability <- na.omit(probability)
        rm(obs, prob)
        
        #JW we could add Tjur metric here
        
        ### Prepare output
        probability$Model <- model.name
        deviance$Model <- model.name
        beta$Model <- model.name
        
        output$deviance <- bind_rows(output$deviance, deviance)
        output$deviance.residual   <- deviance.maxpost
        output$deviance.null <- deviance.primitive
        
        output$probability <- bind_rows(output$probability, probability)
        output$parameters <- bind_rows(output$parameters, beta) #not quite sure about beta here
        # Store the site, samples, community, and input data
        output$sites <- sites
        output$samples <- samples
        output$occur.taxa <- occur.taxa
        output$inf.fact <- inf.fact
        output$env.cond <- env.cond

        # Store priors for community parameters
        output$mu.alpha.comm.pripar <- data$mu_alpha_comm_pripar
        output$sigma.alpha.comm.pripar <- data$sigma_alpha_comm_pripar
        output$mu.beta.comm.pripar <- data$mu_beta_comm_pripar
        output$sigma.beta.comm.pripar <- data$sigma_beta_comm_pripar
        
        # Store community parameters
        output$mu.alpha.comm <- res.extracted[["mu_alpha_comm"]]
        output$mu.alpha.comm.maxpost <- mu.alpha.comm.maxpost
        
        output$sigma.alpha.comm <- res.extracted[["sigma_alpha_comm"]]
        output$sigma.alpha.comm.maxpost <- sigma.alpha.comm.maxpost
        
        output$mu.beta.comm <- res.extracted[["mu_beta_comm"]]
        output$mu.beta.comm.maxpost <- mu.beta.comm.maxpost
        
        output$sigma.beta.comm <- res.extracted[["sigma_beta_comm"]]
        output$sigma.beta.comm.maxpost <- sigma.beta.comm.maxpost
        
        # Store posterior taxon-specific parameters
        output$alpha.taxa <- res.extracted[["alpha_taxa"]]
        output$alpha.taxa.maxpost <- alpha.taxa.maxpost
        
        output$beta.taxa  <- res.extracted[["beta_taxa"]]
        output$beta.taxa.maxpost <- beta.taxa.maxpost
        return(list(res,output))        
    }
}