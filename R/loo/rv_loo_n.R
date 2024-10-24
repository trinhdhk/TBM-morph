# LOO for model ----
library(brms, quietly = TRUE)
# rstan::rstan_options(threads_per_chain=10)
options(future.globals.maxSize = 25*1024^3)
future::plan('multicore', workers=4, gc=TRUE)
# future::plan('sequential')
set.seed(1208)

K <- 10
model <- readRDS('results/rv_model_n.RDS')
if (!file.exists('.cache/foldid.RDS')){
  foldid <- loo::kfold_split_grouped(K=K,  x=model$data$obj)
  saveRDS(foldid, '.cache/foldid.RDS')
} else {
  foldid <- readRDS('.cache/foldid.RDS')
}

# Create kfold loglik matrix
fulldata <- model$data |> transform(foldid = foldid)
compdata <- fulldata |> subset(!is.na(vol))
# loglik <- matrix(NA_real_, nrow = 1000, ncol = nrow(compdata))

# Run the progress
source('R/model/rv_formulae.R')
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

progressr::with_progress(
  enable = TRUE,
  {
    pr <- progressr::progressor(K)
    fits <- future.apply::future_lapply(
      seq_len(K),
      function(k){
        traink <- subset(fulldata, foldid!=k)
        testk  <- subset(compdata, foldid==k)
        subject_time <- with(traink, paste(obj, timeInt)) 
        n_subject_time <- length(unique(subject_time))
        sigma_fml <- sigma ~ zmuvol #s(zsdvol, k=3)
        rstan::rstan_options(threads_per_chain=10)
        
        # Modeling ----
        
        ## Model specs ----
        rv_fml_n <- 
          surv_fml + 
          bf(rv_fml,
             sigma_fml,
             family=gaussian(link='identity')) +
          set_rescor(FALSE)
        
        
        ## Make dummy code ----
        rv_stancode_n <- 
          make_stancode(
            rv_fml_n,
            prior = c(
              prior(logistic(0,1), class='Intercept', resp='evdeath'),
              # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
              prior(student_t(6, 0, 1), class='b', resp='evdeath'),
              prior(exponential(1), class='sds', resp='evdeath'),
              prior(exponential(1), class='sd', resp='evdeath'),
              
              # prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
              prior(student_t(6, 0, 1), class='b', resp='vol'),
              prior(exponential(1), class='sd', resp='vol'),
              
              prior(std_normal(), class='Intercept', dpar='sigma', resp='vol'),
              prior(std_normal(), class='b', dpar='sigma', resp='vol')
              # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
            ),
            data=traink, 
            stanvars = c(sgtbrms::expose_sgt_stanvar()))
        
        modelk <- 
          brm(
            rv_fml_n,
            prior = c(
              prior(logistic(0,1), class='Intercept', resp='evdeath'),
              # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
              prior(student_t(6, 0, 1), class='b', resp='evdeath'),
              prior(exponential(1), class='sds', resp='evdeath'),
              prior(exponential(1), class='sd', resp='evdeath'),
              
              # prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
              prior(student_t(6, 0, 1), class='b', resp='vol'),
              prior(exponential(1), class='sd', resp='vol'),
              
              prior(std_normal(), class='Intercept', dpar='sigma', resp='vol'),
              prior(std_normal(), class='b', dpar='sigma', resp='vol')
              # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
            ),
            data=traink, 
            stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                         rv_stanvar(rv_stancode_n, 'gaussian', n_subject_time)),
            
            cores=1, chains=1,
            seed=seed,
            iter=5000, warmup=3000, thin=1,
            refresh=200, init_r=0.75,
            control = list(adapt_delta=.78, max_treedepth=12)
          )
        
        # loglik[,compdata$foldid==k] <- 
          # log_lik(modelk, 
          #         newdata=testk, 
          #         ndraws=1000, 
          #         cores=4,
          #         resp='vol')
        ll <- log_lik(modelk, 
                      newdata=testk, 
                      ndraws=1000, 
                      cores=4,
                      resp='vol',
                      allow_new_levels=TRUE)
        rm(modelk)
        gc()
        pr()
        c(id=which(compdata$foldid==k), loglik = ll)
      }, future.seed=TRUE)
  }
) 

saveRDS(fits, 'results/loglik_n.RDS')

saveRDS(loo_n, '.cache/icv_loo_n.RDS')