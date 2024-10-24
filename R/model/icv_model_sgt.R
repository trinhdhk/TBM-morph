# Formula creating ----
source('R/model/icv_formulae.R')

# Modeling ----

## Model specs ----
icv_fml <- 
  bf(icv_fml,
     family=sgtbrms::sgt(link='identity')) +
  surv_fml + 
  set_rescor(FALSE)


## Make dummy code ----
icv_stancode <- 
  make_stancode(
    icv_fml,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol'),
      prior(beta(0.5, 0.5), class='lambdap1half', resp='vol')
    ),
    data=icv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
icv_model <-
  brm(
    icv_fml,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      
      prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(exponential(1), class='sigma', resp='vol'),
      prior(beta(0.5, 0.5), class='lambdap1half', resp='vol')
    ),
    data=icv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 icv.stanvar(icv_stancode, 'sgt')),
    
    save_pars = save_pars(all=TRUE),
    cores=2, chains=2,
    sample_file = '.cache/icv_model_sgt.csv', save_warmup=TRUE,
    iter=6000, warmup=3000, thin=1, seed=seed, refresh=200, init_r=0.75,
    control = list(adapt_delta=.88, max_treedepth=12)
  )

### Save model ----
saveRDS(list(data=icv, model=icv_model), file='results/icv_model_sgt.RDS')
