# Formula creating ----
source('R/model/rv_formulae.R')
sigma_fml <- sigma ~ zmuvol #s(zsdvol, k=3)

# Modeling ----

## Model specs ----
rv_fml_n <- 
  bf(rv_fml,
     sigma_fml,
     family=gaussian(link='identity')) +
  surv_fml + 
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
      
      # prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(student_t(6, 0, 1), class='b', resp='vol'),
      prior(exponential(1), class='sd', resp='vol'),
      
      prior(normal(0, 1), class='Intercept', dpar='sigma', resp='vol'),
      prior(normal(0, 1), class='b', dpar='sigma', resp='vol')
      # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
    ),
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
rv_model_n <-
  brm(
    rv_fml_n,
    prior = c(
      prior(logistic(0,1), class='Intercept', resp='evdeath'),
      # prior(sgtv(0.0, 3.0, 0.5, 10.0, 1.5), class='b', resp='evdeath'),
      prior(student_t(6, 0, 1), class='b', resp='evdeath'),
      prior(exponential(1), class='sds', resp='evdeath'),
      
      # prior(sgtv(0.0, 5.0, 0.5, 10.0, 1.5), class='b', resp='vol'),
      prior(student_t(6, 0, 1), class='b', resp='vol'),
      prior(exponential(1), class='sd', resp='vol'),
      
      prior(normal(0, 1), class='Intercept', dpar='sigma', resp='vol'),
      prior(normal(0, 1), class='b', dpar='sigma', resp='vol')
      # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
    ),
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 rv_stanvar(rv_stancode_n, 'gaussian', n_subject_time)),
    
    algorithm='meanfield',
    save_pars = save_pars(all=TRUE),
    cores=1, chains=1,
    # sample_file = '.cache/rv_model_n3.csv', save_warmup=TRUE,
    # iter=4000, warmup=2000, thin=1, 
    # refresh=200, init_r=0.75,
    tol_rel_obj = .001,
    seed=seed, iter=20000, 
    control = list(adapt_delta=.8, max_treedepth=12)
  )

### Save model ----
saveRDS(list(data=rv, model=rv_model_n), file='results/rv_model_n3.RDS')