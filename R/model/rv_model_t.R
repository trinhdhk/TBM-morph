# Formula creating ----
source('R/model/rv_formulae.R')
sigma_fml <- sigma ~ zmuvol #s(zsdvol, k=3)
# Sys.setenv(STAN_NUM_THREADS=5)

# Modeling ----

## Model specs ----
rv_fml_t <- 
  surv_fml + 
  bf(rv_fml,
     sigma_fml,
     family=student(link='identity')) +
  set_rescor(FALSE)


## Make dummy code ----
rv_stancode_t <- 
  make_stancode(
    rv_fml_t,
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
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar()))

### Run real model ----
rv_model_t <-
  brm(
    rv_fml_t,
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
    data=rv, 
    stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                 rv_stanvar(rv_stancode_t, 't', n_subject_time)),
    
    save_pars = save_pars(all=TRUE),
    cores=1, chains=1,
    sample_file = '.cache/rv_model_t.csv', save_warmup=TRUE,
    seed=seed,
    iter=6000, warmup=3000, thin=1,
    refresh=200, init_r=0.75,
    # algorithm='meanfield',
    # tol_rel_obj = .00001,
    # iter=20000, 
    control = list(adapt_delta=.8, max_treedepth=11)
  )

### Save model ----
saveRDS(list(data=rv, model=rv_model_t), file='results/rv_model_t.RDS')

