# Formula creating ----
cat('Prepare data and recipe\n')
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

## Make stancode for VB model ----

### Prepare stancode recipe
params <- alist(
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
  
  save_pars = save_pars(all=TRUE),
  seed=seed
)

### Make stancode for VB model ----
cat('Make stancode for VB model\n')
rv_model_stancode <- do.call(make_stancode, params)
rv_model_standata <- do.call(make_standata, params)

cat('Run VB model\n')
rv_model_vb <- rstan::vb(
  rstan::stan_model(model_code=rv_model_stancode),
  data = rv_model_standata,
  seed=seed,
  iter=20000,
  tol_rel_obj=1e-3
)

### Run real model ----
cat('Extract estimation\n')
pars <- rv_model_vb@model_pars
pars <- pars[!pars %in% c('lp__')]
init <-
  sapply(
    pars,
    function(par){
      dim <- rv_model_vb@par_dims[[par]] 
      mean <-
        rstan::get_posterior_mean(rv_model_vb, par) |>
        c()
      if (length(dim))
        mean <- array(mean, dim=dim)
      mean
    },
    simplify = F
  )
# init <- rstan::summary(rv_model_vb)$summary[, 'mean'] |>
  # as.list()

cat('Run sampling\n')

rv_model_n <-
  do.call(
    brm, 
    modifyList(
      params, 
      alist(
        cores=1, chains=1,
        sample_file = '.cache/rv_model_n.csv', save_warmup=TRUE,
        iter=2000, 
        warmup=500, thin=1,
        init = list(init), init_r=0.75,
        refresh=200,
        control = list(adapt_delta=.8, max_treedepth=12)
      )
    )
  )
  # brm(
  #   
  #   cores=1, chains=1,
  #   sample_file = '.cache/rv_model_n3.csv', save_warmup=TRUE,
  #   
  #   iter=5000, warmup=3000, thin=1,
  #   refresh=200, init_r=0.75,
  #   # algorithm='meanfield',
  #   # tol_rel_obj = .00001,
  #   # iter=20000, 
  #   control = list(adapt_delta=.8, max_treedepth=12)
  # )

### Save model ----
cat('Save model\n')
saveRDS(list(data=rv, model=rv_model_n), file='results/rv_model_n.RDS')

