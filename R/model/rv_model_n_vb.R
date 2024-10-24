# Formula creating ----
cat('Prepare data and recipe\n')
source('R/model/rv_formulae.R')
sigma_fml <- sigma ~ zmuvol #s(zsdvol, k=3)

# Modeling ----

## Model specs ----
rv_fml_n <- 
  surv_fml + 
  bf(rv_fml,
     sigma_fml,
     family=gaussian(link='identity')) +
  set_rescor(FALSE)

## Prepare stancode recipe ----
params <- alist(
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
    
    prior(normal(0, 1), class='Intercept', dpar='sigma', resp='vol'),
    prior(normal(0, 1), class='b', dpar='sigma', resp='vol')
    # prior(exponential(1), class='sds', dpar='sigma', resp='vol')
  ),
  data=rv,
  save_pars = save_pars(all=TRUE),
  seed=seed,
  init_r=0.75
)


## Make dummy code ----
rv_stancode_n <- 
  do.call(
    make_stancode,
    modifyList(
      params,
      alist(stanvars = c(sgtbrms::expose_sgt_stanvar()))
    )
  )

params <- 
  modifyList(
    params, 
    alist(
      stanvars = c(sgtbrms::expose_sgt_stanvar(), 
                   rv_stanvar(rv_stancode_n, 'gaussian', n_subject_time))
    )
  )

### Make stancode for VB model ----
cat('Make stancode for VB model\n')
rv_model_stancode <- do.call(make_stancode, params)
rv_model_standata <- do.call(make_standata, params)

cat('Run VB model\n')
obj <-  rstan::stan_model(model_code=rv_model_stancode |> gsub('std_normal_lpdf(z_1[1])','0',x=_,fixed=T))
rv_model_vb <- rstan::vb(
  obj,
  data = rv_model_standata,
  seed=seed,
  iter=20000,
  init_r=0,
  grad_samples = 20,
  adapt_engaged = FALSE,
  eta=.1
)

### Run real model ----
dcat('Extract estimation\n')
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
        sample_file = '.cache/rv_model_n_vb.csv', save_warmup=TRUE,
        iter=2000, 
        warmup=500, thin=1,
        # init = list(init),
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
saveRDS(list(data=rv, model=rv_model_n), file='results/rv_model_n_vb.RDS')

