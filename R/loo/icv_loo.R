# LOO for model ----
library(brms, quietly = TRUE)
rstan::rstan_options(threads_per_chain=10)
options(future.globals.maxSize = 10*1024^3)
future::plan('multicore', workers=10, gc=TRUE)

## Load model ----
model_n <- readRDS('results/icv_model_n.RDS')
model_t <- readRDS('results/icv_model_t.RDS')
model_gt <- readRDS('results/icv_model_symgt.RDS')

compdata <- model_n$data |>
  subset(!is.na(vol))

foldid <-
  loo::kfold_split_grouped(x=model_n$data$obj)

foldid[is.na(model_n$data$vol)] <- 11

model_n$name <- 'n'
model_t$name <- 't'
model_gt$name <- 'symgt'

## Run LOO ----
myloo <- purrr::partial(kfold,
  resp='vol', folds='grouped', group='obj',
  # Ksub = 1:10,
  # joint='group',
  newdata = model_n$data, recompile=FALSE, chains=1
)

# loos <- future.apply::future_lapply(
#   list(model_n, model_t, model_gt),
#   function(model){
#     rstan::rstan_options(threads_per_chain=10)
#     thisloo = myloo(model$model)
#     saveRDS(thisloo, paste0('.cache/icv_loo_',model$name,'.RDS'))
#     thisloo
#   },
#   future.seed = TRUE
# )
# 
# names(loos) <- c('n', 't', 'symgt')
# 
# saveRDS(loos, 'results/icv_loo.RDS')

cat('Normal \n')
loo_n <- myloo(x = model_n$model)
saveRDS(loo_n, '.cache/icv_loo_n.RDS')

cat('T \n')
loo_t <- myloo(x = model_t$model)
saveRDS(loo_t, '.cache/icv_loo_t.RDS')

cat('GT \n')
loo_symgt <- myloo(x = model_gt$model)
saveRDS(loo_symgt, '.cache/icv_loo_symgt.RDS')


saveRDS(list(loo_n = loo_n,
             loo_t = loo_t,
             loo_symgt = loo_symgt
             ), file = 'results/icv_loo.RDS')

