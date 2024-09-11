library(brms)
library(ggplot2)

icv_pred = conditional_effects(
  icv_model_n$model, 
  'arm',
  conditions = expand.grid(genotype=c('CC', 'CT', 'TT'),
                           hiv = c(F,T), mrc=1:3))$vol.vol_arm
icv_pred <- dplyr::select(icv_pred,-dplyr::ends_with('__'))
icv_pred$.id <- 1:nrow(icv_pred)

icv_pred <- 
  rbind(
    icv_pred |> dplyr::mutate(mri_wk = 0, t = 0),
    icv_pred |> dplyr::mutate(mri_wk = (60/7)/icv_model_n$data$sd_wk[[1]], t = 1)
  ) |> tidybayes::add_epred_draws(icv_model_n$model, ndraws=500, resp='vol') 

icv_pred <- 
  icv_pred |>
  #dplyr::filter(.category=='vol') |>
  dplyr::mutate(.epred = .epred * icv_model_n$data$sdvol[[1]] + icv_model_n$data$muvol[[1]]) |>
  dplyr::group_by(t) |>
  tidyr::pivot_wider(id_cols=c(.id, .draw, mrc, genotype, arm, hiv), names_from='t', values_from='.epred', names_prefix='t.') |>
  dplyr::mutate(volChange=exp(t.1 - t.0))
ggplot(icv_pred |> filter(arm%in%c(0,1)) |> mutate(arm=factor(arm), hiv=factor(hiv, labels=c("HIV-", "HIV+"))),
       aes(y=volChange, x=genotype, color=arm)) + 
  ggdist::stat_pointinterval(position=position_dodge(.2)) + facet_grid(hiv~mrc)
