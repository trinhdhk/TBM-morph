library(tidytable)
library(brms)

seed = 1208

# Helper function ----

rv_stanvar = function(model_code, family, n_subject_time){
  pattern <- "// likelihood including constants(.*?)// priors including constants"
  likelihood <- regmatches(model_code, regexec(pattern, model_code))[[1]][2]
  # likelihood <- gsub('if (!prior_only) {', '', likelihood, fixed=TRUE)
  likelihood <- substr(likelihood, 1, nchar(likelihood)-4)
  
  # adding subject-time list to stan
  subject_time <-
    brms::stanvar(
      x=n_subject_time,
      name='n_subject_time',
      block='data'
    )
  
  surv_lpmf <- 
    brms::stanvar(
      name='surv_lpmf',
      scode='
        real surv_lpmf(int[] y, matrix Xc_evdeath, vector mu_evdeath, vector b_evdeath, vector mu_vol, real b_vol, int n_subject_time){
          real lpmf = 0;
          vector[n_subject_time] lp_n;
          array[n_subject_time] int Y_evdeath;
          for (subject_time in 1:n_subject_time){
            int n_region = 72;
            int i_start = (subject_time - 1)*72 + 1;
            int i_stop = subject_time*72;
            Y_evdeath[subject_time] = y[i_start];
            if (y[i_start] != mean( y[i_start:i_stop])) reject("BUG found!");
            lp_n[subject_time] = sum(mu_evdeath[i_start:i_stop] + Xc_evdeath[i_start:i_stop,:] * b_evdeath)/72;
            lp_n[subject_time] += b_vol * sum(mu_vol[i_start:i_stop]) / 72;
          }
          lpmf += bernoulli_logit_lpmf(Y_evdeath | lp_n);
          return lpmf;
        } 
      ',
      block='functions'
    )
  old_lpdf <- 
    switch(
      family,
      gaussian = 'target += normal_lpdf(Yl_vol | mu_vol, sigma_vol)',
      skew_normal = 'target += skew_normal_lpdf(Yl_vol | mu_vol, omega_vol, alpha_vol)',
      t = 'target += student_t_lpdf(Yl_vol | nu_vol, mu_vol, sigma_vol)',
      st = 'for (n in 1:N_vol) {\n      target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, qmhalf_vol[n]);\n    }',
      sgt = 'for (n in 1:N_vol) {\n      target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, p_vol[n], q_vol[n]);\n    }'
    )
  new_lpdf <-
    switch(
      family,
      gaussian = 'target += normal_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      skew_normal = '
      target += skew_normal_lpdf(Yl_vol[Jobs_vol] | mu_vol[Jobs_vol], omega_vol[Jobs_vol], alpha_vol);
      mu_vol += omega_vol * delta_vol * sqrt(2/pi());',
      t = 'target += student_t_lpdf(Yl_vol[Jobs_vol] | nu_vol[Jobs_vol], mu_vol[Jobs_vol], sigma_vol[Jobs_vol]);',
      st = 'for (n in Jobs_vol) target += constrained_skew_t_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, qmhalf_vol[n]);',
      sgt = 'for (n in Jobs_vol) target += sgt_lpdf(Yl_vol[n] | mu_vol[n], sigma_vol[n], lambdap1half_vol, p_vol[n], q_vol[n]);'
    )
  surv.code <- '
  /*{
    vector[N_1] r2_1_evdeath_1;
    r2_1_evdeath_1[1] = 0;
    r2_1_evdeath_1[2:] = r_1_evdeath_1[2:];
    
    mu_evdeath += (bsp_evdeath[1] + r_1_evdeath_1[J_1_evdeath]) .* (-Yl_vol) +  (bsp_evdeath[1] * exp(r2_1_evdeath_1[J_1_evdeath])) .* mu_vol;
  }*/
  //mu_evdeath += (bsp_evdeath[1] + r_1_evdeath_1[J_1_evdeath]) .* (-Yl_vol) +  (bsp_evdeath[1] * exp(r_1_evdeath_1[J_1_evdeath]-r_1_evdeath_1[J_1_evdeath[1]])) .* mu_vol;
  //mu_evdeath += bsp_evdeath[1] .* (mu_vol[Jmi_vol]-Yl_vol[Jmi_vol]);
  mu_evdeath +=  bsp_evdeath[1] .* (-Yl_vol);
  
  //target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);
  target += surv_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath, mu_vol, bsp_evdeath[1], n_subject_time);
  target += std_normal_lpdf(Ymi_vol);'
  likelihood <- gsub('target += bernoulli_logit_glm_lpmf(Y_evdeath | Xc_evdeath, mu_evdeath, b_evdeath);', '', likelihood, fixed=TRUE)
  likelihood <- gsub(old_lpdf, paste(new_lpdf, surv.code, sep='\n'), likelihood, fixed=TRUE)
  likelihood <- gsub('Yl_vol[Jmi_vol] = Ymi_vol', 'Yl_vol[Jmi_vol] = mu_vol[Jmi_vol]', likelihood, fixed=TRUE)
  
  stopifnot(length(new_lpdf) > 0)
  tdata_code_def <- 'int Jobs_vol[N_vol - Nmi_vol];
  matrix[N_vol - Nmi_vol,K_vol] Xobs_vol;' |>
    brms::stanvar(name='obs_tdata_def', scode=_, block='tdata', position='start')
  tdata_code_end <- '
  {
    int k = 1;
    for (n in 1:N_vol) {
      int mi = 0;
      for (m in 1:Nmi_vol) {
       
        if (Jmi_vol[m] == n) {
          mi = 1;
          break;
        }
      }
      if (mi == 0){
        Jobs_vol[k] = n;
        k += 1;
      }
    }
    Xobs_vol = X_vol[Jobs_vol,:];
  }' |>  brms::stanvar(name='obs_tdata', scode=_, block='tdata', position='end')
  ll_1 <- '/*' |>  brms::stanvar(name='cancel_mid', scode=_, block='model', position='start')
  ll_2 <- '*/' |>  brms::stanvar(name='cancel_out', scode=_, block='likelihood', position='end')
  ll <- likelihood |> brms::stanvar(name='new', scode=_, block='likelihood', position='end')
  c(tdata_code_def, tdata_code_end, ll_1, ll_2, ll, subject_time, surv_lpmf)
}

# Data -----
approx <- 5
source('R/02_prepare_data.R')
rv <- 
  clin_vols |>
  # readRDS('data/imported/clin_vols.RDS') |>
  filter(!grepl('(ventricular)|(entricl)|(Vent)', region)) |>
  mutate(vol = sum(vol), .by=c(timeInt, obj, region)) |>
  select(-number) |>
  unique() 

mu_cd4 <- mean(log2(rv[hiv==TRUE, cd4]), na.rm=TRUE)  
sd_cd4 <- sd(log2(rv[hiv==TRUE, cd4]), na.rm=TRUE)
rv <- mutate(rv,
              weight = scale(weight) |> as.vector(),
              csflym = scale(log10(csflym)) |> as.vector(),
              cd4 = ifelse(hiv, (log2(cd4) - mu_cd4)/sd_cd4, 0)) |>
  mutate(cd4 = replace_na(cd4, mu_cd4)) |>
  select(-ev_death) |>
  rename(ev_death=y) |>
  mutate(strata = paste(arm, as.numeric(hiv), mrc, genotype, sep='-')) |>
  mutate(
    severe=as.numeric(mrc)>1,
    # mrc = as.numeric(mrc)>1, #factor(mrc, ordered=FALSE),
    sd_wk = sd(mri_wk, na.rm=TRUE),
    mu_wk = mean(mri_wk, na.rm=TRUE),
    mri_wk = mri_wk/sd_wk,
    arm = as.numeric(arm)-1) |>
  mutate(across(c(tt_death, tstart, tstop), ~ (.x)/sd_wk)) |>
  mutate(
    muvol = mean(log(vol), na.rm=T),
    sdvol = sd(log(vol), na.rm=T),
    vol = scale(log(vol)) |> as.vector(),
    .by = region
  ) |>
  mutate(
    zmuvol = scale(muvol),
    zsdvol = scale(sdvol)
  ) |>
  arrange(obj, timeInt)

# get the subject time list
subject_time <- with(rv, paste(obj, timeInt)) 
n_subject_time <- length(unique(subject_time))

# Formulae ----
surv_fml <-
  bf(ev_death ~ mi(vol) +
       s(tstop, k = 5, by=hiv) + 
       csflym + cd4 + arm + 
       arm:hiv + csflym:hiv + 
       mo(genotype) + mo(genotype):hiv, 
     family=bernoulli(link='logit')) 

rv_fml <- 
  bf(
    vol | mi() ~ 0+Intercept +
      mri_wk +
      mo(genotype) + mo(mrc) + hiv + 
      mo(mrc):mri_wk + mo(genotype):mri_wk + arm:mri_wk + hiv:mri_wk + 
      arm:mo(mrc):mri_wk + arm:mo(genotype):mri_wk + arm:hiv:mri_wk +
      I(mrc>1):hiv:mri_wk +  arm:I(mrc>1):hiv:mri_wk +
      (1|obj) +
      (1 + mri_wk +
         mo(genotype) + mo(mrc) + hiv + 
         mo(mrc):mri_wk + mo(genotype):mri_wk + arm:mri_wk + hiv:mri_wk + 
         arm:mo(mrc):mri_wk + arm:mo(genotype):mri_wk + arm:hiv:mri_wk + 
         I(mrc>1):hiv:mri_wk +  arm:I(mrc>1):hiv:mri_wk | gr(region, cor=FALSE)
      )
  )






