library(tidytable)
library(survival)
library(discSurv)

if (!exists('approx')) approx <- 10
cat(approx,  '\n')
clin <-  readRDS('data/imported/clinical.RDS') |> 
  rename(subject=usubjid) |>
  mutate(tt_death = tt_death/7)

interval_limits <-
  c(clin[ev_death==1,tt_death]) |> quantile(probs = seq(0.1, 1, by=1/approx)) |>
  c(max(clin$tt_death)+1)
  
interval_tbl <- tidytable(
  timeInt = seq_along(interval_limits),
  tstart = lag(interval_limits, default=0),
  tstop = interval_limits
)
  
clin_disc <- 
  clin |>
  as.data.frame() |>
  contToDisc(timeColumn = 'tt_death',
             intervalLimits=interval_limits)  

clin_disc_long <- 
  dataLong(clin_disc, timeColumn = 'timeDisc', eventColumn = 'ev_death')


vols <- merge(readRDS('data/imported/volumetry.RDS'), select(clin, c(subject, bldate))) |>
  as_tidytable() |>
  mutate(mri_wk =  as.integer(difftime(date, bldate, units='d')) + 751) |>
  filter(mri_wk > -30) |>
  mutate(mri_wk = pmax(mri_wk,0)/7) |>
  filter(mri_wk < 10) |>
  select(-date, -bldate) 

vols_session <- vols |>
  filter(number==4) |>
  select(-number, -name) |>
  as.data.frame() |>
  contToDisc(timeColumn = 'mri_wk',
             intervalLimits=interval_limits) |>
  as_tidytable() |>
  rename(timeInt = timeDisc)

clin_vols <- 
  left_join(clin_disc_long, vols_session, by=c('subject', 'timeInt')) |>
  dplyr::group_by(obj, timeInt) |>
  dplyr::group_modify(
    ~ cbind(.x, name = unique(vols$name))
  ) |>
  ungroup() |>
  dplyr::group_by(obj) |>
  dplyr::group_modify(
    function(x,y){
      if (all(is.na(x$vol))) return(filter(x, FALSE))
      x
    }
  ) |>
  select(-vol) |>
  left_join(select(vols, subject, name, number, mri_wk, vol),
            by = c('subject', 'name', 'mri_wk')) |>
  rename(region=name) |>
  left_join(interval_tbl) |>
  mutate(mri_wk = ifelse(is.na(mri_wk), tstart, mri_wk))

saveRDS(clin_vols, 'data/imported/clin_vols.RDS' )
write.csv(clin_vols, 'data/imported/clin_vols.csv', row.names = FALSE)

