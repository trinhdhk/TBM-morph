library(tidytable)

# 26 TB ----
load('data/clinical/26TB_vad.Rdata')
randolist <- read.csv('data/clinical/26TB_RandlistFinal.csv')
base$RanNo <- as.integer(base$RANNO)
lta4h <- readxl::read_excel('data/clinical/26TB LTA4H genotyping_Mar2022.xlsx') 
lta4h <- lta4h[-1:-2,c(2, 5)]
colnames(lta4h) <- c('USUBJID', 'lta4hresult')
lta4h <- filter(lta4h, grepl('003-', USUBJID))
lta4h$USUBJID <- gsub('26TB-003-', '', lta4h$USUBJID)
tb26 <- inner_join(base, randolist[, c('arm', 'RanNo')], by='RanNo') |>
  filter(grepl('003-', USUBJID)) |>
  left_join(primary_endpoint |> select(USUBJID, TT_DEATH, EV_DEATH)) |>
  mutate(USUBJID=gsub('003-', '', USUBJID)) |>
  inner_join(lta4h) |>
  select(usubjid=USUBJID, arm,
         weight=WEIGHT,
         genotype = lta4hresult,
         bldate = BL_DATE,
         mrc = TBMGRADE,
         csflym = CSFLYMTOT,
         cd4 = CD4,
         na = SODIUM,
         ev_death = EV_DEATH,
         tt_death = TT_DEATH,
         ) |>
  mutate(usubjid=as.integer(usubjid), 
         genotype = ifelse(genotype=='TC', 'CT', genotype) ,
         tt_death=as.integer(tt_death), hiv=TRUE) |>
  filter(usubjid!=215)

# 27 TB ----
load('data/clinical/27allocated.Rdata')
tb27 <- bsl_ran |>
  filter(grepl('003-', usubjid)) |>
  select(usubjid, arm,
         weight=weight,
         genotype=lta4hresult,
         bldate=daterandom2,
         mrc =tbmgrade,
         csflym = csf_lym,
         na = sodium) |>
  inner_join(endpoints_ran2 |> 
               select(usubjid, ev_death = ev_death_12m,
                      tt_death = tt_death_12m)) |>
  mutate(usubjid = gsub('003-', '', usubjid)) |>
  mutate(usubjid=as.integer(usubjid), cd4=0,
         mrc = as.integer(mrc), hiv=FALSE) |>
  filter(usubjid != 2061)
clin_data <- rbind(tb26, tb27) |>
  mutate(genotype=factor(genotype, levels=c('CC','CT', 'TT'), ordered = TRUE))
write.csv(clin_data, 'data/imported/clinical.csv',
          row.names = FALSE)
saveRDS(clin_data, 'data/imported/clinical.RDS')
