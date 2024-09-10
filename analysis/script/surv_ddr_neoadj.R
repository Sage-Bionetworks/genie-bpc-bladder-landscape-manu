# In an August 2024 meeting we discussed whether we could add to the analysis
#   of DDR mutations by looking at the difference between clinical and path
#   staging in those individuals.


library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_output <- here('data', 'survival', 'ddr_neoadj')

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
# dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_alt <- readr::read_rds(here('data', 'genomic','alterations.rds'))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))
dft_tnm <- readr::read_rds(here('data', 'cohort', 'tnm_path_clin_details.rds'))


dft_neoadj <- dft_ca_ind %>% 
  select(record_id, ca_seq, ca_tx_pre_path_stage, stage_dx_iv) %>%
  left_join(
    .,
    select(dft_tnm, record_id, ca_seq, clin_group_clust, path_group_clust),
    by = c('record_id', 'ca_seq')
  )

# Talked to MSK curating team about this a bit - there is no reason people who
#.  are stage IV at dx should have ca_tx_pre_path_stage filled out.
#.  Should be all NA or Unknown.   We will just limit to people not stage IV.

dft_neoadj %<>%
  filter(!(stage_dx_iv %in% "Stage IV")) %>%
  filter(ca_tx_pre_path_stage %in% "Yes")

# We'll merge in the DDR stuff, write this out, then filter down to those that
#   have clinical and pathological staging.  This allows us to more easily look
#   at who is missing (due to missing one of those).


custom_ddr_list <- c(
  'ERCC2', 'ERCC5', 
  'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
  'ATR', 'FANCC'
)

# Change: taking only the first CPT test for DDR alterations.  It makes no sense
#   to get people who had one pop up later on for this analysis.
# This took the count from 128 to 123 (for alterations).
dft_onco_ddr <- dft_cpt %>%
  group_by(record_id, ca_seq) %>%
  arrange(cpt_number) %>%
  slice(1) %>%
  ungroup(.) %>%
  select(sample_id = cpt_genie_sample_id) %>%
  left_join(
    ., dft_alt, by = 'sample_id'
  ) %>%
  filter(hugo %in% custom_ddr_list) %>%
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

# Add in the relevant stuff from CPT data:
dft_onco_ddr <- dft_cpt %>%
  select(
    record_id, ca_seq,
    dx_path_proc_cpt_yrs, # interval from dx to pathology procedure of this CPT. 
    dx_cpt_rep_yrs, # interval from dx to report date of this CPT.
    cpt_genie_sample_id
  ) %>%
  left_join(
    dft_onco_ddr,
    .,
    by = c(sample_id = "cpt_genie_sample_id")
  )

dft_onco_ddr_sum <- dft_onco_ddr %>%
  group_by(record_id, ca_seq) %>%
  summarize(onco_ddr = T) 

dft_neoadj %<>%
  left_join(
    .,
    dft_onco_ddr_sum,
    by = c('record_id', 'ca_seq')
  ) %>%
  replace_na(list(onco_ddr = F))


# Save a version so you can see who's missing
readr::write_rds(
  dft_neoadj,
  here(dir_output, 'ddr_neoadj_unfiltered.rds')
)

dft_neoadj %<>%
  filter(
    !is.na(clin_group_clust) & !is.na(path_group_clust)
  )

lev_gs_change <- c(
  "Stage lower at path",
  "No change in stage",
  "Stage higher at path"
)
  

dft_neoadj %<>%
  mutate(
    group_stage_change_num = case_when(
      path_group_clust == clin_group_clust ~ 0L,
      path_group_clust > clin_group_clust ~ 1L,
      path_group_clust < clin_group_clust ~ -1L,
      T ~ NA_integer_ # should never happen
    ),
    group_stage_change_f = case_when(
      path_group_clust == clin_group_clust ~ lev_gs_change[2],
      path_group_clust > clin_group_clust ~ lev_gs_change[3],
      path_group_clust < clin_group_clust ~ lev_gs_change[1]
    ),
    group_stage_change_f = factor(group_stage_change_f, levels = lev_gs_change),
    onco_ddr_disp = factor(
      case_when(
        onco_ddr ~ "DDR alt.",
        T ~ "No oncogenic DDR"
      )
    )
  )

readr::write_rds(
  dft_neoadj,
  here(dir_output, 'ddr_neoadj.rds')
)

gg_neoadj_ddr_mosaic <- ggplot(
  data = dft_neoadj
) + 
  geom_mosaic(
    aes(x = product(onco_ddr_disp), fill = group_stage_change_f)
  ) + 
  theme_mosaic() + 
  scale_fill_viridis_d(
    name = NULL,
    option = "magma", begin = 0, end = 0.5
  ) +
  theme(
    axis.title = element_blank()
  )

readr::write_rds(
  gg_neoadj_ddr_mosaic,
  here(dir_output, 'gg_neoadj_ddr_mosaic.rds')
)

gg_clin_path_group_stage <- dft_neoadj %>%
  select(record_id, onco_ddr_disp,
         clin_group_clust, path_group_clust) %>%
  pivot_longer(
    cols = contains("group_clust")
  ) %>%
  mutate(name = str_sub(name, 1, 4)) %>%
  ggplot(
    data = .,
  ) + 
  geom_mosaic(
    aes(x = product(name), fill = value)
  ) + 
  theme_mosaic() + 
  scale_fill_viridis_d(
    name = NULL,
    option = "turbo", begin = 0, end = 1,
    guide = guide_legend(reverse = T)
  ) +
  theme(
    axis.title = element_blank(),
    legend.position = "left",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) + 
  facet_wrap(vars(onco_ddr_disp), nrow = 1) 

readr::write_rds(
  gg_clin_path_group_stage,
  here(dir_output, 'gg_clin_path_group_stage.rds')
)


      
      

