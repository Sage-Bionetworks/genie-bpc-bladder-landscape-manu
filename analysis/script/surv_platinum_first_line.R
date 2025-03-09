library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source) # also loads lots of packages.

dir_output <- here('data', 'survival', 'first_line_platinum')

surv_desc_fp <- here('data', 'survival')

read_wrap_clin <- function(p) {
  read_rds(file = here("data", 'cohort', p))
}

dft_pt <- read_wrap_clin("pt.rds")
dft_ca_ind <- read_wrap_clin("ca_ind.rds")
# we're taking the augmented version, which has TMB columns added.  The existing info is the same.
dft_cpt <- read_wrap_clin("cpt_aug.rds")
dft_med_onc <- readr::read_rds(here('data', 'cohort', "med_onc.rds"))
dft_med_onc %<>%
  augment_med_onc_imputed_ecog(., add_numeric = T)
dft_reg <- read_wrap_clin('reg.rds')

dft_lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dft_cases <- dft_lot %>%
  filter(line_therapy %in% 1) %>%
  filter(
    str_detect(regimen_drugs, 'Cisplatin|Carboplatin')
  )


# There was one case that contained both carboplatin and cisplatin.
# It's a zero-day regimen, so we're just going to exclude it.
dft_cases %<>%
  filter(!(str_detect(regimen_drugs, "Carbo") & 
             str_detect(regimen_drugs, "Cisplat"))) 

dft_cases %<>%
  mutate(
    reg_group = case_when(
      str_detect(regimen_drugs, "Carboplatin") ~ "Carboplatin-based",
      str_detect(regimen_drugs, "Cisplatin") ~ "Cisplatin-based",
      T ~ NA_character_
    ),
    reg_group = factor(reg_group)
  )

dft_cases <- left_join(
  dft_cases,
  select(dft_reg, record_id, ca_seq, regimen_number,
         dx_reg_start_int_yrs, tt_os_g_yrs, os_g_status),
  by = c("record_id", "ca_seq", "regimen_number"),
  relationship = "one-to-one"
)

dft_first_cpt <- get_first_cpt(
  ca_ind_dat = dft_ca_ind,
  cpt_dat = dft_cpt
)

dft_cases <- left_join(
  dft_cases,
  dft_first_cpt,
  by = c("record_id", "ca_seq")
)

dft_cases %<>%
  mutate(
    reg_start_cpt_yrs = dx_cpt_rep_yrs - dx_reg_start_int_yrs
  ) %>%
  select(-c(dx_cpt_rep_yrs)) # avoid confusion



# Add some additional variables which might be make for interesting confounders in a model.
dft_timf <- readr::read_rds(
  here('data', 'cohort', 'time_invariant_model_factors.rds')
)

dft_cases %<>%
  left_join(., dft_timf, by = c('record_id', 'ca_seq'))
dft_cases %<>% mutate(
    dob_reg_start_yrs = dob_ca_dx_yrs + dx_reg_start_int_yrs,
    dob_reg_start_days = dob_reg_start_yrs * 365.25, # just need for the next call
  )

dft_latest_med_onc <- get_med_onc_by_timing(
  dat_med_onc = dft_med_onc,
  dat_cutoff = dft_cases,
  var_cutoff = "dob_reg_start_days",
  remove_missing = T
)

dft_cases %<>%
  left_join(
    .,
    select(dft_latest_med_onc, record_id, md_ecog_imp_num),
    by = "record_id"
  )
  

# Get a list of regimens that were in subjects in the case list, and also
#   before the regimen that defines their case.
# This will help us to define whether people had previous therapies,
#   previous platinum therapies, etc.

dft_reg_before_case <- dft_cases %>%
  select(record_id, ca_seq, case_reg_num = regimen_number) %>%
  distinct(.) %>% # certainly redundant.
  left_join(
    .,
    dft_reg,
    by = c('record_id', 'ca_seq')
  ) %>%
  filter(regimen_number < case_reg_num)

dft_reg_before_case %<>% 
  group_by(record_id, ca_seq) %>%
  summarize(
    num_prev_ther = n(),
    num_prev_carbo = sum(str_detect(regimen_drugs, "Carboplatin")),
    num_prev_cis = sum(str_detect(regimen_drugs, "Cisplatin")),
    # makes the reasonable assumption that people don't take carboplatin
    #   and cisplatin in the same regimen.
    num_prev_plat = num_prev_carbo + num_prev_cis,
    num_prev_nonplat = num_prev_ther - num_prev_plat,
    .groups = "drop"
  )

dft_cases %<>%
  left_join(., dft_reg_before_case, by = c("record_id", "ca_seq")) %>%
  tidyr::replace_na(
    replace = list(
      num_prev_ther = 0,
      num_prev_carbo = 0,
      num_prev_cis = 0,
      num_prev_plat = 0,
      num_prev_nonplat = 0
    )
  ) %>%
  mutate(
    bin_prev_ther = num_prev_ther >= 1,
    bin_prev_plat = num_prev_plat >= 1,
    bin_prev_nonplat = num_prev_nonplat >= 1
  )

# We have two levels of regimen_drugs included, choosing cisplatin as a reference
#   due only to it being more common.
dft_cases %<>%
  mutate(
    carboplatin = if_else(
      str_detect(reg_group, "Carboplatin-based"),
      T,
      F
    )
  )

# Save this to be continued in the model script.
readr::write_rds(
  x = dft_cases,
  file = here(dir_output, 'data_first_line_plat.rds')
)
  

      
    
    
  
    






# Create the basic KM plots around this:

dft_cases %<>% remove_trunc_gte_event(
  trunc_var = 'reg_start_cpt_yrs',
  event_var = 'tt_os_g_yrs'
)

surv_obj_cases <- with(
  dft_cases,
  Surv(
    time = reg_start_cpt_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)




gg_first_line_platinum_manu <- plot_one_survfit(
  dat = dft_cases,
  surv_form = surv_obj_cases ~ reg_group,
  plot_title = "OS, from start of first line platinum-based chemotherapy",
  x_exp = 0.1,
  x_breaks = seq(0, 100, by = 0.5)
) + 
  add_quantile(y_value = 0.5, linetype = 81, alpha = 0.75) + 
  coord_cartesian(xlim = c(0,5), ylim = c(0, 1.01))

gg_first_line_platinum_manu <- gg_first_line_platinum_manu + 
  theme(
    plot.title.position = 'panel'
  )

readr::write_rds(
  gg_first_line_platinum_manu,
  here(dir_output, 'gg_first_line_platinum.rds')
)

ggsave(
  gg_first_line_platinum_manu,
  height = 5, width = 8,
  filename = here(
    'output', 'aacr_ss24', 'img',
    '03_first_line_platinum.pdf'
  )
)







