library(fs)
library(purrr)
library(here)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

fs::dir_create('data', 'cohort')

read_wrap <- function(p) {
  read_csv(file = here("data-raw", p), show_col_types = F)
}

dft_pt <- read_wrap("patient_level_dataset.csv")
dft_ca_ind <- read_wrap("cancer_level_dataset_index.csv")
dft_img <- read_wrap("imaging_level_dataset.csv")
dft_med_onc <- read_wrap("med_onc_note_level_dataset.csv")
dft_path <- read_wrap("pathology_report_level_dataset.csv")
dft_reg <- read_wrap("regimen_cancer_level_dataset.csv")
dft_cpt <- read_wrap("cancer_panel_test_level_dataset.csv")
dft_rad <- read_wrap("ca_radtx_dataset.csv")

# A few sanity checks on the data:
if ((dft_pt$record_id %>% duplicated %>% any)) {
  stop("Duplicated records in patient level dataset.")
}
# At the time we started working there was only one index cancer per
#  record id.  Double checking this as we go:
if ((dft_ca_ind %>% 
     count(record_id, ca_seq, sort = T) %>%
     pull(n) %>% is_greater_than(1) %>% any)) {
  stop("Some patients have >1 index cancer - adjust as needed.")
}

# Additional filtering can be done here.  








dft_cohort_keys <- dft_ca_ind %>% select(record_id, ca_seq)
chk_keys_unique <- count(dft_cohort_keys, record_id, ca_seq) %>%
  pull(n) %>% 
  max %>%
  is_in(1)
if (!chk_keys_unique) {
  cli::cli_abort("Duplicate keys found in create_cohort_data.R")
}

key_filt_help <- function(dat) {
  inner_join(
    dft_cohort_keys,
    dat,
    by = c("record_id", "ca_seq"),
    multiple = "all" # default behavior in SQL and dplyr - just silences.
  )
}

n_row_cpt_old <- nrow(dft_cpt)
n_row_reg_old <- nrow(dft_reg)
n_row_rad_old <- nrow(dft_rad)

dft_cpt %<>% key_filt_help(.)
dft_reg %<>% key_filt_help(.)
dft_rad %<>% key_filt_help(.)

cli_alert_info(glue("{n_row_cpt_old-nrow(dft_cpt)} rows removed from dft_cpt for being related to non index cancers"))
cli_alert_info(glue("{n_row_reg_old-nrow(dft_reg)} rows removed from dft_reg for being related to non index cancers"))


# Create additional derived variables.
lev_st_simple <- c("Primary tumor", "Metastatic", "Other")
dft_cpt %<>%
  mutate(
    sample_type_simple_f = case_when(
      is.na(sample_type) ~ NA_character_,
      sample_type %in% "Local recurrence" ~ lev_st_simple[1],
      sample_type %in% "Lymph node metastasis" ~ lev_st_simple[2],
      sample_type %in% "Metastasis site unspecified" ~ lev_st_simple[2],
      sample_type %in% "Not applicable or hematologic malignancy" ~ lev_st_simple[3],
      sample_type %in% "Primary tumor" ~ lev_st_simple[1],
      T ~ NA_character_,
    ),
    sample_type_simple_f = factor(sample_type_simple_f, levels = lev_st_simple)
  )


# Add in a text description of ca_d_site (useless because too many are NOS but whatever)
dft_icd <- readxl::read_xlsx(
  here('data-raw', 'manual', 'icd_topography.xlsx')
)

# there are two cases which weren't found in the downloaded table and we need.  I pulled these from the data guide:i

dft_icd %<>%
  add_row(
    icdo3_code = "C67.3",
    description = "Anterior wall of bladder",
  ) %>%
  add_row(
    icdo3_code = "C67.5",
    description = "Bladder neck",
  ) %>%
  arrange(icdo3_code) %>%
  rename(ca_d_site_txt = description)

dft_ca_ind %<>%
  count(ca_d_site) %>% 
  left_join(
    .,
    dft_icd,
    by = c(ca_d_site = "icdo3_code")
  )






# This is needed so the panels match what we have panel data on:
cli::cli_alert_info(
  "Correcting any entries of 'UHN-OCA-v3' to 'UHN-OCA-V3'."
)
dft_cpt %<>%
  mutate(
    cpt_seq_assay_id = if_else(
      cpt_seq_assay_id %in% 'UHN-OCA-v3', 'UHN-OCA-V3',
      cpt_seq_assay_id
    )
  )










# Write datasets to derived location:
write_wrap <- function(obj, file_name) {
  readr::write_rds(
    x = obj, 
    file = here('data', 'cohort', paste0(file_name, ".rds"))
  )
}

write_wrap(dft_pt, file_name = "pt")
write_wrap(dft_ca_ind, file_name = "ca_ind")
write_wrap(dft_img, file_name = "img")
write_wrap(dft_med_onc, file_name = "med_onc")
write_wrap(dft_path, file_name = "path")
write_wrap(dft_reg, file_name = "reg")
write_wrap(dft_cpt, file_name = "cpt")
write_wrap(dft_rad, file_name = "rad")



    
    
    
