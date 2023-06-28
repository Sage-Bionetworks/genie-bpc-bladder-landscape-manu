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


# generally we'll just be interested in regimens associated with the index ca
dft_reg <- inner_join(
  dft_reg,
  select(dft_ca_ind, record_id, ca_seq),
  by = c("record_id", "ca_seq")
)

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




# Write datasets to derived location:
write_wrap <- function(obj, file_name) {
  readr::write_csv(x = obj, file = here('data', 'cohort', file_name))
}

dft_img <- read_wrap("imaging_level_dataset.csv")
dft_med_onc <- read_wrap("med_onc_note_level_dataset.csv")
dft_path <- read_wrap("pathology_report_level_dataset.csv")
dft_reg <- read_wrap("regimen_cancer_level_dataset.csv")
dft_cpt <- read_wrap("cancer_panel_test_level_dataset.csv")

write_wrap(dft_pt, file_name = "pt.csv")
write_wrap(dft_ca_ind, file_name = "ca_ind.csv")
write_wrap(dft_img, file_name = "img.csv")
write_wrap(dft_med_onc, file_name = "med_onc.csv")
write_wrap(dft_path, file_name = "path.csv")
write_wrap(dft_reg, file_name = "reg.csv")
write_wrap(dft_cpt, file_name = "cpt.csv")


    
    
    
