make_obs_count_df <- function(med_onc_dat,
                              img_dat,
                              reg_dat,
                              ca_ind_dat) {
  
  dft_count_med_onc <- med_onc_dat %>%
    count(record_id) %>%
    rename(n_med_onc = n)
  
  img_type_levs <- c("n_scan_ct",
                     "n_scan_mri",
                     "n_scan_pet_or_pet_ct",
                     "n_scan_pet_mri",
                     "n_scan_bone",
                     "n_scan_other")
  
  dft_count_img <- img_dat %>%
    count(record_id, image_scan_type) %>%
    mutate(
      image_scan_type = case_when(
        image_scan_type %in% "CT" ~ img_type_levs[1],
        image_scan_type %in% "MRI" ~ img_type_levs[2],
        image_scan_type %in% "PET or PET-CT" ~ img_type_levs[3],
        image_scan_type %in% "PET MRI" ~ img_type_levs[4],
        image_scan_type %in% "Bone Scan" ~ img_type_levs[5],
        image_scan_type %in% "Other Nuclear Medicine Scan" ~ img_type_levs[6],
        T ~ "error in image types"
      )
    ) %>%
    pivot_wider(names_from = image_scan_type,
                values_from = n) 
  dft_count_img %<>%
    mutate(
      n_scan_total = (dft_count_img %>% 
                        select(starts_with("n_scan")) %>% 
                        rowSums(na.rm = T))
    )
  
  dft_count_reg <- reg_dat %>% 
    inner_join(
      select(ca_ind_dat, record_id, ca_seq),
      .,
      by = c("record_id", "ca_seq")
    ) %>%
    count(record_id) %>%
    rename(n_regimens = n)
  
  dft_count_all <- full_join(
    dft_count_reg,
    dft_count_img,
    by = "record_id"
  ) %>%
    full_join(.,
              dft_count_med_onc,
              by = "record_id") %>%
    # to get people who have NONE of these observations, we need:
    left_join(
      select(dft_ca_ind, record_id), # leaving off ca_seq for this cohort.
      ., 
      by = "record_id"
    ) %>%
    mutate(
      across(
        .cols = matches("n_.*"),
        .fns = (function (x) if_else(is.na(x), 0, x))
      )
    )
  
  return(dft_count_all)
  
} 