# The prompt for this investigation is looking at second line single agent
#   chemotherapies.  They were interested in taxanes and similar.
# Taxanes are paclitaxel, docetaxel, cabizitaxel

library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_output <- here('data', 'survival', 'line23')

dft_ca_ind <- readr::read_rds(here('data', 'cohort', "ca_ind.rds"))
dft_reg <- readr::read_rds(here('data', 'cohort', "reg.rds"))
dft_lot <- readr::read_rds(here('data', 'dmet', 'lines_of_therapy', 'lot.rds'))
dft_cpt <- readr::read_rds(here('data', 'cohort', "cpt_aug.rds"))

dft_agent_23_line <- dft_lot %>% 
  filter(line_therapy %in% c(2,3)) %>%
  filter(!(regimen_drugs %in% "Investigational Regimen")) %>%
  count(regimen_drugs, line_therapy, sort = T) %>%
  mutate(line_therapy = paste0("line_", line_therapy)) %>%
  pivot_wider(names_from = line_therapy, values_from = n) %>%
  filter(line_2 >= 5 | line_3 >= 5)

readr::write_rds(
  x = dft_agent_23_line,
  file = here(dir_output, 'agent_23_line.rds')
)



# For now I really don't know what the question is.  But we definitely
#   have the numbers to do KM curves for single agent taxanes.

dft_line23_tax <- dft_reg %>%
  filter(regimen_drugs %in% c("Docetaxel", "Paclitaxel")) %>%
  left_join(
    .,
    select(dft_lot, record_id, ca_seq, regimen_number, line_therapy),
    by = c('record_id', 'ca_seq', 'regimen_number')
  ) %>%
  filter(line_therapy %in% c(2,3))

dft_line23_tax %<>%
  arrange(line_therapy, regimen_drugs) %>%
  mutate(
    reg_line = paste0(regimen_drugs, "[", line_therapy, "L]"),
    reg_line = fct_inorder(reg_line)
  ) 

dft_first_cpt <- get_first_cpt(dft_ca_ind, dft_cpt) %>% 
  rename(dx_first_cpt_rep_yrs = dx_cpt_rep_yrs)

dft_line23_tax %<>%
  left_join(
    .,
    dft_first_cpt,
    by = c("record_id", 'ca_seq')
  )

dft_line23_tax %<>%
  mutate(
    # it's expected that many of these will be negative:
    reg_start_first_cpt_rep_yrs = dx_first_cpt_rep_yrs - dx_reg_start_int_yrs
  )

readr::write_rds(
  x = dft_line23_tax,
  file = here(dir_output, 'line23_tax_surv_data.rds')
)




dft_line23_tax %<>% 
  remove_trunc_gte_event(
    trunc_var = 'reg_start_first_cpt_rep_yrs',
    event_var = 'tt_os_g_yrs'
  )

dft_line23_tax %<>% 
  mutate(
    reg_start_first_cpt_rep_yrs = ifelse(
      reg_start_first_cpt_rep_yrs < 0, 0, reg_start_first_cpt_rep_yrs
    )
  )

surv_obj_os_line23_taxane <- with(
  dft_line23_tax,
  Surv(
    time = reg_start_first_cpt_rep_yrs,
    time2 = tt_os_g_yrs,
    event = os_g_status
  )
)

gg_os_line23_taxane <- plot_one_survfit(
  dat = dft_line23_tax,
  surv_form = surv_obj_os_line23_taxane ~ reg_line,
  plot_title = "OS from 2L/3L taxane",
  plot_subtitle = "Adjusted for (independent) delayed entry",
  pal = c("#ee99aa", "#994455", "#6699cc", "#004488"),
  x_breaks = seq(0, 100, by = 0.5),
  risktable_prop = 0.4
)

readr::write_rds(
  x = gg_os_line23_taxane,
  file = here(dir_output, 'gg_line23_taxane.rds')
)

# A little extra but I like the idea:
gg_line23_legend <- tribble(
  ~drug, ~line,
  "Paclitaxel", "2L",
  "Paclitaxel", "3L",
  "Docetaxel", "2L",
  "Docetaxel", "3L"
) %>%
  mutate(
    drug = fct_inorder(drug),
    line = fct_inorder(line),
    dl = fct_inorder(paste0(drug, line))
  ) %>%
  mutate(dl = paste0(drug, line)) %>%
  ggplot(., aes(x = line, y = drug, fill = dl)) +
    geom_tile() + 
    scale_fill_manual(
      values = c("#ee99aa", "#994455", "#6699cc", "#004488")[c(1,3,2,4)]
    ) + 
  theme_void() + 
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11)
  )
  
readr::write_rds(
  x = gg_line23_legend,
  file = here(dir_output, 'gg_line23_legend.rds')
)


  
