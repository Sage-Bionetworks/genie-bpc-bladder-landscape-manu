library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)


sample_main <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_sample.txt'),
  skip = 4
) %>%
  rename_all(tolower)

alt_main <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'alterations.rds')
)

ddr_def <- readr::read_csv(
  here('data-raw', 'manual', 'ddr_gene_papers.csv')
)

ddr_def %<>% mutate(gene = str_trim(gene))

# add the pearl definitions because I don't work at MSK:
ddr_pearl <- readr::read_rds(
  here('data', 'genomic', 'pearl_pathways.rds')
)

ddr_def_pearl_addons <- bind_rows(
  (ddr_pearl %>%
    mutate(source = 'pearl 2015 full') %>%
    select(source, gene) %>%
    distinct(.)),
  (ddr_pearl %>%
     # selecting only the pathways that seem to have caught the attention of 
     #   urothelial carcinoma researchers to date.
     filter(pathway %in% c('NER', 'HR', 'Checkpoint', 'FA', 'BER')) %>%
     mutate(source = 'pearl 2015 select') %>%
     select(source, gene) %>%
     distinct(.))
)

ddr_def <- bind_rows(ddr_def, ddr_def_pearl_addons)

ddr_def_plot <- ddr_def %>%
  filter(!str_detect(source, 'pearl')) %>%
  mutate(year = readr::parse_number(source)) %>%
  arrange(year) %>%
  mutate(source_f = fct_inorder(source))

cli_alert_warning("One more I want to add:  A very small subset covered by all GENIE panels")
  
# first just a basic comparison of what these defs cover:
ddr_def_plot <- ddr_def_plot %>%
  select(-`author affiliation`) %>%
  mutate(included = TRUE) %>%
  complete(source_f, gene, fill = list(included = FALSE))

gg_def_compare <- ggplot(
  ddr_def_plot,
  aes(x = gene, y = source_f, fill = included)
) + 
  geom_tile(color = 'gray60') + 
  scale_fill_manual(values = c("gray90", '#366885')) + 
  scale_y_discrete(expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

gg_def_compare

readr::write_rds(
  gg_def_compare,
  here('data', 'genomic', 'gg_def_compare_heatmap.rds')
)



ddr_def %>%
  select(source, gene) %>%
  nest(.by = source) %>%
  mutate(
    gene_list = purrr::map(.x = data, .f = \(x) list(pull(x,gene)))
  ) %>%
  select(-data) %>% View(.)

first_sample_main <- sample_main %>% 
  group_by(patient_id) %>%
  mutate(
    age_at_seq_report_days_num = as.numeric(age_at_seq_report_days)
  ) %>%
  arrange(seq_year, age_at_seq_report_days_num) %>%
  slice(1) %>%
  ungroup(.)


# Next steps:  Make a function to count the number of first samples (two datasets, sample list and alterations) which have at least one ddr mutation based on a list.  Then process that for all the various definitions and count it up.

