library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dir_out <- here('data', 'genomic', 'ddr_def_compare')

sample_main <- readr::read_tsv(
  here('data-raw', 'genomic', 'main_genie', 'data_clinical_sample.txt'),
  skip = 4
) %>%
  rename_all(tolower)

bladder_samp_mg <- readr::read_rds(
  here('data', 'genomic', 'main_genie', 'bladder_samples_mg.rds')
)

sample_main %<>% filter(sample_id %in% bladder_samp_mg)

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

# One more:  The genes covered in 100% of panels in GENIE BPC, and included in
#    at least one of the four UC papers discussing this.

gp_all <- readr::read_rds(
  here('data', 'genomic', 'gene_panel_all.rds')
)

uc_ddr_author_genes <- ddr_def %>%
  filter(!str_detect(source, 'pearl')) %>%
  pull(gene) %>%
  unique(.)

# There are some of these genes that appear in 0% of our panels:
# gp_all %>%
#   pull(hugo) %>%
#   unique(.) %>%
#   as.character(.) %>%
#   setdiff(uc_ddr_author_genes, .)

gp_ddr_coverage <- gp_all %>%
  group_by(cpt_seq_assay_id) %>%
  summarize(
    gene_overlap = list(intersect(as.character(hugo), # genes in this panel
                                  uc_ddr_author_genes)), # genes in the union of all papers about this.
    .groups = 'drop'
  ) %>% 
  mutate(overlap_size = purrr::map_dbl(.x = gene_overlap, .f = \(x) length(unique(x))))

genie_intersect_panel <- gp_ddr_coverage %>%
  # ad hoc sloppy decision to exclude some of the outlier-tiny panels:
  filter(overlap_size >= 5) %>%
  pull(gene_overlap) %>%
  purrr::reduce(.x = ., .f = intersect)

ddr_def <- bind_rows(
  ddr_def,
  tibble(
    source = 'GENIE int 2024',
    gene = genie_intersect_panel,
  ),
  # also adding in definition we used for the ASCO abstract:
  tibble(
    source = '2025 ASCO GU abstract',
    gene = c('ERCC2', 'ERCC5', 
             'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
             'ATR', 'FANCC')
  )
)





ddr_def_plot <- ddr_def %>%
  filter(!str_detect(source, 'pearl')) %>%
  mutate(year = readr::parse_number(source)) %>%
  arrange(year) %>%
  mutate(source_f = fct_inorder(source))

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
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    legend.position = 'bottom',
    axis.title.y = element_blank()
  )

# gg_def_compare

readr::write_rds(
  gg_def_compare,
  here(dir_out, 'gg_def_compare_heatmap.rds')
)



ddr_def_nest <- ddr_def %>%
  select(source, gene) %>%
  nest(.by = source) %>%
  mutate(
    gene_list = purrr::map(.x = data, .f = \(x) pull(x,gene))
  ) %>%
  select(-data)

first_sample_main <- sample_main %>% 
  group_by(patient_id) %>%
  mutate(
    age_at_seq_report_days_num = as.numeric(age_at_seq_report_days)
  ) %>%
  arrange(seq_year, age_at_seq_report_days_num) %>%
  slice(1) %>%
  ungroup(.)


# Example of using once:    
# count_pts_gene_list(
#   sample_df = first_sample_main,
#   alt_df = alt_main,
#   gene_list = pull(slice(ddr_def_nest, 5), gene_list)[[1]]
# ) %>%
#   count(any_alt)

ddr_def_nest %<>%
  mutate(
    alt_binary = purrr::map(
      .x = gene_list,
      .f = \(x) count_pts_gene_list(
        sample_df = first_sample_main,
        alt_df = alt_main,
        gene_list = x
      )
    ),
    alt_binary_onco = purrr::map(
      .x = gene_list,
      .f = \(x) count_pts_gene_list(
        sample_df = first_sample_main,
        alt_df = filter(alt_main, oncogenic %in% c("Likely Oncogenic", "Oncogenic")),
        gene_list = x
      )
    )
  )

pt_counts_ddr_def <- ddr_def_nest %>%
  select(source, alt_binary, alt_binary_onco) %>%
  pivot_longer(
    cols = c(alt_binary, alt_binary_onco),
    names_to = 'onco_filter',
    values_to = 'pt_df'
  ) %>%
  mutate(
    str = purrr::map_chr(
      .x = pt_df,
      .f = \(x) {x %>%
        summarize(
          n = sum(any_alt),
          d = n()
        ) %>%
        mutate(str = glue("{round((n/d)*100)}% ({n})")) %>%
        pull(str)
      }
    ),
    n_pts = purrr::map_dbl(.x = pt_df, .f = nrow)
  )

pt_counts_ddr_def %<>%
  mutate(
    data = 'GENIE-UC',
    onco_filter = case_when(
      str_detect(onco_filter, 'onco') ~ T,
      T ~ F
    )
  ) %>%
  select(data, onco_filter, source, str, n_pts)
      






# Load in the cohort file from BPC and repeat this process for BPC.
cpt <- readr::read_rds(here('data', 'cohort', 'cpt.rds'))
ca_ind <- readr::read_rds(here('data', 'cohort', 'ca_ind.rds'))

alt_bpc <- readr::read_rds(
  here('data', 'genomic', 'alterations.rds')
)

# just checking:
# setdiff(cpt$cpt_genie_sample_id, sample_main$sample_id)
ddr_def_nest_bpc <- ddr_def_nest %>% select(source, gene_list)

first_sample_bpc <- first_sample_main %>%
  filter(patient_id %in% unique(cpt$record_id))

# Sidenote here:  If you try to subset by sample_id above you get different samples (probably the order is switched after curation or extra detail provided by curation), and a few people go missing as a result.
# Here are the people who differ:
# first_cpt <- get_first_cpt(ca_ind, cpt, include_sample_id = T)
# people_not_in_main_genie <- first_cpt %>%
#   filter(
#     cpt_genie_sample_id %in% setdiff(
#       first_cpt$cpt_genie_sample_id, 
#       first_sample_main$sample_id
#     )
#   ) %>%
#   pull(record_id) 
# filter(first_sample_main, patient_id %in% people_not_in_main_genie)


ddr_def_nest_bpc %<>%
  mutate(
    alt_binary = purrr::map(
      .x = gene_list,
      .f = \(x) count_pts_gene_list(
        sample_df = first_sample_bpc,
        alt_df = alt_bpc,
        gene_list = x
      )
    ),
    alt_binary_onco = purrr::map(
      .x = gene_list,
      .f = \(x) count_pts_gene_list(
        sample_df = first_sample_bpc,
        alt_df = filter(alt_bpc, oncogenic %in% c("Likely Oncogenic", "Oncogenic")),
        gene_list = x
      )
    )
  )

pt_counts_ddr_def_bpc <- ddr_def_nest_bpc %>%
  select(source, alt_binary, alt_binary_onco) %>%
  pivot_longer(
    cols = c(alt_binary, alt_binary_onco),
    names_to = 'onco_filter',
    values_to = 'pt_df'
  ) %>%
  mutate(
    str = purrr::map_chr(
      .x = pt_df,
      .f = \(x) {x %>%
          summarize(
            n = sum(any_alt),
            d = n()
          ) %>%
          mutate(str = glue("{round((n/d)*100)}% ({n})")) %>%
          pull(str)
      }
    ),
    n_pts = purrr::map_dbl(.x = pt_df, .f = nrow)
  )

pt_counts_ddr_def_bpc %<>%
  mutate(
    data = 'GENIE-BPC',
    onco_filter = case_when(
      str_detect(onco_filter, 'onco') ~ T,
      T ~ F
    )
  ) %>%
  select(data, onco_filter, source, str, n_pts)

pt_counts_ddr_def <- bind_rows(
  pt_counts_ddr_def,
  pt_counts_ddr_def_bpc
)

readr::write_rds(
  pt_counts_ddr_def, 
  here(dir_out, 'pt_count_tidy.rds')
)

ft_ddr_def <- pt_counts_ddr_def %>%
  mutate(
    col_title = paste0(data, '|',
                      '(n = ', n_pts, ')', '|',
                      if_else(onco_filter, "Oncogenic", 'Any alt.'))
  ) %>%
  select(source, col_title, str) %>%
  pivot_wider(
    names_from = 'col_title',
    values_from = 'str'
  ) 

new_header <- names(ft_ddr_def) %>% 
  tibble(col_keys = .) %>%
  separate_wider_delim(names = c('one', 'two', 'three'),
                       col_keys, delim = '|', cols_remove = F,
                       too_few = 'align_start')

ft_ddr_def <- ft_ddr_def %>%
  flextable(.) %>%
  flextable::set_header_df(mapping = new_header, key = 'col_keys') %>%
  theme_booktabs(.) %>%
  autofit(.)

readr::write_rds(
  ft_ddr_def,
  here(dir_out, 'ft_ddr_def.rds')
)






# Figure: pull in the main genie panel data and look at coverage of the DDR genes.
panels_main <- count(sample_main, seq_assay_id)
panels_main <- panels_main %>%
  mutate(bpc_panel = seq_assay_id %in% unique(cpt$cpt_seq_assay_id))
gi <- data.table::fread(
  here('data-raw', 'genomic', 'main_genie', 'genomic_information.txt')
)
gi %<>% 
  rename_all(tolower) %>%
  filter(seq_assay_id %in% panels_main$seq_assay_id) %>%
  # the best information I have currently is that includeInPanel = F means the gene is not tested.  Frustratingly we can't seem to get a real answer about that.
  filter(includeinpanel) 

genes_in_ddr_cov_plot <- ddr_def %>% 
  filter(!str_detect(source, "pearl")) %>% 
  pull(gene) %>% 
  unique

ddr_panel_cov <- gi %>%
  filter(hugo_symbol %in% genes_in_ddr_cov_plot)

ddr_panel_cov <- ddr_panel_cov %>% 
  group_by(seq_assay_id, hugo_symbol) %>%
  summarize(tested = T, .groups = 'drop') %>%
  mutate(hugo = factor(hugo_symbol, levels = genes_in_ddr_cov_plot)) %>%
  select(seq_assay_id, hugo, tested) %>%
  complete(seq_assay_id, hugo, fill = list(tested = F))

ddr_panel_cov %<>%
  left_join(., panels_main, by = 'seq_assay_id') %>%
  mutate(assay_str = case_when(
    bpc_panel ~ glue("<span style = 'color:#B12F00;'>{seq_assay_id} (m={n})</span>"),
    T ~ glue("{seq_assay_id} (m={n})")
    )
  ) %>%
  arrange(desc(bpc_panel), desc(n), seq_assay_id) %>% 
  mutate(assay_str = fct_inorder(assay_str),
         assay_str = fct_rev(assay_str)) %>% # rev for plotting.
  rename(gene = hugo)

gg_ddr_panel <- ggplot(
  ddr_panel_cov,
  aes(x = gene, y = assay_str, fill = tested)
) + 
  geom_tile(color = 'gray60') + 
  scale_fill_manual(values = c("gray90", '#366885')) + 
  scale_y_discrete(expand = c(0,0), position = 'right') + 
  scale_x_discrete(expand = c(0,0)) + 
  theme_classic() +
  labs(
    title = "GENIE Urothelial Carinoma panel coverage of DDR genes",
    subtitle = "GENIE BPC panels noted in <span style = 'color:#B12F00;'>orange</span>, m = main GENIE samples."
  ) + 
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_markdown(size = 6),
    legend.position = 'bottom',
    axis.title.y = element_blank(),
    title = element_markdown()
  )

readr::write_rds(
  gg_ddr_panel,
  here(dir_out, 'gg_ddr_panel.rds')
)

            


