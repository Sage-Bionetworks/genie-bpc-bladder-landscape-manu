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

asco_2025_genes <- c('ERCC2', 'ERCC5', 
                     'BRCA1', 'BRCA2', 'RECQL4', 'RAD51C', 'ATM', 
                     'ATR', 'FANCC')

ddr_def <- bind_rows(
  ddr_def,
  tibble(
    source = 'GENIE int 2024',
    gene = genie_intersect_panel,
  ),
  # also adding in definition we used for the ASCO abstract:
  tibble(
    source = '2025 ASCO GU abstract',
    gene = asco_2025_genes,
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

ddr_def_plot %<>%
  group_by(source_f) %>%
  mutate(num_genes = sum(included)) %>%
  mutate(source_f_gene_count = forcats::fct_inorder(
    glue('{source_f} (p={num_genes})')
  )) %>%
  select(-num_genes)

gg_def_compare <- ggplot(
  ddr_def_plot,
  aes(x = gene, y = source_f_gene_count, fill = included)
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
cpt <- readr::read_rds(here('data', 'cohort', 'cpt_aug.rds'))
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

# just going to redo this here:

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

# Just going to redo the names here:
ft_ddr_def <- ddr_def_nest %>% 
  mutate(
    gene_counts = purrr::map_dbl(
      .x = gene_list,
      .f = \(x) {
        length(unique(x))
      }
    ),
    source_f_gene_counts = fct_inorder(
      glue('{source} (p={gene_counts})')
    )
  ) %>%
  select(source, source_f_gene_counts) %>%
  # doing a join and replace here is weird, but ensures we don't mislabel things.
  left_join(
    ft_ddr_def,
    .,
    by = 'source'
  ) %>%
  select(-source) %>%
  rename(source = source_f_gene_counts) %>%
  relocate(source)
  

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

gene_panel_gene_order <- ddr_def %>% 
  group_by(gene) %>%
  summarize(
    in_asco_panel = any(str_detect(source, "ASCO"))
  ) %>%
  arrange(desc(in_asco_panel), gene) %>%
  pull(gene)
    
ddr_panel_cov %<>%
  mutate(gene = factor(gene, levels = gene_panel_gene_order))

gg_ddr_panel <- ggplot(
  ddr_panel_cov,
  aes(x = gene, y = assay_str, fill = tested)
) + 
  geom_tile(color = 'gray60') + 
  scale_fill_manual(values = c("gray90", '#366885')) + 
  scale_y_discrete(expand = c(0,0), position = 'right') + 
  scale_x_discrete(expand = c(0,0), position = 'top') + 
  theme_classic() +
  labs(
    title = "GENIE Urothelial Carinoma panel coverage of DDR genes",
    subtitle = "GENIE BPC panels noted in <span style = 'color:#B12F00;'>orange</span>, m = main GENIE samples."
  ) + 
  theme(
    axis.text.x.top = element_text(size = 6, angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = element_markdown(size = 6),
    legend.position = 'bottom',
    axis.title.y = element_blank(),
    title = element_markdown()
  )

readr::write_rds(
  gg_ddr_panel,
  here(dir_out, 'gg_ddr_panel.rds')
)

            

# Some genes are not covered but alterations are reported (fusions)
# Some genes are tested but never positive (obvious)
genes_observed_or_tested <- unique(c(
  as.character(alt_bpc$hugo),
  as.character(gp_all$hugo)
))

onco_alt_flags <- alt_bpc %>% 
  filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic")) %>%
  group_by(sample_id, hugo) %>%
  summarize(onco_alt = T, .groups = 'drop') %>%
  mutate(hugo = factor(hugo, levels = sort(genes_observed_or_tested))) %>%
  complete(sample_id, hugo, fill = list(onco_alt = F)) %>%
  pivot_wider(
    names_from = 'hugo',
    values_from = 'onco_alt',
    values_fill = F
  )

readr::write_rds(
  onco_alt_flags,
  here('data', 'genomic', 'ddr_def_compare', 'onco_alt_flags.rds')
)

bpc_asco_2025_panel <- first_sample_bpc %>%
  select(sample_id, patient_id) %>%
  left_join(
    .,
    select(onco_alt_flags, sample_id, all_of(asco_2025_genes)),
    by = 'sample_id'
  ) %>%
  mutate(
    across(
      .cols = -c(sample_id, patient_id),
      .fns = \(x) if_else(is.na(x), F, x)
    )
  )

readr::write_rds(
  bpc_asco_2025_panel,
  here('data', 'genomic', 'ddr_def_compare', 'bpc_asco_2025_panel.rds')
)
  

mmr_genes_from_abstracts <- c('MLH1', 'MSH2', 'MSH6', 'PMS1', 'PMS2')

bpc_onco_mmr_flags <- first_sample_bpc %>%
  select(sample_id, patient_id) %>%
  left_join(
    .,
    select(
      onco_alt_flags, sample_id, 
      all_of(mmr_genes_from_abstracts)
    ),
    # this didn't show anything:
    # any_of(setdiff(
    #   sort(pull(filter(ddr_pearl, pathway %in% "MMR"), gene)),
    #   mmr_genes_from_abstracts
    # ))),
    by = 'sample_id'
  ) %>%
  mutate(
    across(
      .cols = -c(sample_id, patient_id),
      .fns = \(x) if_else(is.na(x), F, x)
    )
  )

readr::write_rds(
  bpc_onco_mmr_flags,
  here('data', 'genomic', 'ddr_def_compare', 'bpc_onco_mmr_flags.rds')
)








# One more thing:  Comparison of TMB for DDR/not
first_sample_bpc <- get_first_cpt(ca_ind, cpt, include_sample_id = T) %>%
  select(
    patient_id = record_id, sample_id = cpt_genie_sample_id
  )

tmb_ddr <- count_pts_gene_list(
  first_sample_bpc, 
  alt_df = filter(alt_bpc, oncogenic %in% c("Likely Oncogenic", "Oncogenic")),
  gene_list = asco_2025_genes
) %>%
  rename(ddr = any_alt) %>%
  mutate(ddr = as.logical(ddr))

tmb_mmr <- count_pts_gene_list(
  first_sample_bpc, 
  alt_df = filter(alt_bpc, oncogenic %in% c("Likely Oncogenic", "Oncogenic")),
  gene_list = mmr_genes_from_abstracts
) %>%
  rename(mmr = any_alt) %>%
  mutate(mmr = as.logical(mmr))

tmb_comp <- full_join(
  tmb_ddr,
  tmb_mmr,
  by = c('patient_id', 'sample_id')
)

tmb_comp %<>%
  left_join(
    .,
    select(cpt, sample_id = cpt_genie_sample_id, n_mut, tmb_Mb, tmb_Mb_onco),
    by = 'sample_id'
  )

readr::write_rds(
  tmb_comp,
  here('data', 'genomic', 'ddr_def_compare', 'tmb_comp.rds')
)


# There's a problem with the comparison of DDR to not for TMB (same for MMR).
# Having a DDR alteration strictly bounds the number of mutations at 1, which
#   directly impacts TMB.  So our preditor is a part of the outcome definition, which just won't do.
# We will build a better control, and stack it into the data to show this.
tmb_stack_ddr <- bind_rows(
  (tmb_comp %>%
    mutate(
      grp = case_when(
        ddr ~ "DDR+",
        T ~ "DDR-"
      )
    )),
  (tmb_comp %>%
     filter(n_mut >= 1 & !ddr) %>%
     mutate(
       grp = "DDR-, >0 mut"
     ))
)

tmb_stack_ddr %<>% 
  select(patient_id, grp, tmb_Mb, tmb_Mb_onco) %>%
  mutate(grp = forcats::fct_inorder(grp))

readr::write_rds(
  tmb_stack_ddr,
  here('data', 'genomic', 'ddr_def_compare', 'tmb_stack_ddr.rds')
)


# Repeat for mmr:

tmb_stack_mmr <- bind_rows(
  (tmb_comp %>%
     mutate(
       grp = case_when(
         mmr ~ "MMR+",
         T ~ "MMR-"
       )
     )),
  (tmb_comp %>%
     filter(n_mut >= 1 & !mmr) %>%
     mutate(
       grp = "MMR-, >0 mut"
     ))
)

tmb_stack_mmr %<>% 
  select(patient_id, grp, tmb_Mb, tmb_Mb_onco) %>%
  mutate(grp = forcats::fct_inorder(grp))

readr::write_rds(
  tmb_stack_mmr,
  here('data', 'genomic', 'ddr_def_compare', 'tmb_stack_mmr.rds')
)

  