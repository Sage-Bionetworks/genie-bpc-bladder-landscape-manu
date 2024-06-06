# Note: This script is modified from one generously provided by Konrad Stopsack.
library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

pearl_groups_raw <- readr::read_delim(
  here('data-raw', 'manual', 'dna_repair_genes_pearl.csv'),
  delim = ';'
)

pearl_groups <- pearl_groups_raw

pearl_groups %<>%
  mutate(
    pathway = case_when(
      str_detect(string = DDR, pattern = "Associated") & 
        str_detect(string = Process, pattern = "Checkpoint") ~
        "Checkpoint",
      str_detect(string = DDR, pattern = "Associated") ~
        "Other",
      str_detect(string = DDR, pattern = "probable") ~
        "Probable",
      TRUE ~ str_remove(
        string = Pathway1,
        pattern = "\\s*\\([^\\)]+\\)")
    )
  )

# As stated in Konrad's code the only gene missing a pathway now is H2AFZ.
pearl_groups %<>%
  filter(!is.na(pathway)) %>%
  select(pathway, gene = GeneID) %>% 
  distinct()

# This is the problem we're left with:
# pearl_groups %>% count(gene, sort = T) 

pearl_groups_distinct <- pearl_groups |>
  mutate(
    pathway = fct_relevel(
      factor(pathway),
      # priority of gene->pathway assignment:
      "MMR", "HR", "FA", "NHEJ", "BER", "NER")
  ) |>
  arrange(pathway) |>
  distinct(gene, .keep_all = TRUE) %>%
  mutate(pathway = as.character(pathway))

readr::write_rds(
  pearl_groups_distinct,
  here('data', 'genomic', 'pearl_pathways.rds')
)




