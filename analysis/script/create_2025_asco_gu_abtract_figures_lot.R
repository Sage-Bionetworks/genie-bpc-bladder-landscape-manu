

library(purrr); library(fs); library(here);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dft_lot %<>% filter(line_therapy %in% c(1,2,3))

lot_cat_lev <- c(
  'investigational (masked)',
  'io',
  'platinum chemo',
  'taxane monotherapy',
  'pemetrexed monotherapy',
  'other'
)
  

categorize_lines <- function(reg_vec, factorize = T) {
  reg_vec <- tolower(reg_vec)
  
  lot_cat_lev <- c(
    'investigational (masked)',
    'immunotherapy',
    'platinum chemo',
    'taxane monotherapy',
    'pemetrexed monotherapy',
    'other'
  )
  
  rtn <- case_when(
    # after a long discussion we decided to include nivolumab here even though
    #   it's not in the 2L io figure.
    str_detect(reg_vec, 'investigation') ~ lot_cat_lev[1],
    str_detect(reg_vec, 'atezoliz|pembroliz|nivolum') ~ lot_cat_lev[2],
    str_detect(reg_vec, 'carbo|cis') ~ lot_cat_lev[3],
    reg_vec %in% c('docetaxel', 'paclitaxel') ~ lot_cat_lev[4],
    reg_vec %in% c('docetaxel', 'paclitaxel') ~ lot_cat_lev[5],
    T ~ lot_cat_lev[6]
  )
  
  if (factorize) { 
    rtn <- factor(rtn, levels = lot_cat_lev)
  }
  
  return(rtn)
}


dft_lot %<>% 
  mutate(
    lot_cat = categorize_lines(regimen_drugs)
  )

dft_lot_cat <- dft_lot %>% 
  count(line_therapy, lot_cat) %>%
  group_by(line_therapy) %>%
  mutate(
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = cumsum(lag(fraction, default = 0))
  ) %>%
  ungroup(.)

dft_lot_cat %<>%
  group_by(line_therapy) %>%
  mutate(
    n_people_line = sum(n), # each person can only have one 1L, 2L, etc.
    line_lab = glue('{line_therapy}L (n={n_people_line})')
  ) %>%
  ungroup(.)

library(ggsci)

gg_lot_cat <- ggplot(
  dft_lot_cat, 
  aes(
    ymax = ymax,  ymin = ymin, 
    xmax = 4,  xmin = 3, 
    fill=lot_cat
  )
) +
  theme_void() + 
  geom_rect() +
  coord_polar(theta="y") +
  xlim(c(1, 4)) + 
  facet_wrap(vars(line_lab)) + 
  scale_fill_jama(name = "Category")

dft_lot %>% 
  count(regimen_drugs, lot_cat) %>%
  arrange(n, lot_cat) %>%
  readr::write_csv(
    x = .,
    file = here('output', 'fig', 'asco_gu_2025', 'lot_categories.csv')
  )

ggsave(
  plot = gg_lot_cat,
  filename = here('output', 'fig', 'asco_gu_2025', 'lot_cat.pdf'),
  width = 7, height = 3
)
       
  
  

# library(webr) # for donut - it's cool but not quite what I need.
