library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

dft_lot <- readr::read_rds(
  here('data', 'dmet', 'lines_of_therapy', 'lot.rds')
)

dft_lot %<>% filter(line_therapy %in% c(1, 2, 3))


# before we do this, we need to fix the row which has carbo and cis both.
# it happens to have a 0 day duration of cisplatin, so I'm declaring this a
#   mistake.

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


style_help <- function(ft) {
  ft %>%
    theme_box(.) %>%
    fontsize(size = 6, part = "all") %>%
    valign(part = "all", valign = "top") %>%
    padding(part = 'body', padding = 1) %>%
    merge_v(j = c(1)) %>%
    autofit(.)
}

prop_print <- dft_lot_cat %>%
  mutate(pct = round(fraction * 100, 1)) %>%
  select(line_therapy, lot_cat, n, pct) %>%
  arrange(line_therapy, desc(n)) %>%
  flextable(.) %>%
  style_help(.)

reg_table <- dft_lot %>%
  count(lot_cat, regimen_drugs, name = "exposures") %>%
  arrange(lot_cat, desc(exposures)) %>%
  flextable(.) %>%
  style_help

reg_table_split <- dft_lot %>%
  count(line_therapy, lot_cat, regimen_drugs, name = "exposures") %>%
  arrange(line_therapy, lot_cat, desc(exposures)) %>%
  flextable(.) %>%
  style_help(.) %>%
  merge_v(j = c(1, 2))

save_as_docx(
  values = list(
    `n/pct for Line of therapy plots` = prop_print,
    `Exposure category totals (1L through 3L)` = reg_table,
    `Exposure categories (split by line)` = reg_table_split
  ),
  path = here('output', 'other', 'asco_gu_2025_lot_categories.docx')
)


gg_lot_cat <- ggplot(
  dft_lot_cat,
  aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 4,
    xmin = 3,
    fill = lot_cat
  )
) +
  theme_void() +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(1, 4)) +
  facet_wrap(vars(line_lab)) +
  scale_fill_jama(name = "Category")


# Request (also I like this generally):  Do labels inside or outside each donut slice.

dft_lot_cat %<>%
  mutate(
    lot_cat_short = case_when(
      lot_cat %in% "investigational (masked)" ~ "inv",
      lot_cat %in% "immunotherapy" ~ "io",
      lot_cat %in% "cisplatin-based" ~ "cis",
      lot_cat %in% "carboplatin-based" ~ "carbo",
      lot_cat %in% "taxane monotherapy" ~ "tax",
      lot_cat %in% "other" ~ "oth",
      T ~ "CAT ERROR"
    ),
    y_avg = (ymin + ymax) / 2
  )


gg_lot_cat_2 <- ggplot(
  dft_lot_cat,
  aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 4,
    xmin = 3,
    fill = lot_cat
  )
) +
  theme_void() +
  geom_rect() +
  geom_text(
    aes(x = 5, y = y_avg, color = lot_cat, label = lot_cat_short),
    size = 3,
    vjust = 0.5,
    hjust = 0.5
  ) +
  coord_polar(theta = "y") +
  xlim(c(1, 5)) +
  facet_wrap(vars(line_lab)) +
  scale_fill_jama(name = "Category") +
  scale_color_jama(name = "Category")


dft_lot %>%
  count(regimen_drugs, lot_cat) %>%
  arrange(n, lot_cat) %>%
  readr::write_csv(
    x = .,
    file = here('output', 'fig', 'asco_gu_2025', 'lot_categories.csv')
  )

ggsave(
  plot = gg_lot_cat_2,
  filename = here('output', 'fig', 'asco_gu_2025', 'lot_cat_2.pdf'),
  width = 7,
  height = 3
)

ggsave(
  plot = gg_lot_cat,
  filename = here('output', 'fig', 'asco_gu_2025', 'lot_cat.pdf'),
  width = 7,
  height = 3
)

# library(webr) # for donut - it's cool but not quite what I need.
