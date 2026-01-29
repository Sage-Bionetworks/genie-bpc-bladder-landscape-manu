# A consort-like diagram to explain the flow of patients into our cohort.
library(purrr)
library(fs)
library(here)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

# Pull the patient counts from our ddr report definitions table.
pt_counts_ddr_def <- readr::read_rds(
  here('data', 'genomic', 'ddr_def_compare', 'pt_count_tidy.rds')
)

pt_ct_simple <- pt_counts_ddr_def %>%
  group_by(data) %>%
  slice(1) %>%
  ungroup(.)

# These come from lots of different places unfortunately.  Two of them should be fixed so we're just hard coding them in.
coh_flow <- tribble(
  ~step,
  ~n,

  "Main GENIE",
  205739, # 18.2-consortium specifically

  "Main GENIE \\n Oncotree valid for \\n BPC Bladder",
  pt_ct_simple$n_pts[2],

  "BPC Bladder (curated)",
  716,

  "BPC Bladder UC (FAS)",
  pt_ct_simple$n_pts[1]
)

mermaid_input <- paste0(
  "graph TD;\n",
  "  A[\"<center>Main GENIE<br>n = <b>",
  format(coh_flow$n[1], big.mark = ","),
  "</b></center>\"] --> B[\"<center>Main GENIE<br>Oncotree valid for<br>BPC Bladder<br>n = <b>",
  format(coh_flow$n[2], big.mark = ","),
  "</b></center>\"];\n",
  "  B --> C[\"<center>BPC Bladder (curated)<br>n = <b>",
  format(coh_flow$n[3], big.mark = ","),
  "</b></center>\"];\n",
  "  C --> D[\"<center>BPC Bladder UC (FAS)<br>n = <b>",
  format(coh_flow$n[4], big.mark = ","),
  "</b></center>\"];\n"
)

library(plotly)
library(webshot)

coh_mermaid <- DiagrammeR::mermaid(mermaid_input)

# did a manual shot from here.  Just awful I hate this.
