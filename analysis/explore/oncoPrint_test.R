# Make a complex heatmap oncoprint.
library(magrittr)
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

mat = read.table(
  textConnection(
    "s1,s2,s3
  g1,snv;indel,snv,indel
  g2,,snv;indel,snv
  g3,snv,,indel;snv"
  ),
  row.names = 1,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
)
mat = as.matrix(mat)
mat

get_type_fun = function(x) strsplit(x, ";")[[1]]
get_type_fun(mat[1, 1])
get_type_fun(mat[1, 2])

# ... and here the details of this nonsense leave me.

mat = read.table(
  system.file(
    "extdata",
    package = "ComplexHeatmap",
    "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"
  ),
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]

col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(
      x,
      y,
      w - unit(2, "pt"),
      h - unit(2, "pt"),
      gp = gpar(fill = "#CCCCCC", col = NA)
    )
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(
      x,
      y,
      w - unit(2, "pt"),
      h - unit(2, "pt"),
      gp = gpar(fill = col["HOMDEL"], col = NA)
    )
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(
      x,
      y,
      w - unit(2, "pt"),
      h - unit(2, "pt"),
      gp = gpar(fill = col["AMP"], col = NA)
    )
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(
      x,
      y,
      w - unit(2, "pt"),
      h * 0.33,
      gp = gpar(fill = col["MUT"], col = NA)
    )
  }
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  HOMDEL = alter_graphic("rect", fill = col["HOMDEL"]),
  AMP = alter_graphic("rect", fill = col["AMP"]),
  MUT = alter_graphic("rect", height = 0.33, fill = col["MUT"])
)

column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
heatmap_legend_param = list(
  title = "Alternations",
  at = c("HOMDEL", "AMP", "MUT"),
  labels = c("Deep deletion", "Amplification", "Mutation")
)
oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)

# get the sample order from cbio.
sample_order = scan(
  paste0(
    system.file("extdata", package = "ComplexHeatmap"),
    "/sample_order.txt"
  ),
  what = "character"
)
op <- oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  row_order = 1:nrow(mat),
  column_order = sample_order,
  # These are almost comically bad options to set - removing samples with no alterations is the columns.
  # Same for the empty rows - removing genes that are never altered.
  # but we'll do it for demo purposes.
  remove_empty_columns = T,
  remove_empty_rows = T,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param
)
