# Description:  Calculate the panel size and TMB for each sample.
#   This code is adapated from Alex Baras's code used in early genie releases.


library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

# There was an error with 2.4 for dbplyr when calling getBM() below. 
# This fixed it:
# devtools::install_version("dbplyr", version = "2.3.4")

library(biomaRt)
library(VariantAnnotation) # bioconductor.

dft_cpt <- readr::read_rds(
  here('data', 'cohort', 'cpt.rds')
)

# Loading the oncogenic annotated version - still contains non-oncogenic variants.
dft_maf <- readr::read_tsv(
  here('data', 'genomic', 'mut_onco.txt')
)

dft_gp_all <- readr::read_rds(
  here('data', 'genomic', 'gene_panel_all.rds')
)
vec_relevant_panel <- dft_gp_all %>% 
  pull(cpt_seq_assay_id) %>%
  unique(.)

# This is a more up-to-date version of the bed file, but it shouldn't matter
#   because we expect the panel coverage to be stable over time.
dft_bed <- readr::read_tsv(
  here('data-raw', 'genomic', 'genie_combined_14_1.bed')
)
dft_bed %<>% 
  filter(SEQ_ASSAY_ID %in% vec_relevant_panel)





# Part 1:  Get the data for each panel.
# I left this as lists for simplicity but it may be cleaner to do list columns
#  in a tibble in the future.

list_bed <- split(dft_bed, factor(dft_bed$SEQ_ASSAY_ID))

list_bed <- lapply(
  list_bed,
  function(x) {
    GR <- GRanges(
      seqnames = Rle(paste0("chr", x$Chromosome)), 
      ranges = IRanges(
        start = x$Start_Position, 
        end = x$End_Position
      )
    )
    seqlevels(GR) <- sort(seqlevels(GR))
    return(GR)
  }
)

# get GRanges for each bed limited to coding in grch37 as per ensemb biomart
# connect to ensembl for human'
ens <- useMart(
  "ensembl", 
  dataset = "hsapiens_gene_ensembl", 
  host = "https://grch37.ensembl.org"
)
# get all protein coding exons
geneAnnot <- getBM(
  attributes = c(
    "ensembl_gene_id", 
    "ensembl_transcript_id", 
    "genomic_coding_start", 
    "genomic_coding_end", 
    "ensembl_exon_id", 
    "external_gene_name", 
    "strand", 
    "start_position", 
    "end_position", 
    "exon_chrom_start", 
    "exon_chrom_end", 
    "chromosome_name"
  ),
  filters = c(
    "biotype", 
    "chromosome_name"
  ),
  values = list(
    A = "protein_coding", 
    B = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
  ),
  mart = ens
)

# restrict to records with none NA genomic coding coordinates
idx <- !(is.na(geneAnnot$genomic_coding_start) | is.na(geneAnnot$genomic_coding_start))
codingGR <- GRanges(
  seqnames = Rle(paste0("chr", geneAnnot$chromosome_name[idx])), 
  ranges = IRanges(
    start = geneAnnot$genomic_coding_start[idx], 
    end = geneAnnot$genomic_coding_end[idx]
  )
)
seqlevels(codingGR) <- sort(seqlevels(codingGR))

# contrust GRanges for coding sequence only
list_bed_exon <- lapply(
  list_bed, 
  function(x) {
    GR <- reduce(
      pintersect(
        findOverlapPairs(
          reduce(x), 
          reduce(codingGR)
        )
      )
    )
    seqlevels(GR) <- sort(seqlevels(GR))
    return(GR)
  }
)


dft_cpt_addon <- dft_cpt %>%
  dplyr::select(
    SAMPLE_ID = cpt_genie_sample_id,
    SEQ_ASSAY_ID = cpt_seq_assay_id,
    CENTER = institution
  )

dft_panel_stats <- tibble(SEQ_ASSAY_ID = names(list_bed)) %>%
  mutate(
    list_bed = list_bed,
    list_bed_exon = list_bed_exon
  )
dft_panel_stats %<>%
  mutate(
    bps = purrr::map_dbl(
      .x = list_bed,
      .f = (function(x) {
        # reduce takes care of overlaps and adjacent ranges.
        # width gets the size of each range - sum is the normal sum.
        sum(width(reduce(x)))
      })
    ),
    bpsExon = purrr::map_dbl(
      .x = list_bed_exon,
      .f = (function(x) {
        sum(width(x))
      })
    )
  )

dft_panel_stats %<>% 
  dplyr::select(-c(list_bed, list_bed_exon))


# Part 2: TMB using the panel data.

tumor_aggr <- dft_maf %>%
  filter(!(Variant_Classification %in% "Silent")) %>%
  dplyr::count(Tumor_Sample_Barcode, name = "numberMutations")

dft_cpt_addon <- left_join(
  dft_cpt_addon,
  tumor_aggr,
  by = c(SAMPLE_ID = "Tumor_Sample_Barcode")
) %>%
  mutate(numberMutations = if_else(is.na(numberMutations), 0, numberMutations)) 

# Update: I'm also interested in what happens we limit to oncogenic mutations.
# This is the same steps as above, with one extra filter and a new name.
tumor_aggr_onco <- dft_maf %>% 
  filter(ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")) %>%
  # Doing the filter for silent mutations should be entirely redundant here.
  dplyr::count(Tumor_Sample_Barcode, name = "numberMutations_onco")

dft_cpt_addon <- left_join(
  dft_cpt_addon,
  tumor_aggr_onco,
  by = c(SAMPLE_ID = "Tumor_Sample_Barcode")
) %>%
  mutate(
    numberMutations_onco = if_else(
      is.na(numberMutations_onco), 
      0, 
      numberMutations_onco
    )
  ) 

dft_cpt_addon %<>%
  left_join(
    .,
    dplyr::select(as_tibble(dft_panel_stats), SEQ_ASSAY_ID, bps, bpsExon),
    by = 'SEQ_ASSAY_ID'
  )

dft_cpt_addon %<>%
  mutate(
    tmb_Mb = numberMutations / (bpsExon / 10^6),
    tmb_Mb_onco = numberMutations_onco / (bpsExon / 10^6)
  )

# Note:  The code we're adapting used <2 as low, 2-16 as mid and >16 as high.





# Part 3:  Add these new calculations back into the cpt file.
dft_cpt_addon %<>%
  dplyr::select(
    cpt_genie_sample_id = SAMPLE_ID,
    n_mut = numberMutations,
    panel_bp = bps,
    panel_bp_exon = bpsExon,
    tmb_Mb,
    tmb_Mb_onco
  )

dft_cpt_aug <- left_join(
  dft_cpt,
  dft_cpt_addon,
  by = 'cpt_genie_sample_id',
  relationship = "one-to-one"
)

readr::write_rds(
  x = dft_cpt_aug,
  file = here('data', 'clin', 'dft_cpt_aug.rds')
)



