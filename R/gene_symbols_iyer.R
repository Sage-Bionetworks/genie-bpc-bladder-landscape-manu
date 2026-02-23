gene_symbols_iyer <- function() {
  # pulled from his 2025 paper https://doi.org/10.1200/PO-24-00938
  unique(
    c(
      # MMR pathway:
      'MLH1',
      'MSH2',
      'MSH6',
      'PMS1',
      'PMS2',
      # NER:
      'ERCC2',
      'ERCC3',
      'ERCC4',
      'ERCC5',
      # HR:
      'BRCA1',
      'MRE11A',
      'NBN',
      'RAD50',
      'RAD51',
      'RAD51B',
      'RAD52',
      'RAD54L',
      # FA:
      'BRCA2',
      'BRIP1',
      'FANCA',
      'FANCC',
      'PALB2',
      'RAD51C',
      # Checkpoint:
      'ATM',
      'ATR',
      'CHEK1',
      'CHEK2',
      'MDC1',
      # Others:
      'POLE',
      'MUTYH',
      'PARP1',
      'REQCQL4'
    )
  )
}
