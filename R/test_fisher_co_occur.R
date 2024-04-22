test_fisher_co_occur <- function(
    dat,
    ignore_cols = character(0),
    top = 5,
    alpha = 0.05,
    axis_order = NULL
) {
  dat_top <- get_binary_feature_pos(
    dat,
    ignore_cols = ignore_cols
  )
  
  dat_top %<>% 
    mutate(rank = 1:n()) %>%
    filter(rank <= top)
  
  if (!is.null(axis_order)) {
    if (any( !(dat_top$feature %in% axis_order))) {
      cli_abort("axis_order does not contain all the top genes - needs fixing!")
    }
    
    dat_top %<>%
      mutate(.temp = factor(feature, levels = axis_order)) %>%
      arrange(.temp) %>%
      select(-.temp) 
  }
  
  
  dat_top %<>%
    mutate(
      # just matching the maftools syntax to start:
      feat_lab = paste0(feature, " [", num_pos, "]")
    )
  
  trimmed_dat <- dat %>%
    select(
      -all_of(ignore_cols)
    ) %>%
    select(
      all_of(dat_top$feature)
    )
  
  test_skel <- expand_grid(
    var1 = colnames(dat),
    var2 = colnames(dat)
  ) %>%
    # Trim this list to the upper triangle:
    mutate(
      var1 = ordered(var1, levels = dat_top$feature),
      var2 = ordered(var2, levels = dat_top$feature)
    ) %>%
    filter(
      var1 < var2
    ) %>%
    mutate(
      var1 = as.character(var1),
      var2 = as.character(var2)
    )
  
  
  
  test_skel %<>%
    mutate(
      cont_tab = purrr::map2(
        .x = var1,
        .y = var2,
        .f = (function(v1,v2) {
          cont_tab_variable_pair(
            dat = trimmed_dat,
            var1 = v1,
            var2 = v2
          )
        })
      )
    ) %>%
    unnest(cont_tab)
  
  test_skel %<>%
    mutate(
      fish_result = purrr::pmap(
        .l = list(
          cell_11 = ct_11,
          cell_10 = ct_10,
          cell_01 = ct_01,
          cell_00 = ct_00
        ),
        .f = fisher_test_helper
      )
    ) %>%
    unnest(fish_result)
  
  # just preferences:
  test_skel %<>%
    rename(
      odds_ratio = estimate
    ) %>%
    relocate(
      p.value, .after = conf.high
    ) 
  
  test_skel %<>%
    left_join(
      .,
      select(dat_top, feature, feat_lab),
      by = c(var1 = "feature")
    ) %>%
    rename(var1_lab = feat_lab) %>%
    left_join(
      .,
      select(dat_top, feature, feat_lab),
      by = c(var2 = "feature")
    ) %>%
    rename(var2_lab = feat_lab) %>%
    mutate(
      var1_lab = ordered(var1_lab, levels = dat_top$feat_lab),
      var2_lab = ordered(var2_lab, levels = dat_top$feat_lab)
    ) %>%
    relocate(var1_lab, .after = var2) %>%
    relocate(var2_lab, .after = var1_lab)
  
  
  
  return(test_skel)
  
}
