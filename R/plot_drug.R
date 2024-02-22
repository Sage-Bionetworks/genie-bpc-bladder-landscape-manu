
plot_drug <- function(
    dat_drug_sub,
    pal = rep(c("#EBAC23", "#B80058", "#008CF9", "#006E00", "#00BBAD", "#D163E6", "#B24502", "#FF9287", "#5954D6", "#00C6F8", "#878500", "#00A76C", "#BDBDBD"),4),
    plot_title = NULL,
    plot_subtitle = "Colors represent regimen coding",
    rect_offset = 0.4,
    return_plotly = F
) {
  if (length(unique(dat_drug_sub$record_id)) > 1) {
    cli::cli_abort(
      "Function only intended to work on one person.  More than one exists in the input data"
    )
  }
  
  dat_drug_sub %<>%
    mutate(
      end_show = if_else(
        # if the end date is the same as the start date, perturb by 0.1 to get
        # a displayed interval.
        abs(dx_drug_end_or_lastadm_int - dx_drug_start_int) < 0.5,
        dx_drug_end_or_lastadm_int + 1,
        dx_drug_end_or_lastadm_int,
        dx_drug_end_or_lastadm_int
      )
    )
  
  dat_drug_sub$end_show
  
  dat_drug_sub %<>%
    arrange(regimen_number, dx_drug_start_int) %>%
    mutate(
      agent_f = forcats::fct_inorder(agent),
      agent_f = forcats::fct_rev(agent_f),
      agent_f_pos = as.numeric(agent_f),
      regimen_number_f = factor(regimen_number),
      start_yr = dx_drug_start_int / 365.25,
      end_yr = end_show / 365.25
    )
  
  if (is.null(plot_title)) {
    
    rec <- unique(dat_drug_sub$record_id)
    cs <- unique(dat_drug_sub$ca_seq)
    coh <- unique(dat_drug_sub$cohort)
    
    plot_title <- glue(
      "{rec}, ca_seq = {cs}, {coh} cohort."
    )
  }
  
  # Allows plotly later:
  dat_drug_sub %<>%
    mutate(
      interactive_text = glue(
        "cohort: {dat_drug_sub$cohort}
        subject: {dat_drug_sub$record_id}
        ca_seq: {dat_drug_sub$ca_seq}
        regimen_number: {dat_drug_sub$regimen_number}
        drug: {dat_drug_sub$agent_f}
        start: {round(dat_drug_sub$start_yr, 3)} yrs ({dat_drug_sub$dx_drug_start_int} days) after dx.
        end: {round(dat_drug_sub$dx_drug_end_or_lastadm_int/365.25, 3)} yrs ({dx_drug_end_or_lastadm_int} days) after dx."
      )
    )
  
  gg <- ggplot(
    dat_drug_sub,
    aes(
      xmin = start_yr, 
      xmax = end_yr, 
      ymin = agent_f_pos-rect_offset,
      ymax = agent_f_pos+rect_offset,
      fill = regimen_number_f,
      text = interactive_text
    ),
  ) + 
    geom_rect(size = 2) +
    theme_bw() +
    scale_fill_manual(values = pal) + 
    scale_x_continuous(
      minor_breaks = seq(0, 100, by = (1/6))
    ) + 
    scale_y_continuous(
      breaks = 1:length(levels(dat_drug_sub$agent_f)),
      labels = levels(dat_drug_sub$agent_f)
    ) + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      y = NULL,
      x = "Years from dx"
    ) + 
    theme(
      plot.title.position = "plot"
    )
  
  if (return_plotly) {
    return(plotly::ggplotly(gg, tooltip = 'text'))
  } else {
    return(gg)
  }
}