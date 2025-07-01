impute_vivi <- function(data, sample, replicate, protein, grouping, condition, treatment_condition, reference_condition, intensity, cutoff_MAR, cutoff_MNAR){
  
  data_complete <- data %>% 
    filter( {{ condition }} %in% c(treatment_condition, reference_condition)) %>% 
    distinct({{ grouping }}, {{ condition }}, {{ replicate }}, {{ intensity }} ) %>% 
    complete( {{ grouping }}, {{ condition }},
              {{ replicate }}, fill = list(area = NA)) ###### Check, intensity changed from area
  
  result <- data %>%
    filter( {{ condition }} %in% c(treatment_condition, reference_condition)) %>%
    distinct( {{ protein }}, {{ condition }}, {{ grouping }}, {{ replicate }}) %>%
    group_by( {{ protein }}, {{ condition }},  {{ grouping }}) %>%
    mutate(n_observations = n()) %>%
    ungroup() %>%
    mutate(missingness = "complete") %>%
    mutate(missingness = ifelse(n_observations < cutoff_MAR, "MAR", missingness)) %>%
    mutate(missingness = ifelse(n_observations < cutoff_MNAR, "MNAR", missingness)) %>%
    distinct({{ grouping }}, {{ condition }}, n_observations, missingness)
  
  result <- data_complete %>% 
    left_join(result, by = c(rlang::as_name(rlang::enquo(grouping)), rlang::as_name(rlang::enquo(condition)))) %>% 
    left_join((data %>% distinct({{ protein }}, {{ grouping }})),
              by = c(rlang::as_name(rlang::enquo(grouping)))) %>% 
    mutate(n_observations = ifelse(is.na(n_observations), 0, n_observations)) %>%
    mutate(missingness = ifelse(is.na(missingness), "MNAR", missingness)) %>%
    mutate(comparison = paste(treatment_condition, "vs", reference_condition, sep = "_")) %>%
    mutate(new_sample_id = paste({{ condition }}, {{ replicate }}, sep = "_")) %>%
    group_by( {{ protein }}, {{ condition }},  {{ grouping }}) %>%
    mutate(mean = mean({{ intensity }}, na.rm = TRUE)) %>%
    mutate(sd = sd({{ intensity }}, na.rm = TRUE)) %>%
    mutate(min = min(mean)) %>%
    ungroup() %>% 
    mutate(peptide_id = paste( {{ protein }}, {{ grouping }}, sep = "_"))
  
  complete_peptides <- result %>%
    filter(missingness == "complete") %>%
    pull(peptide_id) %>%
    unique()
  
  MNAR_peptides_to_remove <- result %>%
    filter(missingness == "MNAR") %>%
    filter(!peptide_id %in% complete_peptides) %>%
    pull(peptide_id) %>%
    unique()
  
  result <- result %>%
    filter(!peptide_id %in% MNAR_peptides_to_remove) %>%
    group_by({{ grouping }}) %>%
    mutate(mean = ifelse(is.na(mean), mean({{ intensity }}, na.rm = TRUE), mean)) %>%
    mutate(sd = ifelse(is.na(sd), mean(sd, na.rm = TRUE), sd)) %>%
    mutate(min = ifelse(is.na(min), min(mean), min)) %>%
    ungroup() %>% 
    mutate(mean_calc = ifelse(missingness == "MNAR", min - 3, mean)) %>%
    group_by(new_sample_id, {{ protein }}, {{ condition }}, {{ grouping }}, {{ intensity }}, missingness, comparison) %>%
    do(imputed_intensity = suppressWarnings(rnorm(1, mean = .$mean_calc, sd = .$sd))) %>%
    mutate(imputed_intensity = ifelse(is.na({{ intensity }}), imputed_intensity, {{ intensity }})) %>% 
    filter(!(is.na({{ intensity }}) & missingness %in% c("complete", "MAR"))) %>% 
    ungroup()
  
  return(result)
  
}


# Find largest cluster of peptide differential abundances for each protein

group_peptides_by_similarity_to_median_differential_abundance_in_largest_cluster <-
  function(data,
           protein,
           peptide,
           peptide_intensity,
           condition,
           treatment_condition,
           reference_condition,
           differential_abundance_cutoff = 1) {
    
    data <- data %>% 
      group_by({{ peptide }}, {{ condition }}, {{ protein }}) %>% 
      mutate(mean_precursor_intensity = mean({{ peptide_intensity }})) %>% 
      ungroup() %>% 
      group_by({{ peptide }}, {{ protein }}) %>%
      mutate(differential_abundance = mean_precursor_intensity[{{ condition }} == treatment_condition][1] - mean_precursor_intensity[{{ condition }} == reference_condition][1]) %>% 
      ungroup()
    
    data <- data %>% 
      group_by({{ protein }}) %>%
      distinct({{ protein }}, {{ peptide }}, differential_abundance) %>% 
      filter(n() >= 2) %>% 
      mutate(bandwidth = bw.nrd0(differential_abundance)) %>%
      mutate(n = nmodes(differential_abundance, bw = bandwidth[1])) %>%
      ungroup() %>%
      mutate(n = ifelse(n == 1, 2, n)) %>% 
      mutate(n = ifelse(n > 4, 4, n)) %>% 
      group_by({{ protein }}, n) %>% 
      do(clusters = hclust(dist(.$differential_abundance))) %>% 
      mutate(cluster_assignment = list(cutree(clusters, n))) %>% 
      ungroup() %>% 
      left_join((data %>% distinct({{ peptide }}, {{ protein }})), by = names(dplyr::select(data, {{ protein }}))) %>% 
      group_by({{ protein }}) %>% 
      mutate(row = row_number()) %>% 
      mutate(cluster = cluster_assignment[[1]][row]) %>% 
      ungroup() %>% 
      group_by({{ protein }}, cluster) %>% 
      mutate(peptides_in_cluster = n()) %>% 
      ungroup() %>% 
      group_by({{ protein }}) %>% 
      mutate(max_peptides_in_cluster = max(peptides_in_cluster)) %>% 
      ungroup() %>% 
      group_by({{ protein }}, peptides_in_cluster) %>% 
      mutate(largest_cluster_selected_randomly = ifelse((peptides_in_cluster == max_peptides_in_cluster) & (max(cluster) != min(cluster)), 
                                                        TRUE, 
                                                        FALSE)) %>% 
      mutate(max_peptides_in_cluster = ifelse((peptides_in_cluster == max_peptides_in_cluster) & (max(cluster) != min(cluster)), 
                                              max_peptides_in_cluster + cluster - 1, 
                                              max_peptides_in_cluster)) %>% 
      ungroup() %>% 
      group_by({{ protein }}) %>% 
      mutate(peptide_in_largest_cluster = ifelse(peptides_in_cluster == max_peptides_in_cluster, TRUE, FALSE)) %>% 
      ungroup() %>% 
      distinct({{ protein }}, {{ peptide }}, peptide_in_largest_cluster, largest_cluster_selected_randomly) %>% 
      left_join((data), by = c(names(dplyr::select(data, {{ protein }})), names(dplyr::select(data, {{ peptide }}))))
    
    data <- data %>% 
      filter(peptide_in_largest_cluster == TRUE) %>% 
      group_by({{ protein }}) %>% 
      mutate(median_differential_abundance_in_largest_cluster = median(differential_abundance)) %>% 
      ungroup() %>% 
      distinct({{ protein }}, median_differential_abundance_in_largest_cluster) %>% 
      left_join(data, by = names(dplyr::select(data, {{ protein }}))) %>% 
      mutate(predicted_type_0_peptide = ifelse(differential_abundance > (median_differential_abundance_in_largest_cluster - differential_abundance_cutoff) & differential_abundance < (median_differential_abundance_in_largest_cluster + differential_abundance_cutoff), TRUE, FALSE))
    
    return(data)
    
  }

# Do Grubbs test on the differential abundance of each peptide not in the largest cluster


calculate_peptide_profile_significance <- function(data, 
                                                   protein, 
                                                   peptide, 
                                                   peptide_intensity,
                                                   pep_corr_to_prot,
                                                   condition, 
                                                   treatment_condition, 
                                                   reference_condition,
                                                   structural_group){
  
  data <- data %>% 
    group_by({{ peptide }}, {{ condition }}, {{ protein }}) %>% 
    mutate(mean_precursor_intensity = mean({{ peptide_intensity }})) %>% 
    ungroup() %>% 
    group_by({{ peptide }}, {{ protein }}) %>%
    mutate(differential_abundance = mean_precursor_intensity[{{ condition }} == treatment_condition][1] - mean_precursor_intensity[{{ condition }} == reference_condition][1]) %>% 
    ungroup()
  
  data <- data %>%
    dplyr::select({{ protein }}, {{ structural_group }}) %>%
    dplyr::distinct() %>%
    dplyr::group_by({{ protein }}) %>%
    dplyr::mutate(current_peptide = {{ structural_group }}) %>%
    complete({{ structural_group }}, current_peptide) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = as.numeric({{ structural_group }} == current_peptide)) %>% 
    distinct( {{ protein }}, current_peptide, {{ structural_group }}, group) %>% 
    left_join(data %>% select(c({{ protein }} ,{{ structural_group }}, {{ pep_corr_to_prot }}, differential_abundance)) %>% distinct(), by = c("current_peptide" = names(dplyr::select(data, {{ structural_group }})), names(dplyr::select(data, {{ protein }})))) %>%
    dplyr::group_by({{ protein }}, {{ structural_group }}) %>%
    mutate(keep = ifelse(((group == 0)*({{ pep_corr_to_prot }} == 0)), 0, 1)) %>%
    mutate(keep = ifelse(((group == 1)*({{ pep_corr_to_prot }} == 1)), 0, keep)) %>%
    filter(keep == 1) %>%
    mutate(keep = max(group == 1)) %>%
    filter(keep == 1) %>%
    select(-c(keep)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by({{ structural_group }}) %>%
    filter(n() >= 4) %>%
    ungroup() %>%
    dplyr::group_by({{ protein }}, {{ structural_group }}) %>%
    do(model = grubbs.test(.$differential_abundance)) %>%
    mutate(tested_value = str_extract(model$alternative, "[\\-]?[0-9]+\\.[0-9]+")) %>%
    mutate(tested_value = as.numeric(tested_value)) %>%
    mutate(peptide_group_pval_init = model$p.value) %>%
    dplyr::ungroup() %>%
    left_join(data %>% select(c({{ protein }} ,{{ structural_group }}, differential_abundance)) %>% distinct(), by = c(names(dplyr::select(data, {{ structural_group }})), names(dplyr::select(data, {{ protein }})))) %>%
    mutate(peptide_group_pval = ifelse((abs(tested_value - differential_abundance) < 1e-5), peptide_group_pval_init, 1)) %>% 
    dplyr::group_by({{ protein }}) %>%
    dplyr::mutate(peptide_group_per_protein_adj_pval = p.adjust(peptide_group_pval, method = "BH")) %>%
    dplyr::ungroup() %>%
    right_join(data, by = c(names(dplyr::select(data, {{ structural_group }})), names(dplyr::select(data, {{ protein }})), "differential_abundance")) %>%
    arrange(peptide_group_per_protein_adj_pval)
  
  data <- data %>%
    distinct() %>%
    mutate(significant = (peptide_group_pval < 0.05 ) * ({{ pep_corr_to_prot }} == 0))
  
  return(data)
}


# Adapted from protti, to visualize the final results

peptide_profile_plot_2 <- function(data,
                                   sample,
                                   peptide,
                                   intensity_log2,
                                   grouping,
                                   targets,
                                   pep_corr_to_prot,
                                   protein_abundance_plot = FALSE,
                                   interactive = FALSE,
                                   export = FALSE,
                                   export_name = "peptide_profile_plots") {
  . <- NULL
  n_samples <- length(unique(dplyr::pull(data, {{ sample }})))
  protti_colours <- "placeholder" # assign a placeholder to prevent a missing global variable warning
  utils::data("protti_colours", envir = environment()) # then overwrite it with real data
  if (missing(targets)) stop("Please provide at least one target to plot!")
  if (!("all" %in% targets)) {
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity_log2 }}, {{ grouping }}, {{ pep_corr_to_prot }}) %>%
      tidyr::drop_na({{ intensity_log2 }}) %>%
      dplyr::filter({{ grouping }} %in% targets) %>%
      split(dplyr::pull(., !!ensym(grouping)))
  }
  if ("all" %in% targets) {
    groups <- length(unique(dplyr::pull(data, {{ grouping }})))
    message("Splitting into ", groups, " groups and returning ", groups, " plots.")
    input <- data %>%
      dplyr::distinct({{ sample }}, {{ peptide }}, {{ intensity_log2 }}, {{ grouping }}) %>%
      tidyr::drop_na({{ intensity_log2 }}) %>%
      split(dplyr::pull(., !!ensym(grouping)))
  }
  pb <- progress::progress_bar$new(
    total = length(input),
    format = " Creating plots [:bar] :current/:total (:percent) :eta"
  )
  plot <- purrr::map2(
    .x = input,
    .y = names(input),
    .f = ~ {
      pb$tick()
      ggplot2::ggplot(.x, ggplot2::aes({{ sample }}, {{ intensity_log2 }}, group = {{ peptide }}, color = {{ pep_corr_to_prot }})) +
        ggplot2::geom_point() +
        ggplot2::geom_line(size = 1) +
        ggplot2::labs(
          title = paste("Peptide profiles:", .y),
          x = "Sample",
          y = "Intensity [log2]",
          col = NULL
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.title.x = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(size = 15),
          axis.text.x = ggplot2::element_text(size = 12, angle = 75, hjust = 1),
          axis.title.y = ggplot2::element_text(size = 15),
          legend.title = ggplot2::element_text(size = 15),
          legend.text = ggplot2::element_text(size = 15),
          strip.text.x = ggplot2::element_text(size = 15),
          strip.text = ggplot2::element_text(size = 15),
          strip.background = ggplot2::element_blank()
        ) +
        ggplot2::scale_color_manual(values =  colorset)
    }
  )
  if (interactive == FALSE) {
    if (export == TRUE) {
      grDevices::pdf(
        file = paste0(export_name, ".pdf"),
        width = 10 * ceiling(n_samples / 10), # more samples need wider plot
        height = 6
      )
      suppressWarnings(print(plot))
      grDevices::dev.off()
    } else {
      return(plot)
    }
  }
  if (interactive == TRUE) {
    if (length(targets) > 1) {
      stop(strwrap("Please only provide one target for interactive plots!",
                   prefix = "\n", initial = ""
      ))
    }
    if (export == TRUE) {
      stop(strwrap("Interactive plots cannot be exported! Please decide if you
want an interactive plot or if you want to export your plots.",
                   prefix = "\n", initial = ""
      ))
    }
    plotly::ggplotly(plot[[1]], tooltip = c("x", "y", "group"))
  }
}


barcode_plot_new <- function(data,
                             start_position,
                             end_position,
                             protein_length,
                             coverage = NULL,
                             colouring = NULL,
                             protein_id = NULL,
                             facet = NULL,
                             cutoffs = NULL) {
  # Check if there is more than one protein even though protein_id was specified.
  if (!missing(protein_id)) {
    if (length(unique(dplyr::pull(data, {{ protein_id }}))) > 1) {
      stop("If data contains information of multiple proteins use the facet argument, not the protein_id argument")
    }
  }
  # Check if there are more than 20 proteins for faceting.
  if (!missing(facet)) {
    if (length(unique(dplyr::pull(data, {{ facet }}))) > 20) {
      n_proteins <- length(unique(dplyr::pull(data, {{ facet }})))
      twenty_proteins <- unique(dplyr::pull(data, {{ facet }}))[1:20]
      data <- data %>%
        dplyr::filter({{ facet }} %in% twenty_proteins)
      warning(paste(
        "Only the first 20 proteins from", rlang::as_name(enquo(facet)),
        "have been used for plotting since there are", n_proteins,
        "proteins. Consider mapping over subsetted datasets."
      ))
    }
  }
  # Apply fold change  and significance cutoff if fold change is provided
  if (!missing(cutoffs)) {
    fc_name <- names(cutoffs)[1]
    sig_name <- names(cutoffs)[2]
    fc <- cutoffs[1]
    sig <- cutoffs[2]
    
    colouring <- sym("change")
    
    data <- data %>%
      dplyr::mutate({{ colouring }} := ifelse(((!!ensym(fc_name) >= fc | !!ensym(fc_name) <= -fc) & !!ensym(sig_name) <= sig), "Changed", "Unchanged")) %>%
      dplyr::mutate({{ colouring }} := forcats::fct_rev(ifelse(is.na({{ colouring }}), "Unchanged", {{ colouring }}))) %>%
      dplyr::arrange({{ colouring }})
  }
  # Add coverage to protein ID name if present.
  if (!missing(coverage) & !missing(facet)) {
    data <- data %>%
      mutate({{ facet }} := paste0({{ facet }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  if (!missing(coverage) & !missing(protein_id)) {
    data <- data %>%
      mutate({{ protein_id }} := paste0({{ protein_id }}, " (", round({{ coverage }}, digits = 1), "%)"))
  }
  # Create plot
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(
      ymin = -2.5,
      ymax = 2.5,
      xmax = {{ end_position }} / {{ protein_length }} * 100,
      xmin = ({{ start_position }} - 1) / {{ protein_length }} * 100,
      fill = {{ colouring }}
    ),
    size = 0.7
    ) +
    ggplot2::scale_fill_manual(values = c(
      "#5680C1", "#d9363c", "#B96DAD", "#64CACA", "#81ABE9", "#F6B8D1", "#99F1E4", "#9AD1FF", "#548BDF", "#A55098", "#3EB6B6",
      "#87AEE8", "#CA91C1", "#A4E0E0", "#1D4F9A", "#D7ACD2", "#49C1C1"
    )) +
    ggplot2::scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = "Protein Sequence", title = {
      if (!missing(protein_id)) unique(dplyr::pull(data, {{ protein_id }}))
    }) +
    {
      if (!missing(facet)) ggplot2::facet_wrap(rlang::new_formula(NULL, rlang::enquo(facet)))
    } +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.ticks.x = element_blank(),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA)
    )
}

