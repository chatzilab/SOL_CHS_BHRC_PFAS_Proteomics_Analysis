# Table helpers ------------------------------------------------------------

collapse_terms <- function(df, id_col, term_col, n = 8, sep = "; ") {
  df %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(.data[[term_col]])) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      !!term_col := paste0(
        unique(.data[[term_col]])[seq_len(min(
          n,
          dplyr::n_distinct(.data[[term_col]])
        ))],
        collapse = sep
      ),
      .groups = "drop"
    )
}


# PWAS ---------------------------------------------------------------------

run_pfas_omics <- function(
  df,
  pfas_var,
  omics_vars,
  covariates,
  dir_output,
  dir_figures = dir_output,
  fdr_cutoff = 0.2,
  p_cutoff = 0.05,
  colors = c(
    "q-value < 0.2" = "red",
    "p-value < 0.05" = "orange",
    "Null" = "gray"
  )
) {
  res <- owas(
    df = df,
    var = pfas_var,
    omics = omics_vars,
    covars = covariates,
    var_exposure_or_outcome = "exposure",
    family = "gaussian"
  )

  res <- res %>%
    dplyr::mutate(
      estimate_beta = estimate,
      se_beta = se,
      ci_beta_low = estimate_beta - 1.96 * se_beta,
      ci_beta_high = estimate_beta + 1.96 * se_beta,
      estimate_pct = (2^estimate_beta - 1) * 100,
      ci_pct_low = (2^ci_beta_low - 1) * 100,
      ci_pct_high = (2^ci_beta_high - 1) * 100,
      pct_change_ci = sprintf(
        "%.2f (%.2f, %.2f)",
        estimate_pct,
        ci_pct_low,
        ci_pct_high
      ),
      beta_ci = sprintf(
        "%.4f (%.4f, %.4f)",
        estimate_beta,
        ci_beta_low,
        ci_beta_high
      )
    )

  sig_fdr <- res %>%
    dplyr::filter(adjusted_pval < fdr_cutoff) %>%
    dplyr::count()
  sig_p <- res %>% dplyr::filter(p_value < p_cutoff) %>% dplyr::count()

  message(glue::glue(
    "{pfas_var}: {sig_fdr$n} associations with q < {fdr_cutoff}, ",
    "{sig_p$n} with p < {p_cutoff}"
  ))

  output_file <- fs::path(dir_output, paste0(pfas_var, "_proteomics.csv"))
  readr::write_csv(res, output_file)

  res <- res %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        p_value < p_cutoff & adjusted_pval < fdr_cutoff ~ "q-value < 0.2",
        p_value < p_cutoff & adjusted_pval >= fdr_cutoff ~ "p-value < 0.05",
        TRUE ~ "Null"
      ),
      label = ifelse(Significance != "Null", as.character(feature_name), NA)
    )

  plot <- ggplot(
    data = res,
    aes(
      x = estimate_pct,
      y = -log10(p_value),
      col = Significance,
      label = label
    )
  ) +
    geom_point() +
    geom_text_repel(max.overlaps = 15) +
    scale_colour_manual(values = colors) +
    theme_minimal() +
    xlab(expression("% change per doubling of PFAS")) +
    ylab(expression("-log"[10] * "(p-value)")) +
    labs(
      title = paste0(
        "Volcano plot for ",
        pfas_var,
        " associations with proteomics"
      ),
      subtitle = glue::glue(
        "log2(proteomics) ~ beta1*log2({pfas_var}) + covariates"
      ),
      caption = paste("Covariates:", paste(covariates, collapse = ", "))
    ) +
    theme(plot.caption = element_text(hjust = 0))

  plot_file <- fs::path(dir_figures, paste0(pfas_var, "_volcano.png"))
  ggsave(filename = plot_file, plot = plot, width = 8, height = 6)

  list(
    results = res,
    plot = plot,
    output_file = output_file,
    plot_file = plot_file
  )
}


run_pfas_mixture <- function(
  df,
  exp_vars,
  omics_vars,
  covariates,
  dir_output,
  dir_figures = dir_output,
  q = 3,
  confidence_level = 0.95,
  fdr_cutoff = 0.10,
  p_cutoff = 0.05,
  colors = c(
    "FDR < 0.10" = "red",
    "p < 0.05 only" = "orange",
    "Null" = "grey70"
  )
) {
  required_vars <- unique(c(exp_vars, omics_vars, covariates))
  missing_vars <- setdiff(required_vars, names(df))

  if (length(missing_vars) > 0) {
    stop(
      "These variables are missing from df: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  res <- owas_qgcomp(
    df = df,
    expnms = exp_vars,
    omics = omics_vars,
    covars = covariates,
    q = q,
    confidence_level = confidence_level
  )

  mixture_label <- paste(exp_vars, collapse = " + ")

  res <- res %>%
    dplyr::mutate(
      mixture = mixture_label,
      n_exposures = length(exp_vars),
      estimate_beta = psi,
      ci_beta_low = lcl_psi,
      ci_beta_high = ucl_psi,
      estimate_pct = (2^estimate_beta - 1) * 100,
      ci_pct_low = (2^ci_beta_low - 1) * 100,
      ci_pct_high = (2^ci_beta_high - 1) * 100,
      beta_ci = sprintf(
        "%.3f (%.3f, %.3f)",
        estimate_beta,
        ci_beta_low,
        ci_beta_high
      ),
      pct_change_ci = sprintf(
        "%.1f%% (%.1f%%, %.1f%%)",
        estimate_pct,
        ci_pct_low,
        ci_pct_high
      ),
      Significance = dplyr::case_when(
        adjusted_pval < fdr_cutoff ~ paste0(
          "FDR < ",
          format(fdr_cutoff, nsmall = 2)
        ),
        p_value < p_cutoff ~ paste0(
          "p < ",
          format(p_cutoff, nsmall = 2),
          " only"
        ),
        TRUE ~ "Null"
      ),
      label = ifelse(
        Significance == "Null",
        NA_character_,
        as.character(feature)
      )
    )

  significant_results <- res %>%
    dplyr::filter(adjusted_pval < fdr_cutoff | p_value < p_cutoff)

  message(glue::glue(
    "{length(exp_vars)}-PFAS mixture: ",
    "{sum(res$adjusted_pval < fdr_cutoff, na.rm = TRUE)} proteins with FDR < {fdr_cutoff}, ",
    "{sum(res$p_value < p_cutoff, na.rm = TRUE)} proteins with p < {p_cutoff}"
  ))

  output_file <- fs::path(
    dir_output,
    paste0(
      "mixture_",
      paste(exp_vars, collapse = "_"),
      "_proteomics_results.csv"
    )
  )
  readr::write_csv(res, output_file)

  plot <- ggplot(
    res,
    aes(
      x = estimate_beta,
      y = -log10(p_value),
      color = Significance,
      label = label
    )
  ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey50"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey50"
    ) +
    geom_point(alpha = 0.85, size = 2) +
    geom_text_repel(
      max.overlaps = 20,
      size = 3,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey60",
      na.rm = TRUE
    ) +
    scale_color_manual(values = colors) +
    theme_minimal(base_size = 12) +
    labs(
      title = "PFAS mixture-proteomics associations",
      subtitle = glue::glue(
        "Mixture = {mixture_label}\n",
        "Effect estimate (psi) represents change in log2 protein level per one-quantile increase in the PFAS mixture"
      ),
      x = expression("Mixture effect estimate (" * psi * ")"),
      y = expression(-log[10] * "(p-value)"),
      color = NULL,
      caption = paste("Adjusted for:", paste(covariates, collapse = ", "))
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0)
    )

  plot_file <- fs::path(
    dir_figures,
    paste0("mixture_", paste(exp_vars, collapse = "_"), "_volcano.png")
  )
  ggsave(filename = plot_file, plot = plot, width = 8, height = 6, dpi = 300)

  list(
    results = res,
    significant_results = significant_results,
    plot = plot,
    output_file = output_file,
    plot_file = plot_file
  )
}


# Functional annotation -----------------------------------------------------

annotate_pwas_function <- function(
  pwas_annot,
  q_col = "adjusted_pval",
  q_cutoff = 0.2,
  uniprot_col = "UniProt",
  feature_col = "feature_name",
  olink_col = "OlinkID",
  n_go_terms = 8,
  n_reactome_terms = 8
) {
  prot_tbl <- pwas_annot %>%
    dplyr::distinct(
      .data[[feature_col]],
      .data[[uniprot_col]],
      .data[[olink_col]]
    )

  prot_tbl <- prot_tbl %>%
    dplyr::mutate(
      ENTREZID = mapIds(
        org.Hs.eg.db,
        keys = .data[[uniprot_col]],
        keytype = "UNIPROT",
        column = "ENTREZID",
        multiVals = "first"
      ),
      SYMBOL = mapIds(
        org.Hs.eg.db,
        keys = .data[[uniprot_col]],
        keytype = "UNIPROT",
        column = "SYMBOL",
        multiVals = "first"
      ),
      GENENAME = mapIds(
        org.Hs.eg.db,
        keys = .data[[uniprot_col]],
        keytype = "UNIPROT",
        column = "GENENAME",
        multiVals = "first"
      )
    )

  go_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(na.omit(prot_tbl$ENTREZID)),
    keytype = "ENTREZID",
    columns = c("GOALL", "ONTOLOGYALL")
  )

  go_bp_ids <- go_raw %>%
    dplyr::filter(!is.na(GOALL), ONTOLOGYALL == "BP") %>%
    dplyr::distinct(ENTREZID, GOALL) %>%
    dplyr::rename(GOID = GOALL)

  go_terms <- AnnotationDbi::select(
    GO.db,
    keys = unique(go_bp_ids$GOID),
    keytype = "GOID",
    columns = c("GOID", "TERM")
  ) %>%
    dplyr::distinct(GOID, TERM)

  go_bp_collapsed <- go_bp_ids %>%
    dplyr::left_join(go_terms, by = "GOID") %>%
    dplyr::filter(!is.na(TERM)) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(
      go_bp_terms = paste(
        unique(TERM)[seq_len(min(n_go_terms, dplyr::n_distinct(TERM)))],
        collapse = "; "
      ),
      .groups = "drop"
    )

  react_path_list <- AnnotationDbi::mget(
    unique(na.omit(prot_tbl$ENTREZID)),
    reactomeEXTID2PATHID,
    ifnotfound = NA
  )

  react_long <- enframe(
    react_path_list,
    name = "ENTREZID",
    value = "PATHID"
  ) %>%
    unnest_longer(PATHID, values_to = "PATHID") %>%
    dplyr::filter(!is.na(PATHID))

  react_names <- tibble(
    PATHID = unique(react_long$PATHID),
    reactome_term = unname(AnnotationDbi::mget(
      unique(react_long$PATHID),
      reactomePATHID2NAME,
      ifnotfound = NA
    ))
  ) %>%
    dplyr::filter(!is.na(reactome_term))

  react_collapsed <- react_long %>%
    dplyr::left_join(react_names, by = "PATHID") %>%
    dplyr::filter(!is.na(reactome_term)) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(
      reactome_terms = paste(
        unique(reactome_term)[seq_len(min(
          n_reactome_terms,
          dplyr::n_distinct(reactome_term)
        ))],
        collapse = "; "
      ),
      .groups = "drop"
    )

  prot_annot <- prot_tbl %>%
    dplyr::left_join(go_bp_collapsed, by = "ENTREZID") %>%
    dplyr::left_join(react_collapsed, by = "ENTREZID") %>%
    dplyr::mutate(
      .annot_text = str_to_lower(paste(
        go_bp_terms,
        reactome_terms,
        sep = " ; "
      )),
      category = dplyr::case_when(
        str_detect(
          .annot_text,
          "innate|toll|tlr|lipopolysaccharide|lps|cd14|complement|neutrophil|cytokine|interleukin|tumor necrosis"
        ) ~ "Innate immunity",
        str_detect(
          .annot_text,
          "cholesterol|lipid|lipoprotein|fatty acid|triglyceride|apolipoprotein"
        ) ~ "Lipid metabolism",
        str_detect(
          .annot_text,
          "coagulation|platelet|fibrin|hemostas|thrombo"
        ) ~ "Coagulation / hemostasis",
        str_detect(
          .annot_text,
          "oxidative|reactive oxygen|glutathione|peroxid|antioxid"
        ) ~ "Oxidative stress",
        str_detect(
          .annot_text,
          "extracellular matrix|collagen|integrin|cell adhesion"
        ) ~ "ECM / adhesion",
        TRUE ~ "Other / unclear"
      )
    ) %>%
    dplyr::select(-.annot_text)

  pwas_annot %>%
    dplyr::left_join(
      prot_annot %>%
        dplyr::select(
          tidyselect::all_of(c(feature_col, uniprot_col, olink_col)),
          ENTREZID,
          SYMBOL,
          GENENAME,
          go_bp_terms,
          reactome_terms,
          category
        ),
      by = c(
        setNames(feature_col, feature_col),
        setNames(uniprot_col, uniprot_col),
        setNames(olink_col, olink_col)
      )
    )
}


# WGCNA models --------------------------------------------------------------

run_me_model <- function(me_name, pfas_name, data, covars) {
  fml <- paste(me_name, "~", pfas_name, "+", paste(covars, collapse = " + "))
  fit <- lm(as.formula(fml), data = data)

  broom::tidy(fit, conf.int = TRUE) %>%
    dplyr::filter(term == pfas_name) %>%
    dplyr::mutate(
      module = str_remove(me_name, "^ME"),
      pfas = pfas_name
    ) %>%
    dplyr::select(
      module,
      pfas,
      estimate,
      std.error,
      conf.low,
      conf.high,
      p.value
    )
}


run_me_qgcomp <- function(me_name, data, pfas_names, covars, q = 4) {
  fml <- paste(me_name, "~", paste(c(pfas_names, covars), collapse = " + "))

  fit <- qgcomp.noboot(
    f = as.formula(fml),
    expnms = pfas_names,
    data = data,
    q = q,
    family = gaussian()
  )

  psi <- fit$psi
  se <- sqrt(fit$var.psi)

  tibble(
    module = str_remove(me_name, "^ME"),
    mixture = "PFAS mixture",
    estimate = psi,
    std.error = se,
    conf.low = psi - 1.96 * se,
    conf.high = psi + 1.96 * se,
    p.value = 2 * pnorm(abs(psi / se), lower.tail = FALSE),
    q = q
  )
}


run_me_qgcomp_weights <- function(me_name, data, pfas_names, covars, q = 4) {
  fml <- paste(me_name, "~", paste(c(pfas_names, covars), collapse = " + "))

  fit <- qgcomp.noboot(
    f = as.formula(fml),
    expnms = pfas_names,
    data = data,
    q = q,
    family = gaussian()
  )

  dplyr::bind_rows(
    tibble(
      module = str_remove(me_name, "^ME"),
      pfas = names(fit$pos.weights),
      weight = as.numeric(fit$pos.weights),
      direction = "Positive"
    ),
    tibble(
      module = str_remove(me_name, "^ME"),
      pfas = names(fit$neg.weights),
      weight = as.numeric(fit$neg.weights),
      direction = "Negative"
    )
  ) %>%
    dplyr::filter(!is.na(weight), weight > 0)
}


# Enrichment ----------------------------------------------------------------

map_proteins_to_entrez <- function(protein_vec, from_type = "SYMBOL") {
  protein_vec <- unique(na.omit(protein_vec))

  mapped <- clusterProfiler::bitr(
    protein_vec,
    fromType = from_type,
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )

  mapped %>%
    dplyr::distinct(.data[[from_type]], ENTREZID, .keep_all = TRUE)
}


run_reactome_enrichment_entrez <- function(
  df,
  universe_df,
  group_name,
  entrez_col = "ENTREZID",
  p_cutoff = 0.05,
  padj_cutoff = 0.05,
  min_genes = 3
) {
  input_genes <- df %>%
    dplyr::filter(!is.na(.data[[entrez_col]]), .data[[entrez_col]] != "") %>%
    dplyr::mutate(ENTREZID_chr = as.character(.data[[entrez_col]])) %>%
    dplyr::distinct(ENTREZID_chr) %>%
    dplyr::pull(ENTREZID_chr)

  universe_genes <- universe_df %>%
    dplyr::filter(!is.na(.data[[entrez_col]]), .data[[entrez_col]] != "") %>%
    dplyr::mutate(ENTREZID_chr = as.character(.data[[entrez_col]])) %>%
    dplyr::distinct(ENTREZID_chr) %>%
    dplyr::pull(ENTREZID_chr)

  if (length(input_genes) < min_genes) {
    return(tibble(
      enrichment_group = group_name,
      n_input = length(input_genes),
      n_universe = length(universe_genes)
    ))
  }

  reactome_obj <- ReactomePA::enrichPathway(
    gene = input_genes,
    universe = universe_genes,
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = TRUE
  )

  reactome_df <- as.data.frame(reactome_obj)

  if (nrow(reactome_df) == 0) {
    return(tibble(
      enrichment_group = group_name,
      n_input = length(input_genes),
      n_universe = length(universe_genes)
    ))
  }

  reactome_df %>%
    as_tibble() %>%
    dplyr::mutate(
      enrichment_group = group_name,
      n_input = length(input_genes),
      n_universe = length(universe_genes),
      GeneRatio_num = map_dbl(
        str_split(GeneRatio, "/"),
        ~ as.numeric(.x[1]) / as.numeric(.x[2])
      ),
      BgRatio_num = map_dbl(
        str_split(BgRatio, "/"),
        ~ as.numeric(.x[1]) / as.numeric(.x[2])
      ),
      RichFactor = GeneRatio_num / BgRatio_num,
      neg_log10_p = -log10(pvalue),
      neg_log10_fdr = -log10(p.adjust),
      sig_p = pvalue < p_cutoff,
      sig_fdr = p.adjust < padj_cutoff,
      Description_clean = Description %>%
        str_replace_all(" \\(.*?\\)", "") %>%
        str_replace_all("_", " "),
      Description_short = str_wrap(Description_clean, width = 45)
    ) %>%
    dplyr::arrange(p.adjust, pvalue)
}


make_ukb_protein_sets <- function(
  ukb,
  phenotypes,
  adjp_cutoff = 0.05,
  id_col = c("UniProt", "OlinkID"),
  direction_choice = c("any", "risk", "protective")
) {
  id_col <- match.arg(id_col)
  direction_choice <- match.arg(direction_choice)

  ukb_f <- ukb %>%
    dplyr::filter(.data$phenotype %in% phenotypes) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(P.Value), !is.na(HR)) %>%
    dplyr::mutate(
      ukb_direction = dplyr::case_when(
        HR > 1 ~ "risk",
        HR < 1 ~ "protective",
        TRUE ~ "null"
      )
    ) %>%
    dplyr::filter(P.Value <= adjp_cutoff)

  if (direction_choice != "any") {
    ukb_f <- ukb_f %>% dplyr::filter(.data$ukb_direction == direction_choice)
  }

  sets <- ukb_f %>%
    dplyr::group_by(phenotype) %>%
    dplyr::summarise(proteins = list(unique(.data[[id_col]])), .groups = "drop")

  list(
    ukb_filtered = ukb_f,
    sets_table = sets,
    set_list = rlang::set_names(sets$proteins, sets$phenotype)
  )
}


run_fgsea_by_pfas <- function(
  pfas_name,
  pwas_results,
  protein_list,
  beta_col = "estimate_beta",
  p_col = "p_value"
) {
  df <- pwas_results %>%
    dplyr::filter(
      pfas == pfas_name,
      !is.na(UniProt),
      !is.na(.data[[beta_col]]),
      !is.na(.data[[p_col]]),
      .data[[p_col]] > 0
    )

  ranks <- sign(df[[beta_col]]) * (-log10(df[[p_col]]))
  names(ranks) <- df$UniProt

  ranks <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])
  ranks <- sort(ranks, decreasing = TRUE)

  fgsea::fgseaMultilevel(
    pathways = protein_list,
    stats = ranks,
    minSize = 15,
    maxSize = 5000
  ) %>%
    dplyr::mutate(pfas = pfas_name)
}


code_from_years_edu <- function(years_edu) {
  suppressWarnings(as.integer(as.numeric(years_edu)))
}


bhrc_edu_3level <- function(years_edu_code) {
  dplyr::case_when(
    years_edu_code <= 11 ~ "Less than high school",
    years_edu_code == 12 ~ "High school graduate",
    years_edu_code >= 13 ~ "Some college and beyond",
    TRUE ~ NA_character_
  ) %>%
    factor(
      levels = c(
        "Less than high school",
        "High school graduate",
        "Some college and beyond"
      )
    )
}
