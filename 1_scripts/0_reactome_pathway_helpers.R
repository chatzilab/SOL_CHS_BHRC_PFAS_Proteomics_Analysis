reactome_broad_pathway_pattern <- function() {
  paste(
    c(
      "Metabolism$",
      "Immune System$",
      "Signal Transduction$",
      "Disease$",
      "Gene expression \\(Transcription\\)$",
      "Transport of small molecules$",
      "Developmental Biology$",
      "Cell Cycle$"
    ),
    collapse = "|"
  )
}

get_reactome_pathways_for_proteins <- function(
  protein_df,
  broad_pathway_pattern = reactome_broad_pathway_pattern()
) {
  empty_protein_pathways <- protein_df %>%
    dplyr::select(protein, protein_label) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      entrez_id = NA_character_,
      pathway_id = NA_character_,
      pathway_name = NA_character_,
      is_broad_pathway = NA
    )

  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    return(empty_protein_pathways)
  }

  entrez_map <- tryCatch(
    map_proteins_to_entrez(protein_df$protein_label, from_type = "SYMBOL") %>%
      dplyr::transmute(
        protein_label = SYMBOL,
        entrez_id = as.character(ENTREZID)
      ),
    error = function(e) tibble::tibble()
  )

  if (nrow(entrez_map) == 0) {
    return(empty_protein_pathways)
  }

  protein_entrez_df <- protein_df %>%
    dplyr::select(protein, protein_label) %>%
    dplyr::distinct() %>%
    dplyr::left_join(entrez_map, by = "protein_label")

  input_genes <- protein_entrez_df %>%
    dplyr::filter(!is.na(entrez_id), entrez_id != "") %>%
    dplyr::distinct(entrez_id) %>%
    dplyr::pull(entrez_id)

  if (length(input_genes) == 0) {
    return(empty_protein_pathways)
  }

  reactome_obj <- tryCatch(
    ReactomePA::enrichPathway(
      gene = input_genes,
      organism = "human",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      readable = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(reactome_obj)) {
    return(empty_protein_pathways)
  }

  reactome_df <- as.data.frame(reactome_obj)

  if (nrow(reactome_df) == 0) {
    return(empty_protein_pathways)
  }

  reactome_gene_pathways <- reactome_df %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      pathway_id = ID,
      pathway_name = Description,
      entrez_id = geneID
    ) %>%
    tidyr::separate_rows(entrez_id, sep = "/") %>%
    dplyr::mutate(
      is_broad_pathway = stringr::str_detect(
        pathway_name,
        broad_pathway_pattern
      )
    ) %>%
    dplyr::distinct(entrez_id, pathway_id, pathway_name, is_broad_pathway)

  protein_entrez_df %>%
    dplyr::left_join(reactome_gene_pathways, by = "entrez_id")
}

select_reactome_pathways_to_keep <- function(
  pathway_assignments,
  n_pathways_to_show = c(5, 10)
) {
  kept_pathways <- pathway_assignments %>%
    dplyr::filter(!is.na(pathway_name), !is_broad_pathway) %>%
    dplyr::distinct(protein, pathway_id, pathway_name) %>%
    dplyr::count(pathway_id, pathway_name, name = "n_selected_proteins") %>%
    dplyr::filter(n_selected_proteins > 1) %>%
    dplyr::arrange(dplyr::desc(n_selected_proteins), pathway_name) %>%
    dplyr::slice_head(n = n_pathways_to_show[[2]])

  if (nrow(kept_pathways) < n_pathways_to_show[[1]]) {
    kept_pathways <- pathway_assignments %>%
      dplyr::filter(!is.na(pathway_name)) %>%
      dplyr::distinct(protein, pathway_id, pathway_name) %>%
      dplyr::count(pathway_id, pathway_name, name = "n_selected_proteins") %>%
      dplyr::filter(n_selected_proteins > 1) %>%
      dplyr::arrange(dplyr::desc(n_selected_proteins), pathway_name) %>%
      dplyr::slice_head(n = n_pathways_to_show[[2]])
  }

  kept_pathways
}

assign_reactome_pathways_to_proteins <- function(
  protein_df,
  pathway_assignments,
  kept_pathways
) {
  if (
    nrow(pathway_assignments) == 0 ||
      all(is.na(pathway_assignments$pathway_name)) ||
      nrow(kept_pathways) == 0
  ) {
    return(
      protein_df %>%
        dplyr::mutate(
          reactome_pathway_id = NA_character_,
          reactome_pathway = "Other/unique Reactome pathways"
        )
    )
  }

  selected_pathway_assignments <- pathway_assignments %>%
    dplyr::inner_join(kept_pathways, by = c("pathway_id", "pathway_name")) %>%
    dplyr::group_by(protein) %>%
    dplyr::arrange(
      dplyr::desc(n_selected_proteins),
      pathway_name,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      reactome_pathway_id = pathway_id,
      reactome_pathway = pathway_name
    ) %>%
    dplyr::select(protein, reactome_pathway_id, reactome_pathway)

  protein_df %>%
    dplyr::left_join(selected_pathway_assignments, by = "protein") %>%
    dplyr::mutate(
      reactome_pathway = dplyr::case_when(
        protein %in%
          pathway_assignments$protein[
            !is.na(pathway_assignments$pathway_name)
          ] &
          is.na(reactome_pathway) ~ "Other/unique Reactome pathways",
        TRUE ~ dplyr::coalesce(
          reactome_pathway,
          "Other/unique Reactome pathways"
        )
      ),
      reactome_pathway_id = dplyr::if_else(
        reactome_pathway == "Other/unique Reactome pathways",
        NA_character_,
        reactome_pathway_id
      )
    )
}

finalize_reactome_pathway_assignments <- function(
  protein_df,
  pathway_assignments,
  kept_pathways
) {
  current_kept_pathways <- kept_pathways

  repeat {
    pathway_assignments_final <- assign_reactome_pathways_to_proteins(
      protein_df,
      pathway_assignments,
      current_kept_pathways
    )

    final_pathway_counts <- pathway_assignments_final %>%
      dplyr::filter(!is.na(reactome_pathway_id)) %>%
      dplyr::distinct(protein, reactome_pathway_id, reactome_pathway) %>%
      dplyr::count(
        reactome_pathway_id,
        reactome_pathway,
        name = "n_selected_proteins"
      )

    next_kept_pathways <- current_kept_pathways %>%
      dplyr::select(pathway_id, pathway_name) %>%
      dplyr::inner_join(
        final_pathway_counts,
        by = c(
          "pathway_id" = "reactome_pathway_id",
          "pathway_name" = "reactome_pathway"
        )
      ) %>%
      dplyr::filter(n_selected_proteins > 1) %>%
      dplyr::arrange(dplyr::desc(n_selected_proteins), pathway_name)

    if (nrow(next_kept_pathways) == nrow(current_kept_pathways)) {
      return(list(
        pathways = pathway_assignments_final,
        pathways_to_keep = next_kept_pathways
      ))
    }

    if (nrow(next_kept_pathways) == 0) {
      return(list(
        pathways = assign_reactome_pathways_to_proteins(
          protein_df,
          pathway_assignments,
          next_kept_pathways
        ),
        pathways_to_keep = next_kept_pathways
      ))
    }

    current_kept_pathways <- next_kept_pathways
  }
}

annotate_reactome_pathways <- function(
  protein_df,
  n_pathways_to_show = c(5, 10),
  broad_pathway_pattern = reactome_broad_pathway_pattern()
) {
  pathway_assignments <- get_reactome_pathways_for_proteins(
    protein_df,
    broad_pathway_pattern = broad_pathway_pattern
  ) %>%
    dplyr::arrange(protein_label, pathway_name)

  pathways_to_keep <- select_reactome_pathways_to_keep(
    pathway_assignments,
    n_pathways_to_show = n_pathways_to_show
  )

  finalize_reactome_pathway_assignments(
    protein_df,
    pathway_assignments,
    pathways_to_keep
  )
}

make_reactome_pathway_colors <- function(pathway_levels) {
  distinct_pathway_palette <- c(
    "#0072B2",
    "#56B4E9",
    "#332288",
    "#6699CC",
    "#4477AA",
    "#CC79A7",
    "#AA4499",
    "#882255",
    "#CC6677",
    "#AA3377",
    "#661100",
    "#8C564B",
    "#6A51A3",
    "#4D4D4D",
    "#999999"
  )

  pathway_colors <- stats::setNames(
    rep(distinct_pathway_palette, length.out = length(pathway_levels)),
    pathway_levels
  )
  pathway_colors["Other/unique Reactome pathways"] <- "grey55"
  pathway_colors
}
