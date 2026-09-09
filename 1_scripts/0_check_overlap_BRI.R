library(fs)
library(here)
library(dplyr)
library(readr)
library(stringr)

here::i_am("1_scripts/0_check_overlap_BRI.R")

source(fs::path("1_scripts", "!directories.R"))

#################
# BRI
#################
BRI <- readxl::read_xlsx(
  fs::path(
    dir_home |> dirname() |> dirname() |> dirname(),
    "Chatzi Projects and Analyses",
    "BRI",
    "0_Original_Data",
    "BRI_Olink_with_demographics.xlsx"
  ),
  sheet = 1
)

BRI_prote_names <- BRI$Assay |> unique()
length(BRI_prote_names) # 1472

#################
# SOLAR CHS BHRC
#################
proteomics_l <- readRDS(fs::path(
  dir_data,
  "ShARP_BHRC_Shared",
  "ShARP_BHRC_Data_Shared",
  "OLINK_proteomics_bridged_BHRC_SOLAR_CHS.rds"
))

cohort_names <- c("SOLAR", "CHS", "BHRC")
proteomics_metadata_cols <- c("LABID", "study", "Batch")

overlap_proteins <- proteomics_l[cohort_names] %>%
  purrr::map(~ setdiff(names(.x), proteomics_metadata_cols)) %>%
  purrr::reduce(intersect)

proteomics_l <- proteomics_l[c(cohort_names, "protein_annotation")]

proteomics_l[cohort_names] <- proteomics_l[cohort_names] %>%
  purrr::map(
    ~ .x %>%
      dplyr::select(
        tidyselect::any_of(proteomics_metadata_cols),
        tidyselect::all_of(overlap_proteins)
      )
  )

proteomics_l$protein_annotation <- proteomics_l$protein_annotation[
  cohort_names
] %>%
  purrr::map(~ .x %>% dplyr::filter(meta_id %in% overlap_proteins))

BHRC_prote_names <- proteomics_l$protein_annotation$BHRC$meta_Assay
length(BHRC_prote_names)

#################
# Intersect
#################
intersect(BRI_prote_names, BHRC_prote_names) |> length() # 1427
