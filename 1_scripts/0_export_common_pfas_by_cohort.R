options(scipen = 999)

library(fs)
library(here)
library(dplyr)
library(readr)
library(stringr)

here::i_am("1_scripts/0_export_common_pfas_by_cohort.R")

source(fs::path("1_scripts", "!directories.R"))

SOLAR_pfas_names <- c(
  "pfos_log2",
  "pfoa_log2",
  "pfna_log2",
  "pfda_log2",
  "pfhxs_log2",
  "pfhps_log2"
)
CHS_pfas_names <- c(
  "ma_pfos_log2",
  "ma_pfoa_log2",
  "ma_pfna_log2",
  "ma_pfda_log2",
  "ma_pfhxs_log2",
  "ma_pfhps_log2"
)
BHRC_pfas_names <- c(
  "PFOS_log2",
  "PFOA_log2",
  "PFNA_log2",
  "PFDA_log2",
  "PFHxS_log2",
  "PFHpS_log2"
)

common_pfas_names <- c(
  "PFOS_log2",
  "PFOA_log2",
  "PFNA_log2",
  "PFDeA_log2",
  "PFHxS_log2",
  "PFHpS_log2"
)
common_pfas_original_names <- stringr::str_remove(common_pfas_names, "_log2$")

add_original_scale_pfas <- function(df) {
  df %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::all_of(common_pfas_names),
        ~ 2^.x,
        .names = "{stringr::str_remove(.col, '_log2$')}"
      )
    ) %>%
    dplyr::select(
      ID,
      tidyselect::all_of(common_pfas_original_names),
      tidyselect::all_of(common_pfas_names)
    )
}

dir_pfas_burden_score <- fs::path(dir_temp_data, "PFAS_burden_score")
fs::dir_create(dir_pfas_burden_score)

# SOLAR ---------------------------------------------------------------------

load(fs::path(
  dir_data,
  "ShARP_SOLAR_Shared",
  "ShARP_SOLAR_Data_Shared",
  "Solar_Analysis_20251003.Rda"
))

SOLAR <- analysis_df
rm(analysis_df)

SOLAR_common_pfas <- SOLAR %>%
  dplyr::group_by(id) %>%
  dplyr::filter(visit == min(visit)) %>%
  dplyr::ungroup() %>%
  dplyr::select(ID = id, tidyselect::all_of(SOLAR_pfas_names)) %>%
  dplyr::rename_with(
    ~common_pfas_names,
    tidyselect::all_of(SOLAR_pfas_names)
  ) %>%
  add_original_scale_pfas()

readr::write_csv(
  SOLAR_common_pfas,
  fs::path(dir_pfas_burden_score, "SOLAR_PFAS_burden_score_input.csv")
)


# CHS -----------------------------------------------------------------------

CHS <- readRDS(fs::path(
  dir_data,
  "ShARP_CHS_Shared",
  "ShARP_CHS_Data_Shared",
  "CHS MetaAir MetaChem Cleaned Redcap and Exposures Outcome Data long V4.rds"
))

CHS_common_pfas <- CHS %>%
  dplyr::filter(study == "MetaAir") %>%
  dplyr::select(ID = id, tidyselect::all_of(CHS_pfas_names)) %>%
  dplyr::rename_with(
    ~common_pfas_names,
    tidyselect::all_of(CHS_pfas_names)
  ) %>%
  add_original_scale_pfas()

readr::write_csv(
  CHS_common_pfas,
  fs::path(dir_pfas_burden_score, "CHS_PFAS_burden_score_input.csv")
)


# BHRC ----------------------------------------------------------------------

BHRC <- readRDS(fs::path(
  dir_data,
  "ShARP_BHRC_Shared",
  "ShARP_BHRC_Data_Shared",
  "CCHC_usc250_long_12172025.rds"
))

BHRC_pfas <- readRDS(fs::path(
  dir_data,
  "ShARP_BHRC_Shared",
  "ShARP_BHRC_Data_Shared",
  "PFAS_CCHC.rds"
))

BHRC_common_pfas <- BHRC %>%
  dplyr::filter(visit == 1) %>%
  dplyr::select(RRID, ID = LABID, visit) %>%
  dplyr::left_join(
    BHRC_pfas %>%
      dplyr::select(RRID, tidyselect::all_of(BHRC_pfas_names)),
    by = "RRID"
  ) %>%
  dplyr::arrange(ID, visit) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  dplyr::select(ID, tidyselect::all_of(BHRC_pfas_names)) %>%
  dplyr::rename_with(
    ~common_pfas_names,
    tidyselect::all_of(BHRC_pfas_names)
  ) %>%
  add_original_scale_pfas()

readr::write_csv(
  BHRC_common_pfas,
  fs::path(dir_pfas_burden_score, "BHRC_PFAS_burden_score_input.csv")
)


tibble::tibble(
  cohort = c("SOLAR", "CHS", "BHRC"),
  n_rows = c(
    nrow(SOLAR_common_pfas),
    nrow(CHS_common_pfas),
    nrow(BHRC_common_pfas)
  ),
  output_file = fs::path(
    dir_pfas_burden_score,
    c(
      "SOLAR_PFAS_burden_score_input.csv",
      "CHS_PFAS_burden_score_input.csv",
      "BHRC_PFAS_burden_score_input.csv"
    )
  )
) %>%
  print()

message(
  "Upload these PFAS burden score input files to ",
  "https://pfasburden.shinyapps.io/app_pfas_burden/",
  " to calculate the PFAS burden score."
)
