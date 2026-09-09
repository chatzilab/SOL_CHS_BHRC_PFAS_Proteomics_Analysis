# Home directory
dir_home <- here::here() |>
  dirname() |>
  dirname() |>
  fs::path()

# Project directory
dir_project <- here::here() |>
  fs::path()

# Analysis-ready data shared across projects
dir_data <- fs::path(dir_home, "ShARP_Data_Shared")

# Project data
dir_temp_data <- fs::path(dir_project, "0_data")

# Results and reports
dir_report <- fs::path(dir_project, "2_reports")

# Figures
dir_figures <- fs::path(dir_project, "3_figures")

# BRI SLP data
dir_slp_data <- fs::path(
  dir_home |> dirname() |> dirname() |> dirname(),
  "Chatzi Projects and Analyses",
  "BRI",
  "1_Data_Cleaning",
  "0_cleaned_data"
)
