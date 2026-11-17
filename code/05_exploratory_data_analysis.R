# -----------------  Loading the required packages -----------------

packages <- c("duckdb", 
              "lubridate", 
              "tidyverse", 
              "dplyr",
              "table1", 
              "broom", 
              "arrow", 
              "rvest", 
              "readr", 
              "fst", 
              "data.table", 
              "collapse", 
              "tictoc",
              "yaml",
              "rprojroot",
              "comorbidity")

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE) 
    }
}

sapply(packages, install_if_missing)
rm(packages, install_if_missing)

# -----------------  End loading the required packages



# -----------------  Loading local config -----------------

# Find project root
project_root <- find_root(rprojroot::has_dir("config"))

# Read YAML config
config <- yaml::read_yaml(file.path(project_root, "config", "config.yaml"))
global_config <- yaml::read_yaml(file.path(project_root, "config", "global_config.yaml"))

# Assign config values to R variables
tables_location <- config$tables_location
project_location <- config$project_location
site <- config$institution
site_time_zone <- config$time_zone
file_type <- config$file_type

# Study start and end
STUDY_START <- as.POSIXct(global_config$study_start, tz=site_time_zone) # Start time
STUDY_END <- as.POSIXct(global_config$study_end, tz=site_time_zone) # End time


COVID_START <- as.POSIXct(global_config$covid_start, tz=site_time_zone) # Start time
COVID_END <- as.POSIXct(global_config$covid_end, tz=site_time_zone) # End time

# Device names that are included as NIV
NIV_NAMES <- global_config$niv_devices 

rm(config, global_config)

# -----------------  End loading local config -----------------



# -----------------  Loading required analysis ready tables -----------------

# Names of the three analytic tables expected to be loaded
required_tables <- c(
  "baseline_chars",
  "vary_chars",
  "outcomes_chars"
)

# Folder where the saved tables are located
tables_location <- file.path(project_location, "private_tables")

# File type to load
file_type <- "csv"

# List all files with the selected file type in the target folder
table_files <- list.files(
  path = tables_location,
  pattern = paste0("^.*\\.", file_type, "$"),
  full.names = TRUE
)

# Extract table names from file names by removing path and file extension
available_tables <- table_files %>%
  basename() %>%
  str_remove(paste0("\\.", file_type, "$"))

# Identify any required tables that are missing from the folder
missing_tables <- setdiff(required_tables, available_tables)

# Stop if any required tables are missing
if (length(missing_tables) > 0) {
  stop(
    "Error: Missing required tables: ",
    paste(missing_tables, collapse = ", ")
  )
}

# Keep only files corresponding to the required tables
required_files <- table_files[available_tables %in% required_tables]

# Function to load one table depending on the file type
load_one_table <- function(file) {
  if (file_type == "csv") {
    read_csv(file, show_col_types = FALSE)
  } else if (file_type == "parquet") {
    arrow::open_dataset(file)
  } else if (file_type == "fst") {
    fst::read_fst(file)
  } else {
    stop("Unsupported file format: ", file_type)
  }
}

# Load each required table and assign it to an object with the same name
for (file in required_files) {
  object_name <- basename(file) %>%
    str_remove(paste0("\\.", file_type, "$")) %>%
    make.names()
  
  assign(object_name, load_one_table(file), envir = .GlobalEnv)
}

rm(
  required_tables,
  table_files,
  available_tables,
  missing_tables,
  required_files,
  load_one_table,
  file,
  object_name
)

# -----------------  End loading required analysis ready tables -----------------



# -----------------  Start defining global variables and outlier thresholds -----------------

# Outlier thresholds
outlier_thresholds <- read_csv(paste0(project_location, "/outlier-thresholds/project_outlier_thresholds.csv"), 
                               show_col_types=FALSE)

# -----------------  End defining global variables and outlier thresholds -----------------



# -----------------  Processing data for persona period analysis -----------------

###Estimating 
risk_df <- vary_chars



