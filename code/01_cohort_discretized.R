# Sarah Goldfarb
# 02/04/2026

{ # -------  Setup
  
  cat("Setup...\n")
  
  { # -----------------  Loading the required packages -----------------
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
                  "rprojroot")
    
    install_if_missing <- function(package) {
      if (!require(package, character.only = TRUE)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = TRUE)
      }
    }
    
    sapply(packages, install_if_missing)
    rm(packages, install_if_missing)
    
  } # -----------------  End loading the required packages
  
  { # -----------------  Loading local config -----------------
    
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
    
    # Cutoff for describing hospital as IMC capable
    IMC_CAPABLE_CUTOFF <- global_config$imc_capable_cutoff
    
    # Study start and end
    STUDY_START <- as.POSIXct(global_config$study_start, tz=site_time_zone) # Start time
    STUDY_END <- as.POSIXct(global_config$study_end, tz=site_time_zone) # End time
    
    
    COVID_START <- as.POSIXct(global_config$covid_start, tz=site_time_zone) # Start time
    COVID_END <- as.POSIXct(global_config$covid_end, tz=site_time_zone) # End time
    
    # Device names that are included as NIV
    NIV_NAMES <- global_config$niv_devices 
    
    rm(config, global_config)
  } # -----------------  End loading local config
  
  { # -----------------  Creating output folders if missing -----------------
    if (!dir.exists(paste0(project_location, "/", site, "_project_output"))) {
      dir.create(paste0(project_location,"/", site, "_project_output"))
    }
    if (!dir.exists(paste0(project_location, "/", site, "_project_output/sensitivity_analysis"))) {
      dir.create(paste0(project_location,"/", site, "_project_output/sensitivity_analysis"))
    }
    if (!dir.exists(paste0(project_location, "/private_tables"))) {
      dir.create(paste0(project_location, "/private_tables"))
    }
  } # -----------------  End creating output folders if missing
  
  { # -----------------  Loading required CLIF tables -----------------
    
    # Tables that should be set to TRUE for this project
    required_tables <- c("hospitalization", 
                         "adt", 
                         "respiratory_support", 
                         "medication_admin_continuous", 
                         "crrt_therapy",
                         "vitals",
                         "labs",
                         "code_status")
    
    # List all CLIF files in the directory
    clif_table_filenames <- list.files(path = tables_location, 
                                       pattern = paste0("^clif_.*\\.", file_type, "$"), 
                                       full.names = TRUE)
    
    # Create a lookup table for required files based on table_flags
    required_files <- clif_table_filenames[
      # Remove all file name components before and including clif
      sub(".*clif_", "", 
          # Remove file name
          sub(paste0("\\.", file_type, "$"),"",clif_table_filenames)) 
      %in% required_tables
    ]
    
    # Check if all required files are present
    missing_tables <- setdiff(required_files, clif_table_filenames)
    if (length(missing_tables) > 0) {
      stop(paste("Error: Missing required tables:", paste(missing_tables, collapse = ", ")))
    }
    
    # Define the cast function to convert large_string to string
    cast_large_utf8_to_utf8 <- function(x) {
      # Only applies to Arrow objects that have a schema()
      # (open_dataset() returns an Arrow Dataset / query)
      sch <- arrow::schema(x)
      
      large_cols <- vapply(
        sch$fields,
        function(f) f$type$ToString() %in% c("large_string", "large_utf8"),
        logical(1)
      )
      
      cols_to_cast <- sch$names[large_cols]
      if (length(cols_to_cast) == 0) return(x)
      
      x %>% mutate(across(all_of(cols_to_cast), ~ arrow::cast(.x, arrow::utf8())))
    }
    
    # Read the required files into a list of data frames
    if (file_type == "parquet") {
      data_list <- lapply(required_files, open_dataset)
      # Apply cast function to normalize Arrow string types right after import
      data_list <- lapply(data_list, cast_large_utf8_to_utf8)
    } else if (file_type == "csv") {
      data_list <- lapply(required_files, read_csv)
    } else if (file_type == "fst") {
      data_list <- lapply(required_files, read.fst)
    } else {
      stop("Unsupported file format")
    }
    
    # Assign the data frames to variables based on their file names
    for (i in seq_along(required_files)) {
      # Extract the base name of the file (without extension)
      object_name <- str_remove(basename(required_files[i]), paste0("\\.", file_type, "$"))
      # Make the object name valid for R (replace invalid characters with underscores)
      object_name <- make.names(object_name)
      # Assign the tibble to a variable with the name of the file
      assign(object_name, data_list[[i]])
    }
    
    # Clean space
    rm(data_list, i, missing_tables, 
       object_name, required_files, required_tables, cast_large_utf8_to_utf8)
    
  } # -----------------  End loading required CLIF tables
  
  { # -----------------  Defining global variables and outlier thresholds
    # Outlier thresholds
    outlier_thresholds <- read_csv(paste0(project_location, "/outlier-thresholds/project_outlier_thresholds.csv"), show_col_types=FALSE)
    
    ED_ALLOWED_LOCATIONS <- c("ed", "procedural", "other") # Lists locations that would be considered still part of ED time
    
    HOSPITAL_BLOCK_TIME_HRS <- 6 # Maximum number of hours allowed between hospitalizations to make them blocked together
    
    
    # pH limit to include
    PH_CUTOFF <- 7.1
    
    # SF limit to include
    SF_CUTOFF <- 120
    
    
    # Age limits for study, defined by outlier_thresholds file
    AGE_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "age_at_admission"]
    AGE_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "age_at_admission"]
  } # -----------------  End defining global variables and outlier thresholds
  
  cat("End Setup!\n")
  
}# -------  End setup