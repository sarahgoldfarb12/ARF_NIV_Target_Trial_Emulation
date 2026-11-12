# Lizi Shao 
# 11/02/2026

{ # -------  Setup
  
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
    
    # COVID start and end
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
    # only these two tables are used for the ICU/IMC ratio pipeline.
    required_tables <- c("hospitalization", "adt")
    
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
    missing_tables <- setdiff(
      paste0("clif_", required_tables, ".", file_type),
      basename(required_files)
    )
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
  
  
}# -------  End setup


{ # ------- Creating hospital blocks 
  
  # Collect CLIF hopsitalization data as a data table
  hospitalization_information <- as.data.table(
    dplyr::collect(clif_hospitalization)
  )
  
  # Collect CLIF hopsitalization data as a data table
  setorder(hospitalization_information, patient_id, admission_dttm)
  
  # Create column taht stores the time of previous discharge
  hospitalization_information[, prev_discharge := shift(discharge_dttm),
                              by = patient_id]
  
  # Create column that stores the gap in hours between current admission and prior discharge
  hospitalization_information[, gap_hours := as.numeric(
    difftime(admission_dttm, prev_discharge, units = "hours")
  )]
  
  # Create column that stores the gap in hours between current admission and prior discharge
  hospitalization_information[, new_block :=
                                is.na(prev_discharge) | gap_hours > HOSPITAL_BLOCK_TIME_HRS]
  
  # Compute cumulative block number per patient (e.g., first hospital block =1 , second = 2, etc.)
  hospitalization_information[, block_num :=
                                cumsum(tidyr::replace_na(new_block, TRUE)), by = patient_id]
  
  # Determine end of hospitalization time for discharge location
  hospitalization_information[, end_ts :=
                                dplyr::coalesce(discharge_dttm, admission_dttm)]
  
  # Remove helper columns
  hospitalization_information[, c("prev_discharge", "gap_hours", "new_block") := NULL]
  
  # preserve numeric versus POSIXct type for all-missing groups.
  max_or_na <- function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)
  min_or_na <- function(x) if (all(is.na(x))) x[NA_integer_] else min(x, na.rm = TRUE)
  
  hospital_blocks <- hospitalization_information[, .(
    age_at_admission = min_or_na(age_at_admission),
    block_start_admit = min_or_na(admission_dttm),
    block_end_discharge = max_or_na(discharge_dttm),
    num_hospitalizations = .N,
    discharge_location = discharge_category[which.max(end_ts)],
    has_missing_discharge = anyNA(discharge_dttm)
  ), by = .(patient_id, block_num)]
  
  hospital_blocks[, hospital_block_id := paste0(
    patient_id, "_", format(block_start_admit, "%Y%m%d%H%M%S")
  )]
  
  # Apply study overlap AFTER creating blocks.
  # This retains admissions before STUDY_START that extend into the window.
  # Keep uncertain discharge endpoints for patients who don't yet have a discharge dttm?
  hospital_blocks <- hospital_blocks[
    block_start_admit < STUDY_END &
      (has_missing_discharge | block_end_discharge > STUDY_START)
  ]
  
  # Create key for each hospitalization_id to link to hospitalization block
  hospital_block_key <- unique(
    hospitalization_information[
      hospital_blocks, on = .(patient_id, block_num), nomatch = 0
    ][, .(hospitalization_id, patient_id, hospital_block_id)]
  )
  
  # Everything up until this point should be identical with existing general ETL pipeline

{ # ------- ICU/IMC ADT intervals and unique patient-days
  
  # NEW: use the hospital on each ADT interval, including interhospital transfers.
  capacity_adt <- clif_adt |>
    dplyr::select(hospitalization_id, hospital_id, location_category,
                  in_dttm, out_dttm) |>
    dplyr::collect() |>
    as.data.table()
  
  capacity_adt <- merge(
    capacity_adt, hospital_block_key,
    by = "hospitalization_id", all = FALSE, sort = FALSE
  )
  
  # Remove again intervals definitively outside the study window.
  capacity_adt <- capacity_adt[
    (is.na(in_dttm) | in_dttm < STUDY_END) &
      (is.na(out_dttm) | out_dttm > STUDY_START)
  ]
  
  hospital_ids <- sort(unique(capacity_adt$hospital_id[
    !is.na(capacity_adt$hospital_id)
  ]))
  
  # Limit categories to ICU and stepdown
  capacity_adt[, location_lower := tolower(trimws(location_category))]
  capacity_adt <- capacity_adt[location_lower %chin% c("icu", "stepdown")]
  
  # For every adt record, create the following variables:
  # maximum of in_dttm of ICU or stepdown or STUDY_START (for those triaged before STUDY_START)
  capacity_adt[, stay_start := pmax(in_dttm, STUDY_START)]
  # maximum for discharge
  capacity_adt[, stay_end := pmin(out_dttm, STUDY_END)]
  # the date of the ICU or stepdown admission
  capacity_adt[, first_date := as.Date(stay_start, tz = site_time_zone)]
  # the date of the ICU or stepdown discharge
  capacity_adt[, last_date := as.Date(stay_end, tz = site_time_zone)]
  # If someone is discharged at midnight, their end date is midnight - 1. This applies to discharge only. 
  capacity_adt[, end_midnight := as.POSIXct(
    as.character(last_date), tz = site_time_zone
  )]
  capacity_adt[stay_end == end_midnight, last_date := last_date - 1L]
  
  # Collapse identical date ranges before expanding; no rows for block gaps.
  day_ranges <- unique(capacity_adt[, .(
    patient_id, hospital_id, location_lower, first_date, last_date
  )])
  
  # Expands multiple rows (# of days) per group, similar to reframe() in dplyr
  patient_days <- day_ranges[, .(
    calendar_date = seq.Date(first_date, last_date, by = "day")
  ), by = .(patient_id, hospital_id, location_lower, first_date, last_date)]
  
  # Exact counting key: patient * hospital * calendar date * ICU/IMC.
  patient_days <- unique(patient_days[, .(
    patient_id, hospital_id, calendar_date, location_lower
  )])
  
  daily_census <- patient_days[, .(n_patients = .N),
                               by = .(hospital_id, calendar_date, location_lower)]
  
  # Estimate capacity using each hospital/unit's maximum daily patient count.
  hospital_capacity <- daily_census[, .(
    capacity = max(n_patients)
  ), by = .(
    hospital_id, location_lower
  )]
  
  hospital_capacity <- hospital_capacity[order(hospital_id, location_lower)]
}

{ # ------- Hospital capacity estimates and daily ICU/IMC ratio
  
  daily_census_wide <- dcast(daily_census, hospital_id + calendar_date ~ location_lower,
                             value.var = "n_patients")
  
  hospital_capacity_wide <- dcast(hospital_capacity, hospital_id ~ location_lower,
                                  value.var = "capacity")
  
  setnames(hospital_capacity_wide, old = c("icu", "stepdown"), new = c("max_icu", "max_stepdown"))
  
  daily_capacity_ratio <- merge(daily_census_wide, hospital_capacity_wide, by = "hospital_id", all.x = T)
  
  daily_capacity_ratio[, capacity_ratio := ((icu)/(max_icu))/((stepdown)/(max_stepdown))]

}

# Save separate capacity outputs; leave the analytical cohort outputs intact.
capacity_output_dir <- file.path(project_location, paste0(site, "_project_output"),
                                 "capacity")
dir.create(capacity_output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(capacity_wide, file.path(capacity_output_dir, "hospital_capacity.csv"))
fwrite(daily_capacity_ratio, file.path(capacity_output_dir, "daily_capacity_ratio.csv"))
