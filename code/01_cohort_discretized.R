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
    required_tables <- c("patient", # age, sex
                         "vitals", # height (for bmi), weight (for bmi), temperature, hr, map, spo2 (for sf)
                         "adt",  # unit_location
                         "hospital_diagnosis", # elixhuaser comorbidities
                         "respiratory_support", # fio2 (for sf, pf), resp_device
                         "medication_admin_continuous", # anti_hypertensive_drip, naloxone
                         "labs", # pao2 (for pf), co2, ph, albumin, sodium, potassium, bicarb, wbc, lactate, plt, tbili, cr
                         "patient_assessments", # GCS (for non resp sofa)
                         "code_status" # code_status
                         )

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
  
  { # -----------------  Loading cohort data, hosp key, global variables, outlier thresholds
    # No pandemic cohort is final cohort
    final_cohort <- read_csv(paste0(project_location,"/private_tables/no_pandemic_one_encounter_per_patient.csv"), 
                             show_col_types=FALSE)
    
    hospital_block_key <- read_csv(paste0(project_location, "/private_tables/hospital_block_key.csv"), 
                                   show_col_types=FALSE)
    
    # Add the discharge information to the final cohort
    final_cohort <- final_cohort |>
      left_join(hospital_block_key |>
                  select(hospital_block_id, 
                         block_start=block_start_admit, 
                         block_end=block_end_discharge,
                         discharge_location) |>
                  distinct(),
                by = "hospital_block_id")
    
    # Outlier thresholds
    outlier_thresholds <- read_csv(paste0(project_location, "/outlier-thresholds/project_outlier_thresholds.csv"), 
                                   show_col_types=FALSE)
    
  } # -----------------  End defining global variables and outlier thresholds
  
  { # -----------------  Subsetting CLIF tables
    
    
    required_tables <- c("patient", # age, sex
                         "vitals", # height (for bmi), weight (for bmi), temperature, hr, map, spo2 (for sf)
                         "adt",  # unit_location
                         "hospital_diagnosis", # elixhuaser comorbidities
                         "respiratory_support", # fio2 (for sf, pf), resp_device
                         "medication_admin_continuous", # anti_hypertensive_drip, naloxone
                         "labs", # pao2 (for pf), co2, ph, albumin, sodium, potassium, bicarb, wbc, lactate, plt, tbili, cr
                         "patient_assessments", # GCS (for non resp sofa)
                         "code_status" # code_status
    )
    
    cat("Subsetting CLIF tables...\n")
    
    hospital_block_key_obj <- arrow::arrow_table(
      hospital_block_key |>
        select(patient_id,
               hospitalization_id,
               hospital_block_id) |>
        filter(hospital_block_id %in% final_cohort$hospital_block_id))

    cat("---Patient starting...\n")
    clif_adt <- clif_adt |>
      select(patient_id, birth_date, sex_category) |>
      inner_join(hospital_block_key_obj, by = c("patient_id")) |>
      compute()
    cat("---Patient complete!\n") 
    
    cat("---Vitals starting...\n")
    clif_vitals <- clif_vitals |>
      filter(vital_category %in% c("temp_c","heart_rate","sbp","dbp","map","spo2", "height_kg", "weight_kg")) |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Vitals complete!\n")
    
    cat("---ADT starting...\n")
    clif_adt <- clif_adt |>
      select(hospitalization_id, in_dttm, out_dttm, location_name, location_category) |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---ADT complete!\n")
    
    
  } # -----------------  End subsetting CLIF tables
  
  cat("End Setup!\n")
  
}# -------  End setup


{ # -----------------  Defining baseline table and varying characteristics 
  baseline_chars <- final_cohort |>
    select(
      hospital_block_id, # already saved
      year, # already saved
      hospital_id=first_hospital_id, # already saved
      patient_id, # already saved
      t0, # time enrolled in trial (4 hours after starting NIV)
      age,# clif_patient
      sex, # clif_patient
      bmi, # clif_vitals
      albumin, # clif_labs
      chronic_lung, # clif_hospital_diagnosis
      pulmonary_vasc, # clif_hospital_diagnosis
      renal_failure, # clif_hospital_diagnosis
      liver_disease, # clif_hospital_diagnosis
      met_cancer, # clif_hospital_diagnosis
      elixhauser_index # clif_hospital_diagnosis
    )
  
  varying_chars <- c("unit_location", # adt
                     "temperature", # clif_vitals
                     "heart_rate", # clif_vitals
                     "map", # clif_vitals
                     "resp_device", # clif_respiratory_support
                     "sf_ratio", # clif_vitals, clif_respiratory_support
                     "pf_ratio", # clif_labs, clif_respiratory_support
                     "co2_v_a", # clif_labs
                     "ph_v_a", # clif_labs
                     "sodium", # clif_labs
                     "potassium", # clif_labs
                     "bicarb", # clif_labs
                     "wbc", # clif_labs
                     "lactate", # clif_labs
                     "anti_hypertensive_drip", # clif_medication_admin_continuous
                     "pressors", # clif_medication_admin_continuous
                     "naloxone", # clif_medication_admin_continuous, clif_medication_admin_intermittent
                     "non_resp_sofa", # clif_labs, clif_vitals, clif_medication_admin_continuous, clif_patient_assessments
                     "code_status" # clif_code_status
                     )
} # -----------------  End defining baseline table and varying characteristics 

