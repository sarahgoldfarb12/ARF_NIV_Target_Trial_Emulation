# Sarah Goldfarb
# 02/02/2026

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
    
    # Device names that are included as NIV
    NIV_NAMES <- global_config$niv_devices 
    
    rm(config, global_config)
  } # -----------------  End loading local config
  
  { # -----------------  Creating output folders if missing -----------------
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
                         "labs")
    
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
    
    # Age limits for study, defined by outlier_thresholds file
    AGE_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "age_at_admission"]
    AGE_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "age_at_admission"]
  } # -----------------  End defining global variables and outlier thresholds
}# -------  End setup

{ # -------  Creating and merging hospital blocks
{ # -------  Creating hospital blocks
  # Select only hospitalizations that fit the below criteria
  # Only set the starting time admission (study end cutoff below as to not remove ends of hospital blocks)
  clif_hospitalization <- clif_hospitalization |>
    filter(age_at_admission >= AGE_MIN, 
           admission_dttm >= STUDY_START)|>
    compute()
  
  # Collect CLIF hopsitalization data as a data table
  hospitalization_information <- as.data.table(clif_hospitalization)
  
  # Order by patient ID and admission date/time
  setorder(hospitalization_information, patient_id, admission_dttm)
  
  # Create column taht stores the time of previous discharge
  hospitalization_information[, prev_discharge := shift(discharge_dttm, n = 1), by = patient_id]
  
  # Create column that stores the gap in hours between current admission and prior discharge
  hospitalization_information[, gap_hours := as.numeric(difftime(admission_dttm, prev_discharge, units = "hours"))]
  
  # Flag a hospitalization if it represents a new hospital block
  hospitalization_information[, new_block := is.na(prev_discharge) | gap_hours > HOSPITAL_BLOCK_TIME_HRS]
  
  # Compute cumulative block number per patient (e.g., first hospital block =1 , second = 2, etc.)
  hospitalization_information[, block_num := cumsum(replace_na(new_block, TRUE)), by = patient_id]
  
  # Determine end of hospitalization time for discharge location
  hospitalization_information[, end_ts := coalesce(discharge_dttm, admission_dttm)]
  
  
  #Remove helper columns
  hospitalization_information[, c("prev_discharge", "gap_hours", "new_block") := NULL]
  
  # Create functions to return NA (if it is NA) or the min/max of a group
  max_or_na <- function(x) if(all(is.na(x))) as.POSIXct(NA) else max(x, na.rm = TRUE)
  min_or_na <- function(x) if(all(is.na(x))) as.POSIXct(NA) else min(x, na.rm = TRUE)
  
  # Store information on grouped hospitalizations (hospital blocks)
  hospital_blocks <- hospitalization_information[, .(
    # Stores the patient age on first admission within block
    age_at_admission = min_or_na(age_at_admission),
    # Stores the start of the hospital block time
    block_start_admit = min_or_na(admission_dttm),
    # Stores the end of the hospital block time
    block_end_discharge = max_or_na(discharge_dttm),
    # Number of hospitalizations within a block
    num_hospitalizations = .N,
    # Stores the discharge location from last admission within block
    discharge_location = discharge_category[which.max(end_ts)] 
  ), by = .(patient_id, block_num)]
  
  # Create unique hospital block IDs
  hospital_blocks[, hospital_block_id := paste0(patient_id, "_", format(block_start_admit, "%Y%m%d%H%M%S"))]
  
  # Filter blocks within study period (already filtered for start and min age above)
  hospital_blocks <- hospital_blocks[block_start_admit <= STUDY_END &
                                       age_at_admission <= AGE_MAX]
  
  # Create key for each hospitalization_id to link to hospitalization block
  hospital_block_key <- hospitalization_information[hospital_blocks, 
                                                    on = .(patient_id, block_num), 
                                                    nomatch = 0][
                                                      , .(hospitalization_id, patient_id, hospital_block_id,
                                                          block_num, num_hospitalizations,
                                                          block_start_admit, block_end_discharge, discharge_location)
                                                    ]
} # -----------------  End creating hospital blocks

{ # -------  Merging hospital blocks with CLIF tables
  hospital_block_key_obj <- arrow::arrow_table(
    hospital_block_key |>
      select(patient_id,
             hospitalization_id,
             hospital_block_id,
             block_start_admit,
             block_end_discharge,
             discharge_location))
  
  clif_hospitalization <- clif_hospitalization |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id", "patient_id")) |>
    compute()
  
  clif_adt <- clif_adt |>
    select(hospitalization_id, hospital_id, hospital_type, in_dttm, out_dttm, location_name, location_category, location_type) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  
  clif_respiratory_support <- clif_respiratory_support |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  
  clif_medication_admin_continuous <- clif_medication_admin_continuous |>
    filter(admin_dttm >= as.POSIXct(STUDY_START)) |>
    filter(med_group == "vasoactives") |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  
  clif_crrt_therapy <- clif_crrt_therapy |>
    select(hospitalization_id, crrt_mode_category, recorded_dttm) |> 
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  
  clif_labs <- clif_labs |>
    filter(lab_category %in% c("ph_arterial", "ph_venous")) |>
    # ensure same datatype
    mutate(hospitalization_id = as.character(hospitalization_id)) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  
  clif_vitals <- clif_vitals |>
    filter(vital_category %in% c("spo2")) |>
    # ensure same datatype
    mutate(hospitalization_id = as.character(hospitalization_id)) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
    
} # -----------------  End merging hospital blocks with CLIF tables
} # -------  End creating and merging hospital blocks

# Keep track of cohort size
cohort_size <- data.frame(
  step = integer(),
  description = character(),
  n_encounters = integer(),
  n_patients = integer(),
  stringsAsFactors = FALSE
)

{ # -------  Grouping ED blocks
  final_cohort <-  clif_adt |>
    # Select only adt events within inclusion criteria (adult, study window, etc)
    semi_join(clif_hospitalization, by = "hospital_block_id")
  
  # Set final cohort as a data table
  setDT(final_cohort)
  
  # Precompute string operations
  final_cohort[,location_lower := tolower(location_category)]
  final_cohort[,ed_equivalent := location_lower %chin% ED_ALLOWED_LOCATIONS]
  
  # Set order within the data table
  setorder(final_cohort, hospital_block_id, in_dttm)
  
  # Identify and select only blocks that start in the ed
  ed_starting_blocks <- final_cohort[, .(starts_with_ed = first(location_lower) == "ed"), 
                                     by = hospital_block_id][starts_with_ed == TRUE, hospital_block_id]
  final_cohort <- final_cohort[hospital_block_id %chin% ed_starting_blocks]
  
  # Add the first and last hospital information, diagnosis information
  final_cohort[, `:=` (
    first_hospital_id = first(hospital_id),
    first_hospital_type = first(hospital_type),
    last_hospital_id = last(hospital_id),
    last_hospital_type = last(hospital_type),
    # Take the primary diagnosis associated with the last hospitalization 
    # E.g., if there are two hospitalizations, take the one affiliated with the second
    # Note: if multiple primary diagnosis within the last hospitalization in a block, will take first (alphanumerically) diagnosis of the last (chronilogically) hospitalization wtihin the block
    primary_diagnosis = last(primary_diagnosis), 
    diagnosis_code_format = last(diagnosis_code_format)
  ), by = hospital_block_id]
  
  # Add triage location info and row numbers
  final_cohort[, `:=`(
    # Identify the index of first traige location (e.g., ward, icu, stepdown); set to NA if dishcarge
    triage_idx = fcoalesce(match(TRUE, !ed_equivalent), NA_integer_),
    row_num = .I - .I[1] + 1
  ), by = hospital_block_id]
  
  # Identify the first triage location; set to discharge if index for it does not exist
  final_cohort[, `:=`(
    triage_location = ifelse(
      is.na(triage_idx),
      "discharge",
      location_category[triage_idx]
    ),
    triage_location_raw = ifelse(
      is.na(triage_idx),
      "discharge",
      location_name[triage_idx]
    )
  ), by = hospital_block_id]
  
  
  # Only keep rows before triage
  final_cohort <- final_cohort[is.na(triage_idx) | row_num < triage_idx]
  
  # Remove trailing non-ED locations
  final_cohort[, max_ed_out := max(out_dttm[location_lower == "ed"], na.rm = TRUE), by = hospital_block_id]
  final_cohort <- final_cohort[out_dttm <= max_ed_out]
  
  # Summarize
  final_cohort <- final_cohort[, {
    .(patient_id = patient_id[1L],
      start_ed = min(in_dttm),
      stop_ed = max(out_dttm),
      ed_hours = as.numeric(difftime(max(out_dttm), min(in_dttm), units="hours")),
      triage_location = triage_location[1L],
      triage_location_raw = triage_location_raw[1L],
      first_hospital_id = first_hospital_id[1L],
      first_hospital_type = first_hospital_type[1L],
      last_hospital_id = last_hospital_id[1L],
      last_hospital_type = last_hospital_type[1L],
      primary_diagnosis = primary_diagnosis[1L],
      diagnosis_code_format = diagnosis_code_format[1L]
    )
  }, by = hospital_block_id]
  
  # Update cohort size 
  cohort_size <- cohort_size |>
    add_row(
      step = 0,
      description = "All adult ED encounters within date range",
      n_encounters = nrow(final_cohort),
      n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
    )
} # -------  End grouping ED block