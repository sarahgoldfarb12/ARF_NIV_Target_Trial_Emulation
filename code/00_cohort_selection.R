# Sarah Goldfarb
# 02/02/2026

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

{ # -------  Creating and merging hospital blocks
{ # -------  Creating hospital blocks
  
  cat("Creating hospital blocks...\n")
  
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
  
  cat("End creating hospital blocks!\n")
} # -----------------  End creating hospital blocks

{ # -------  Merging hospital blocks with CLIF tables
  
  cat("Merging hospital blocks with CLIF tables...\n")
  
  hospital_block_key_obj <- arrow::arrow_table(
    hospital_block_key |>
      select(patient_id,
             hospitalization_id,
             hospital_block_id,
             block_start_admit,
             block_end_discharge,
             discharge_location))
  
  cat("---Hospitalization starting...\n")
  clif_hospitalization <- clif_hospitalization |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id", "patient_id")) |>
    compute()
  cat("---Hospitalization complete!\n")
  
  cat("---ADT starting...\n")
  clif_adt <- clif_adt |>
    select(hospitalization_id, hospital_id, hospital_type, in_dttm, out_dttm, location_name, location_category, location_type) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---ADT complete!\n")
  
  cat("---Resp support starting...\n")
  clif_respiratory_support <- clif_respiratory_support |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---Resp support complete!\n")
  
  cat("---Meds starting...\n")
  clif_medication_admin_continuous <- clif_medication_admin_continuous |>
    filter(admin_dttm >= as.POSIXct(STUDY_START)) |>
    filter(med_group == "vasoactives") |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---Meds complete!\n")
  
  cat("---CRRT starting...\n")
  clif_crrt_therapy <- clif_crrt_therapy |>
    select(hospitalization_id, crrt_mode_category, recorded_dttm) |> 
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---CRRT complete!\n")
  
  cat("---Labs starting...\n")
  clif_labs <- clif_labs |>
    filter(lab_category %in% c("ph_arterial", "ph_venous")) |>
    # ensure same datatype
    mutate(hospitalization_id = as.character(hospitalization_id)) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---Labs complete!\n")
  
  cat("---Vitals starting...\n")
  clif_vitals <- clif_vitals |>
    filter(vital_category %in% c("spo2")) |>
    # ensure same datatype
    mutate(hospitalization_id = as.character(hospitalization_id)) |>
    inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
    compute()
  cat("---Vitals complete!\n")
  
  cat("End merging hospital blocks with CLIF tables!\n")
  
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

{ # -------  Grouping ED blocks (also selecting only for encounters that start in ED)
  
  cat("Grouping ED blocks...\n")
  
  cat("---Formatting cohort into data table...\n")
  final_cohort <-  clif_adt |>
    # Select only adt events within inclusion criteria (adult, study window, etc)
    semi_join(clif_hospitalization, by = "hospital_block_id") |>
    collect()
  
  # Set final cohort as a data table
  setDT(final_cohort)
  
  # Precompute string operations
  final_cohort[,location_lower := tolower(location_category)]
  final_cohort[,ed_equivalent := location_lower %chin% ED_ALLOWED_LOCATIONS]
  
  # Set order within the data table
  setorder(final_cohort, hospital_block_id, in_dttm)
  cat("---Formatting cohort into data table complete!\n")
  
  # Identify and select only blocks that start in the ed
  cat("---Selecting for only encounters that start in the ED...\n")
  ed_starting_blocks <- final_cohort[, .(starts_with_ed = first(location_lower) == "ed"), 
                                     by = hospital_block_id][starts_with_ed == TRUE, hospital_block_id]
  final_cohort <- final_cohort[hospital_block_id %chin% ed_starting_blocks]
  cat("---Selecting for only encounters that start in the ED complete!\n")

  cat("---Adding hospital information...\n")
  # Add the first and last hospital information, diagnosis information
  final_cohort[, `:=` (
    first_hospital_id = first(hospital_id),
    first_hospital_type = first(hospital_type),
    last_hospital_id = last(hospital_id),
    last_hospital_type = last(hospital_type)
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
      last_hospital_type = last_hospital_type[1L]
    )
  }, by = hospital_block_id]
  
  cat("---Adding hospital and diagnosis information complete!\n")
  
  # Update cohort size 
  cohort_size <- cohort_size |>
    add_row(
      step = 0,
      description = "All adult ED encounters within date range",
      n_encounters = nrow(final_cohort),
      n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
    )
} # -------  End grouping ED blocks (also selecting only for encounters that start in ED)

{ # -------  Inclusion criteria
  
  { # -----------------  Inclusion criteria 1: NIV started in ED
    cat("Applying Inclusion criteria 1: NIV started in ED...\n")
    
    final_cohort <- arrow::arrow_table(final_cohort) |>  
      # Link hospitalizations to any respiratory support event
      inner_join(
        clif_respiratory_support |>
          select(hospital_block_id, recorded_dttm, device_category),
        by = "hospital_block_id"
      ) |>
      # Select only respiratory support events that are NIV and within a given ED block
      filter(
        tolower(device_category) %in% NIV_NAMES,
        recorded_dttm >= start_ed,
        recorded_dttm <= stop_ed
      ) |>
      distinct() |>
      collect() |>
      group_by(hospital_block_id) |>
      # Select only the first NIV within the ED
      slice_min(recorded_dttm, with_ties = FALSE) |>
      ungroup() |>
      # Rename column to set when NIV is started
      rename(niv_start = recorded_dttm)|>
      # Rename column to set device category as the starting device
      rename(first_niv_device = device_category) |>
      # Create column that defines t = 0 as 4 hours after starting NIV
      mutate(t_0 = niv_start + hours(4))
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 1,
        description = "Including: Where NIV is started in the ED",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    cat("Inclusion criteria 1: NIV started in ED complete!\n")
    
  } # -----------------  End Inclusion criteria 1: NIV started in ED
  
  { # -----------------  Inclusion criteria 2: Continuous NIV for at least 4 hours after starting
    cat("Applying Inclusion criteria 2: Continuous NIV for at least 4 hours after starting...\n")
    
    # Will pull all hospital block IDs where a "non-niv" device is charted between starting NIV and t=0
    non_niv <- arrow::arrow_table(final_cohort |> select(hospital_block_id, niv_start, t_0)) |>  
      # Link hospitalizations to any respiratory support event
      inner_join(
        clif_respiratory_support |>
          select(hospital_block_id, recorded_dttm, device_category),
        by = "hospital_block_id"
      ) |>
      # Select only respiratory support events that are NOT NIV that occur between t = -4 and t = 0
      # Assume NA carries forward from prior device
      filter(
        !is.na(device_category),
        !(tolower(device_category) %in% NIV_NAMES),
        recorded_dttm >= niv_start,
        recorded_dttm <= t_0
      ) |>
      select(hospital_block_id) |>
      distinct() |>
      collect()
    
    # Remove the non-niv hospital blocks from cohort
    final_cohort <- final_cohort |>
      anti_join(non_niv, by="hospital_block_id")
    
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 2,
        description = "Including: Continuous NIV for at least 4 hours after starting",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    # Clean space
    rm(non_niv)
    
    cat("Inclusion criteria 2: Continuous NIV for at least 4 hours after starting complete!\n")
  } # -----------------  End Inclusion criteria 2: Continuous NIV for at least 4 hours after starting
  
} # -------  End inclusion criteria

{ # -------  Exclusion criteria
  
  { # -----------------  Exclusion criteria 1: Leave ED (admitted or discharged) within 4 hours
    cat("Applying Exclusion criteria 1: Leave ED (admitted or discharged) within 4 hours of starting NIV...\n")
    
    final_cohort <- final_cohort |>
      # Selecting only for hospital blocks where they leave the ED AFTER the 4 hours is up
      filter(stop_ed > t_0)
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 3,
        description = "Excluding: Leaving ED within 4 hours of NIV start",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    cat("Exclusion criteria 1: Leave ED (admitted or discharged) within 4 hours of starting NIV complete!\n")
  } # -----------------  End Exclusion criteria 1: Leave ED (admitted or discharged) within 4 hours
  
  { # -----------------  Exclusion criteria 2: Receive IMV before t = 0 (aka 4 hours after starting NIV)
    
    cat("Applying Exclusion criteria 2: Receive IMV before t = 0...\n")
    # Pull all intubations, selecting for first occurrence
    imv_events <- arrow::arrow_table(final_cohort |> select(hospital_block_id)) |>  
      # Link hospitalizations to any IMV event
      inner_join(
        clif_respiratory_support |>
          filter(tolower(device_category) == "imv") |>
          select(hospital_block_id, first_imv_dttm = recorded_dttm),
        by = "hospital_block_id"
      ) |>
      collect() |>
      # Order within each hospitalization
      arrange(hospital_block_id, first_imv_dttm) |>
      # Take the first instance per hospitalization
      group_by(hospital_block_id) |>
      slice(1) |>
      ungroup()
    
    # Pull all trachs, selecting for first occurrence
    trach_events <- arrow::arrow_table(final_cohort |> select(hospital_block_id)) |>  
      # Link hospitalizations to any trach event
      inner_join(
        clif_respiratory_support |>
          filter(tracheostomy == 1) |>
          select(hospital_block_id, first_trach_dttm = recorded_dttm),
        by = "hospital_block_id"
      ) |>
      collect() |>
      # Order within each hospitalization
      arrange(hospital_block_id, first_trach_dttm) |>
      # Take the first instance per hospitalization
      group_by(hospital_block_id) |>
      slice(1) |>
      ungroup()
    
    # Add first imv and trach event times to data
    final_cohort <- final_cohort |>
      left_join(imv_events, by="hospital_block_id") |>
      left_join(trach_events, by="hospital_block_id")
    
    # Remove any instances where imv or trach occurs before t = 0
    final_cohort <- final_cohort |>
      filter(
        # patient is never intubated or they are first intubated after t = 0
        (is.na(first_imv_dttm) | first_imv_dttm > t_0) &
        # AND...
          # patient is never trached or they are first trached after t = 0
          (is.na(first_trach_dttm) | first_trach_dttm > t_0) 
      )
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 4,
        description = "Excluding: IMV (intubation or trach) before t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    # Clean up
    rm(imv_events, trach_events)
    
    cat("Exclusion criteria 2: Receive IMV before t = 0 complete!\n")
      
  } # -----------------  End Exclusion criteria 2: Receive IMV before t = 0 (aka 4 hours after starting NIV)
  
  { # -----------------  Exclusion criteria 3: Receive vasopressors before t = 0 (aka 4 hours after starting NIV)
    
    cat("Applying Exclusion criteria 3: Receive vasopressors before t = 0...\n")
    
    # Identify all "first" pressor events within a hospital block, no matter the time
    pressor_events <- clif_medication_admin_continuous |>
      # Only look at vasoactives 
      filter(
        # Only want vasoactive drugs
        tolower(med_group) == "vasoactives",
        # Remove cases when stopped
        tolower(mar_action_name) != "stopped",
        # Select only rows with valid (existing, numeric, and positive) values
        !is.na(med_dose) & !is.nan(med_dose) & med_dose > 0
      )|>
      # Select only relevant variables
      select(hospital_block_id, first_vasoactive_time=admin_dttm,
             first_vasoactive_name=med_category, 
             first_vasoactive_dose=med_dose, 
             first_vasoactive_unit=med_dose_unit, mar_action_name) |>
      # Only keep entries for hospitalizations we are including
      semi_join(final_cohort, by = "hospital_block_id") |>
      collect() |>
      # Order within each hospitalization
      arrange(hospital_block_id, first_vasoactive_time) |>
      # Take the first instance per hospitalization
      group_by(hospital_block_id) |>
      slice(1) |>
      ungroup()
    
    # Store the first vasoactive times
    final_cohort <- final_cohort |>
      left_join(
        pressor_events |>
          select(hospital_block_id, first_vasoactive_time, first_vasoactive_name, 
                 first_vasoactive_dose, first_vasoactive_unit),
        by = "hospital_block_id"
      )
    
    # Remove any instances where vasoactive occurs before t = 0
    final_cohort <- final_cohort |>
      filter(
        # patient is never on vasoactive or they are first getting vasoactives after t = 0
        (is.na(first_vasoactive_time) | first_vasoactive_time > t_0)
      )
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 6,
        description = "Excluding: Pressors before t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    rm(pressor_events)
    
    cat("Exclusion criteria 3: Receive vasopressors before t = 0 complete!\n")
    
  } # -----------------  End Exclusion criteria 3: Receive vasopressors before t = 0 (aka 4 hours after starting NIV)
  
  { # -----------------  Exclusion criteria 4: Receive CRRT before t = 0 (aka 4 hours after starting NIV)
    
    cat("Applying Exclusion criteria 4: Receive CRRT before t = 0...\n")
    
    # Identify all "first" pressor events within a hospital block, no matter the time
    crrt_events <- clif_crrt_therapy |>
      select(hospital_block_id, first_crrt_mode=crrt_mode_category, first_crrt_dttm=recorded_dttm) |>
      head(10)|>
      collect()
    
    # Merge with current patient data
    final_cohort <- final_cohort |>
      left_join(crrt_events, by="hospital_block_id")
    
    # Remove any instances where crrt occurs before t = 0
    final_cohort <- final_cohort |>
      filter(
        # patient never gets crrt or they first get crrt after t = 0
        (is.na(first_crrt_dttm) | first_crrt_dttm > t_0)
      )
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 7,
        description = "Excluding: CRRT before t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    rm(crrt_events)
    cat("Exclusion criteria 4: Receive CRRT before t = 0 complete!\n")
  } # -----------------  End Exclusion criteria 4: Receive CRRT before t = 0 (aka 4 hours after starting NIV)
  
  { # -----------------  Exclusion criteria 5: Most recent pH <= 7.1 at t = 0
    cat("Applying Exclusion criteria 5: Most recent pH <= 7.1...\n")
    
    severe_ph <- clif_labs |>
      select(hospital_block_id, lab_result_dttm, lab_category, lab_value) |>
      # Add t = 0 timestamps for each lab result
      left_join(arrow::arrow_table(final_cohort |> select(hospital_block_id, t_0)),
                by="hospital_block_id") |>
      # Only select ph values that result before t = 0
      filter(lab_result_dttm < t_0) |>
      collect() |>
      mutate(lab_value = as.numeric(lab_value)) |>
      # Look only at pH labs, filtering out outliers (must do after collect)
      filter((lab_category == "ph_arterial" & 
                lab_value >= outlier_thresholds$lower_limit[outlier_thresholds$variable_name=="ph_arterial"] &
                lab_value <= outlier_thresholds$upper_limit[outlier_thresholds$variable_name=="ph_arterial"])
             | (lab_category == "ph_venous"& 
                  lab_value >= outlier_thresholds$lower_limit[outlier_thresholds$variable_name=="ph_venous"] &
                  lab_value <= outlier_thresholds$upper_limit[outlier_thresholds$variable_name=="ph_venous"])) |>
      # Order within each hospitalization
      arrange(hospital_block_id, lab_result_dttm) |>
      group_by(hospital_block_id) |>
      # Take the LAST instance per hospitalization (already filtered to before t = 0)
      slice(n()) |>
      ungroup() |>
      # Identify all pH values that are "too severe" to include
      filter(lab_value <= PH_CUTOFF)
    
    # Remove blocks with pH too severe
    final_cohort <- final_cohort |>
      anti_join(severe_ph, by="hospital_block_id")
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 8,
        description = "Excluding: Most recent pH <= 7.1 at t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    cat("Exclusion criteria 5: Most recent pH <= 7.1 complete!\n")
    
    # Clean space
    rm(severe_ph)
    
  } # -----------------  End Exclusion criteria 5: Most recent pH <= 7.1 at t = 0
  
  { # -----------------  Exclusion criteria 6: Most recent SF <= 120 at t = 0
    cat("Applying Exclusion criteria 6: Most recent SF <= 120...\n")
    
    recent_sat <- clif_vitals |>
      select(hospital_block_id, vital_dttm=recorded_dttm, vital_category, vital_value) |>
      # Look only at o2 sats
      filter(vital_category == "spo2") |>
      # Add t = 0 timestamps for each measurement
      left_join(arrow::arrow_table(final_cohort |> select(hospital_block_id, t_0)),
                by="hospital_block_id") |>
      # Only select sats that are measured before t = 0
      filter(vital_dttm < t_0) |>
      collect() |>
      mutate(vital_value = as.numeric(vital_value)) |>
      filter(vital_value >= outlier_thresholds$lower_limit[outlier_thresholds$variable_name=="spo2"] &
               vital_value <= outlier_thresholds$upper_limit[outlier_thresholds$variable_name=="spo2"]) |>      
      # Order within each hospitalization
      arrange(hospital_block_id, vital_dttm) |>
      group_by(hospital_block_id) |>
      # Take the LAST instance per hospitalization (already filtered to before t = 0)
      slice(n()) |>
      ungroup() |>
      # Only select sats that can be used for SF ratio (aka, must be less than or equal to 97%)
      filter(vital_value <= 97)
      
    recent_fio2 <-  clif_respiratory_support |>
      select(hospital_block_id, recorded_dttm, device_category, fio2_set) |>
      # Only select relevant patients
      inner_join(arrow::arrow_table(final_cohort),
                by="hospital_block_id") |>
      collect()|>
      # Filter to keep only values within the thresholds (inclusive)
      # Keep NAs too, as assume they stay the same if on same device
      filter(is.na(fio2_set) | 
               (fio2_set >= outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "fio2_set"] 
                & fio2_set <= outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "fio2_set"])) |>
      ## Create ongoing fio2 ##
      # Order by time
      arrange(hospital_block_id, recorded_dttm) |>
      group_by(hospital_block_id) |>
      # Create running variable to keep track of device groupings
      mutate(device_run = rleid(hospital_block_id, device_category)) |>
      ungroup() |>
      arrange(hospital_block_id, recorded_dttm) |>
      # Group by given hospitalization and within a device run (e.g., no carrying through if NC to BiPAP)
      group_by(hospital_block_id, device_run) |>
      # Carry forward FiO2 within same device run (fill NAs with prior value)
      mutate(fio2_set = zoo::na.locf(fio2_set, na.rm = FALSE)) |>
      # Carry backward if first entries missing FiO2 (fill NAs with next value)
      mutate(fio2_set = zoo::na.locf(fio2_set, fromLast = TRUE, na.rm = FALSE)) |>
      ungroup() |>
      # Merge with last sat and require FiO2 to be set BEFORE last sat
      left_join(recent_sat |> select(hospital_block_id, vital_dttm), by="hospital_block_id") |>
      # Select only FiO2s set before the last sat
      # This also removes any FiO2s after t = 0, since last sat is by definition before t = 0
      # Also removes any NA values for vital dttm (any rows without a most recent sat)
      filter(recorded_dttm < vital_dttm) |>
      # Select only the last fio2 before the last sat at t = 0
      # Order within each hospitalization
      arrange(hospital_block_id, recorded_dttm) |>
      group_by(hospital_block_id) |>
      # Take the LAST instance per hospitalization (already filtered to before t = 0)
      slice(n()) |>
      ungroup() |>
      # Selecting only relevant columns
      select(hospital_block_id, recorded_dttm, device_category, fio2_set)
    
    # Calculating SF ratios
    severe_SF <- recent_fio2 |>
      inner_join(recent_sat, by="hospital_block_id") |>
      mutate(sf_ratio = vital_value/fio2_set) |>
      # Only want to keep those with severe SF ratios
      filter(!is.na(sf_ratio),
             sf_ratio <= SF_CUTOFF)
    
    # Merge with final cohort
    final_cohort <- final_cohort |>
      anti_join(severe_SF, by="hospital_block_id")
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 9,
        description = "Excluding: Most recent SF <= 120 at t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    
    cat("Exclusion criteria 6: Most recent SF <= 120 complete!\n")
    
    #rm(recent_sat, recent_fio2, severe_SF)
  } # -----------------  End Exclusion criteria 6: Most recent SF <= 120 at t = 0
  
  { # -----------------  Exclusion criteria 7: Code status "comfort care only" at t = 0
    
    cat("Applying exclusion criteria 7: Code status 'comfort care only'...\n")
    
    # Only keep code status entries for patients included
    clif_code_status <- clif_code_status |>
      filter(patient_id %in% final_cohort$patient_id) |>
      compute()
    
    recent_code_status <- arrow::arrow_table(final_cohort) |>
      # Select relevant columns
      select(hospital_block_id, patient_id, t_0) |>
      # Join with code status
      left_join(clif_code_status, by="patient_id") |>
      # Only keep code status entries within specified time frame
      filter(start_dttm <= t_0) |>
      collect() |>
      arrange(patient_id, start_dttm) |>
      group_by(patient_id) |>
      # Take the LAST instance per patient at t = 0
      # (assume most recent code status carries forward, even from different hospitalization)
      slice(n()) |>
      ungroup()
    
    # Add to final cohort, remove AND patients
    final_cohort <- final_cohort |>
      left_join(recent_code_status |> 
                  select(hospital_block_id, code_status_t0 = code_status_category), 
                by = "hospital_block_id") |>
      # Allow NA, which we can presume is full
      filter(is.na(code_status_t0) | tolower(code_status_t0) != "and")
      
    # Clean space
    rm(recent_code_status)
    
    # Update cohort size 
    cohort_size <- cohort_size |>
      add_row(
        step = 10,
        description = "Excluding: Code status 'comfort care only' at t = 0",
        n_encounters = nrow(final_cohort),
        n_patients = nrow(final_cohort |> select(patient_id) |> distinct())
      )
    cat("Exclusion criteria 7: Code status 'comfort care only' complete!\n")
  } # -----------------  End Exclusion criteria 7: Code status "comfort care only" at t = 0
  
} # -------  End exclusion criteria

{ # -------  Output cohort
  # Cohort size
  write_csv(cohort_size, file.path(paste0(project_location, "/", site,"_project_output/"), paste0(site,"_cohort_size.csv")))
  
  # Send to csv
  write_csv(final_cohort, file.path(paste0(project_location, "/private_tables"), "all_encounters.csv"))
  
  # Set seed for reproducibility
  set.seed(1)
  
  # Randomly select one row per patient
  final_cohort_trimmed <- final_cohort |>
    group_by(patient_id) |>
    slice_sample(n = 1) |>
    ungroup()
  
  # Write
  write_csv(final_cohort_trimmed, file.path(paste0(project_location, "/private_tables"), "one_encounter_per_patient.csv"))
  
} # -------  End output cohort

{ # -------  Output cohort: no covid sensitivity analysis
  
  no_pandemic_cohort <- final_cohort |>
    filter(start_ed < COVID_START | start_ed > COVID_END)
  
  write_csv(no_pandemic_cohort, 
            file.path(paste0(project_location, "/private_tables"), 
                      "no_pandemic_encounters.csv"))
  
  # Output cohort size for non-pandemic
  write_csv(cohort_size |>
              add_row(
                step = 11,
                description = "Sensitivity Analysis: Exclude encounters presenting to ED during COVID pandemic",
                n_encounters = nrow(no_pandemic_cohort),
                n_patients = nrow(no_pandemic_cohort |> select(patient_id) |> distinct())
              ), 
            file.path(paste0(project_location, "/", site,"_project_output/sensitivity_analysis/"), paste0(site,"_no_pandemic_cohort_size.csv")))
  
  # Set seed for reproducibility
  set.seed(1)
  
  # Randomly select one row per patient
  no_pandemic_cohort_trimmed <- no_pandemic_cohort |>
    group_by(patient_id) |>
    slice_sample(n = 1) |>
    ungroup()
  
  # Write
  write_csv(no_pandemic_cohort_trimmed, file.path(paste0(project_location, "/private_tables"), "no_pandemic_one_encounter_per_patient.csv"))
  
} # -------  End output cohort: no covid sensitivity analysis

{ # ------- Saving hospital block info
  hospital_block_key <- hospital_block_key |>
    filter(hospital_block_id  %in% final_cohort$hospital_block_id) |>
    mutate(hospitalization_id = as.character(hospitalization_id))
  
  # Save the subset of hospital blocks
  write_csv(hospital_block_key, paste0(project_location, "/private_tables/hospital_block_key_all_encounters.csv"))

  write_csv(hospital_block_key |>
              filter(hospital_block_id  %in% final_cohort_trimmed$hospital_block_id), 
            paste0(project_location, "/private_tables/hospital_block_key_one_encounter_per_patient.csv"))
  
} # ------- End saving hospital block info