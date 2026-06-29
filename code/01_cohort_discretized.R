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
    required_tables <- c("patient", # sex
                         "hospitalization", # age at admission
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
                                   col_types = cols(patient_id = col_character(),
                                                    hospitalization_id  = col_character(),
                                                    hospital_block_id = col_character()),
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
    
    
    required_tables <- c("patient", # sex
                         "hospitalization", # age at admission
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
    
    cat("---Vitals starting...\n")
    clif_vitals <- clif_vitals |>
      filter(vital_category %in% c("temp_c","heart_rate","map","spo2", "height_cm", "weight_kg")) |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Vitals complete!\n")
    
    cat("---ADT starting...\n")
    clif_adt <- clif_adt |>
      select(hospitalization_id, hospital_id, in_dttm, out_dttm, location_name, location_category) |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---ADT complete!\n")
    
    cat("---Patient starting...\n")
    clif_patient <- clif_patient |>
      inner_join(final_cohort %>% 
                   select(patient_id), by = c("patient_id")) |>
      compute()
    cat("---Patient complete!\n")
    
    cat("---Labs starting...\n")
    clif_labs <- clif_labs |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Labs complete!\n")
    
    cat("---hospitalization starting...\n")
    clif_hospitalization <- clif_hospitalization |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---hospitalization complete!\n")
    
    cat("---hospital_diagnosis starting...\n")
    clif_hospital_diagnosis <- clif_hospital_diagnosis |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---hospital_diagnosis complete!\n")
    
    
  } # -----------------  End subsetting CLIF tables
  
  cat("End Setup!\n")
  
}# -------  End setup


{ # -----------------  Defining baseline table
  
  baseline_chars <- final_cohort |>
    select(
      patient_id, #Saved in final cohort
      hospital_block_id, #final_cohort
      niv_start, #final_cohort
      t_0 #final_cohort
    ) %>% 
    mutate(year = as.integer(year(niv_start))) %>% 
    
    #age at admission at the time of the visit
    left_join(hospital_block_key_obj %>% 
                select(-hospital_block_id) %>% 
                collect(),
              by = "patient_id") %>% 
    left_join(clif_hospitalization %>% 
                select(hospitalization_id, admission_dttm, age_at_admission) %>% 
                collect(), 
              by = "hospitalization_id") %>%
    group_by(patient_id) %>%
    arrange(admission_dttm) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(age = as.integer(age_at_admission)) %>% 
    select(-age_at_admission) %>% 
    
    #albumin value by t=0 (operationalized by collect_dttm. NOTE: Check with Chad collect_dttm OR result_dttm?)
    left_join(clif_labs %>% 
                select(hospitalization_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>% 
                filter(lab_category == "albumin") %>% 
                collect(),
              by = "hospitalization_id") %>% 
    mutate(
      albumin_baseline = if_else(
        !is.na(lab_collect_dttm) & lab_collect_dttm < t_0,
        coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
        NA_real_
      )
    ) %>% 
    group_by(patient_id) %>% 
    arrange(is.na(albumin_baseline),
            desc(lab_collect_dttm)) %>% 
    slice(1) %>% #Defensive coding to take the most recent albumin, in case multiple labs were ordered prior to t0
    ungroup() %>% 
    select(-lab_collect_dttm, -lab_category, -lab_value, -lab_value_numeric)
  
  #Create Elixhauser flags using comorbidity package
  elixhauser_flags <- clif_hospital_diagnosis |>
    collect() |>
    filter(
      poa_present == 1,
      str_detect(diagnosis_code_format, regex("icd.*10", ignore_case = TRUE))
    ) |>
    mutate(
      diagnosis_code = str_remove_all(diagnosis_code, "\\."),
      diagnosis_code = toupper(diagnosis_code)
    ) |>
    select(hospital_block_id, diagnosis_code) |>
    comorbidity::comorbidity(
      id = "hospital_block_id",
      code = "diagnosis_code",
      map = "elixhauser_icd10_quan",
      assign0 = FALSE
    )
  
  #Reference: IMC ARF project code 02 file 
  elixhauser_scores <- elixhauser_flags |>
    mutate(
      elixhauser_index = 
        chf * 7 +
        carit * 5 +
        valv * -1 +
        pcd * 4 +
        pvd * 2 +
        pmin(hypunc + hypc, 1) * 0 +
        para * 7 +
        ond * 6 +
        cpd * 3 +
        pmin(diabunc + diabc, 1) * 0 +
        hypothy * 0 +
        rf * 5 +
        ld * 11 +
        pud * 0 +
        aids * 0 +
        lymph * 9 +
        metacanc * 12 +
        solidtum * 4 +
        rheumd * 0 +
        coag * 3 +
        obes * -4 +
        wloss * 6 +
        fed * 5 +
        blane * -2 +
        dane * -2 +
        alcohol * 0 +
        drug * -7 +
        psycho * 0 +
        depre * -3,
      
      elixhauser_count = 
        chf +
        carit +
        valv +
        pcd +
        pvd +
        pmin(hypunc + hypc, 1) +
        para +
        ond +
        cpd +
        pmin(diabunc + diabc, 1) +
        hypothy +
        rf +
        ld +
        pud +
        aids +
        lymph +
        metacanc +
        solidtum +
        rheumd +
        coag +
        obes +
        wloss +
        fed +
        blane +
        dane +
        alcohol +
        drug +
        psycho +
        depre
    ) |>
    select(hospital_block_id, elixhauser_index, elixhauser_count)
  
  #Baseline conditions
  baseline_chars <- baseline_chars %>% 
    left_join(elixhauser_flags %>% 
                select(hospital_block_id,
                       pcd,
                       rf,
                       ld,
                       metacanc) %>%
                rename(pulmonary_vasc = pcd,
                        renal_failure = rf,
                       liver_disease = ld,
                       metastatic_cancer = metacanc)) %>% 
    left_join(elixhauser_scores |> 
                select(hospital_block_id, elixhauser_index, elixhauser_count),
              by = "hospital_block_id") %>% 
    mutate(elixhauser_index = coalesce(elixhauser_index, 0),
           elixhauser_count = coalesce(elixhauser_count, 0))
  
  # -----------------  Defining varying characteristics table
  
  ##Create scaffold of 2 hour time blocks
  vary_chars <- baseline_chars %>% 
    select(patient_id, hospital_block_id, t_0) %>%
    crossing(time_block = as.integer(seq(-2,22,by = 2))) %>% #Keep t-2 to t0 block to make carrying forward easier. NOTE: Do NOT include this block in analysis.  
    mutate(block_start = t_0 + hours(time_block),
           block_end = t_0 + hours(time_block + 2)) %>% 
    select(-t_0) 
  
  #Location of the patient
  #NOTE: Time complexity is currently O(nm), can restructure later to O(nlogm)
  location <- vary_chars %>% 
    left_join(clif_adt %>% 
                select(hospital_block_id, hospital_id, in_dttm, out_dttm, location_category) %>% 
                collect(), 
              by = "hospital_block_id",
              relationship = "many-to-many") %>% 
    filter(in_dttm <= block_end, #Find the time of the location corresponding to the end of the given 2-hour block. 
           is.na(out_dttm) | out_dttm >= block_end) %>% 
    group_by(hospital_block_id, time_block) %>% 
    slice_max(in_dttm, n = 1, with_ties = F) %>%  
    ungroup() %>% 
    rename(unit_location = location_category)
  
  vary_chars <- vary_chars %>% #Add back dropped rows and apply imputation (carry forward)
    left_join(location %>% select(-patient_id, -block_start, -block_end), by = c("hospital_block_id", "time_block")) %>% 
    group_by(hospital_block_id) %>% 
    arrange(time_block) %>% 
    tidyr::fill(unit_location, .direction = "down") %>% 
    ungroup() %>% 
    select(-in_dttm, -out_dttm)
  
  #Code status
  code_status <- clif_code_status %>% 
    collect() %>%
    inner_join(baseline_chars %>% select(patient_id, t_0),
              by = "patient_id") %>%
    group_by(patient_id) %>%
    arrange(start_dttm) %>% 
    mutate(end_dttm = lead(start_dttm)) %>%
    filter(start_dttm <= t_0 + hours(24),
           is.na(end_dttm) | end_dttm >= t_0) %>%
    ungroup() %>%
    select(patient_id, start_dttm, end_dttm, code_status_category, code_status_name)
  
  code_status <- vary_chars %>% 
    left_join(code_status,
              by = "patient_id",
              relationship = "many-to-many") %>% 
    filter(start_dttm <= block_end, #Find the time of the location corresponding to the end of the given 2-hour block. 
           is.na(end_dttm) | end_dttm >= block_end) %>% 
    group_by(patient_id, time_block) %>% 
    slice_max(start_dttm, n = 1, with_ties = F) %>%  
    ungroup()
  
  vary_chars <- vary_chars %>% #Add back dropped rows and apply imputation (carry forward)
    left_join(code_status %>% select(-patient_id, -block_start, -block_end, -code_status_name, -hospital_id, -unit_location), by = c("hospital_block_id", "time_block")) %>% 
    group_by(hospital_block_id) %>% 
    arrange(time_block) %>% 
    tidyr::fill(code_status_category, .direction = "down") %>% 
    ungroup() 
  
  # TODO later
  # varying_chars <- c("unit_location", # adt
  #                    "temperature", # clif_vitals
  #                    "heart_rate", # clif_vitals
  #                    "map", # clif_vitals
  #                    "resp_device", # clif_respiratory_support
  #                    "sf_ratio", # clif_vitals, clif_respiratory_support
  #                    "pf_ratio", # clif_labs, clif_respiratory_support
  #                    "co2_v_a", # clif_labs
  #                    "ph_v_a", # clif_labs
  #                    "sodium", # clif_labs
  #                    "potassium", # clif_labs
  #                    "bicarb", # clif_labs
  #                    "wbc", # clif_labs
  #                    "lactate", # clif_labs
  #                    "anti_hypertensive_drip", # clif_medication_admin_continuous
  #                    "pressors", # clif_medication_admin_continuous
  #                    "naloxone", # clif_medication_admin_continuous, clif_medication_admin_intermittent
  #                    "non_resp_sofa", # clif_labs, clif_vitals, clif_medication_admin_continuous, clif_patient_assessments
  #                    "code_status" # clif_code_status
  #                    )
  
} # -----------------  End defining baseline table and varying characteristics 

