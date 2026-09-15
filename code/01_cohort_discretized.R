# Sarah Goldfarb & Lizi Shao
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
      inner_join(hospital_block_key_obj, by = c("patient_id", "hospitalization_id")) |>
      compute()
    cat("---hospitalization complete!\n")
    
    cat("---hospital_diagnosis starting...\n")
    clif_hospital_diagnosis <- clif_hospital_diagnosis |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---hospital_diagnosis complete!\n")
    
    cat("---Respiratory device starting...\n")
    clif_respiratory_support <- clif_respiratory_support |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Respiratory device complete!\n")
    
    cat("---Patient assessment starting...\n")
    clif_patient_assessments <- clif_patient_assessments |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Patient assessment complete!\n")
    
    cat("---Continuous medications starting...\n")
    clif_medication_admin_continuous <- clif_medication_admin_continuous |>
      inner_join(hospital_block_key_obj, by = c("hospitalization_id")) |>
      compute()
    cat("---Continuous medications complete!\n")
    
  } # -----------------  End subsetting CLIF tables
  
  cat("End Setup!\n")
  
}# -------  End setup


# -----------------  Defining baseline table

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
              collect() %>% 
              distinct(),
            by = "patient_id") %>% 
  left_join(clif_hospitalization %>% 
              select(hospitalization_id, admission_dttm, age_at_admission) %>% 
              collect() %>% 
              distinct(), 
            by = "hospitalization_id") %>%
  group_by(patient_id) %>%
  arrange(admission_dttm) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(age = as.integer(age_at_admission)) %>% 
  select(-age_at_admission) %>% 
  
  #albumin value by t=0 (operationalized by collect_dttm. NOTE: Check with Chad collect_dttm OR result_dttm?)
  left_join(clif_labs %>% 
              collect() %>% 
              distinct() %>%
              select(hospitalization_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>% 
              filter(str_to_lower(trimws(lab_category)) == "albumin"),
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
  select(-lab_collect_dttm, -lab_category, -lab_value, -lab_value_numeric, -hospitalization_id)

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

##Adding elixhauser index and count back into Baseline conditions
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
  crossing(time_block = as.integer(seq(-2,22,by = 2))) %>% #Allow 
  mutate(block_start = t_0 + hours(time_block),
         block_end = t_0 + hours(time_block + 2)) 

##Join with hospitalization id
vary_chars <- vary_chars %>% 
  left_join(clif_hospitalization %>% 
              select(patient_id, hospitalization_id, discharge_dttm, admission_dttm) %>% 
              collect() %>% 
              mutate(#
                hospitalization_id = as.character(hospitalization_id),
                #Some patients might not have a discharge time (no record of discharge)
                discharge_dttm_join = coalesce(discharge_dttm, 
                                               as.POSIXct("9999-12-31 23:59:59", tz = "UTC"))) %>% 
              distinct(), 
            by = join_by(patient_id, 
                         block_start < discharge_dttm_join,
                         block_end > admission_dttm)) %>%
  group_by(patient_id, time_block) %>%
  #If a time block traverses multiple hospitalizations, select the one closest to the end of the time block. 
  arrange(desc(admission_dttm), .by_group = TRUE) %>% 
  slice(1) %>%
  ungroup() %>% 
  select(-discharge_dttm, -discharge_dttm_join, -admission_dttm) %>% 
  #Some patients leave the ER before t+24, making their hospitalization_id NA for the corresponding time blocks. 
  group_by(patient_id) %>% 
  arrange(time_block, .by_group = TRUE) %>% 
  tidyr::fill(hospitalization_id, .direction = "down") %>% 
  ungroup() 


# ----- Transforming variables into varying characteristics table

##Location of the patient
location <- vary_chars %>% 
  select(hospital_block_id,
         time_block, 
         block_end) %>%
  left_join(clif_adt %>% 
              collect() %>% 
              select(hospital_block_id, hospital_id, in_dttm, out_dttm, location_category) %>% 
              filter(!str_to_lower(trimws(location_category)) %in% c("procedural", "radiology", "dialysis", "other")) %>% 
              distinct(), 
            by = "hospital_block_id",
            relationship = "many-to-many") %>% #Each unique patient has >1 rows, each corresponding to a time block -> specify many-to-many 
  group_by(hospital_block_id, time_block) %>% 
  filter(in_dttm <= block_end, #Find the time of the location corresponding to the end of each 2-hour block. 
         is.na(out_dttm) | out_dttm > block_end) %>% 
  arrange(desc(in_dttm), .by_group = TRUE) %>%
  slice(1) %>% #In case of duplicates
  ungroup() %>% 
  rename(unit_location = location_category)

vary_chars <- vary_chars %>% #the join Adds back dropped rows and apply imputation (carry forward) in the parent table 
  left_join(location %>% select(hospital_block_id, time_block, unit_location), by = c("hospital_block_id", "time_block")) %>% 
  group_by(hospital_block_id) %>% 
  arrange(time_block, .by_group = TRUE) %>% 
  tidyr::fill(unit_location, .direction = "down") %>% 
  ungroup() 


##Code status
code_status <- vary_chars %>% 
  distinct(patient_id, t_0) %>% 
  left_join(clif_code_status %>% 
              collect() %>%
              distinct(),
            by = "patient_id", 
            relationship = "many-to-many") %>%
  group_by(patient_id) %>%
  arrange(start_dttm) %>% 
  mutate(end_dttm = lead(start_dttm)) %>% #Need end_dttm (see above) -> set start_dttm of the next code status as end_dttm of current
  filter(start_dttm <= t_0 + hours(24),
         is.na(end_dttm) | end_dttm >= t_0) %>%
  ungroup() %>%
  select(patient_id, start_dttm, end_dttm, code_status_category, code_status_name)

code_status <- vary_chars %>% 
  select(patient_id, hospital_block_id, time_block, block_end) %>%
  left_join(code_status,
            by = "patient_id",
            relationship = "many-to-many") %>% 
  filter(start_dttm <= block_end, #Find the time of the location corresponding to the end of the given 2-hour block. 
         is.na(end_dttm) | end_dttm > block_end) %>% 
  group_by(hospital_block_id, time_block) %>% 
  arrange(desc(start_dttm), .by_group = TRUE) %>%
  slice(1) %>% 
  ungroup()

vary_chars <- vary_chars %>% #Add back dropped rows and apply imputation (carry forward)
  left_join(code_status %>% select(-patient_id, -block_end, -code_status_name), by = c("hospital_block_id", "time_block")) %>% 
  group_by(hospital_block_id) %>% 
  arrange(time_block, .by_group = TRUE) %>% 
  tidyr::fill(code_status_category, .direction = "down") %>% 
  select(-start_dttm, -end_dttm) %>%
  ungroup() 

#Borrow code status from admission to t_0 for remaining leading missing code_status_category
code_status_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_code_status %>% 
              collect() %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(start_dttm, .by_group = TRUE) %>%
  mutate(end_dttm = lead(start_dttm)) %>%
  filter(start_dttm >= admission_dttm & start_dttm <= t_0) %>%
  arrange(desc(start_dttm), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(patient_id, hospital_block_id, code_status_category_t0 = code_status_category)

vary_chars <- vary_chars %>% 
  left_join(code_status_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(code_status_category = if_else(
    is.na(code_status_category) & cumsum(!is.na(code_status_category)) == 0,
    code_status_category_t0,
    code_status_category
  )) %>%
  ungroup() %>%
  select(-code_status_category_t0)


##Temperature
avg_temp <- vary_chars %>% 
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "temp_c",
                     vital_value >= 32 & vital_value <= 44) %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(recorded_dttm >= block_start & recorded_dttm < block_end) %>% #Rule change: find the average of all values in the given time block.
  summarise(temp = if_else(all(is.na(vital_value)), NA_real_, mean(vital_value, na.rm = TRUE)),
            .groups = "drop") 

vary_chars <- vary_chars %>% 
  left_join(avg_temp, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(temp, .direction = "down") %>%
  ungroup()


#Borrow temp from admission to t_0 for remaining leading missing temp
temp_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "temp_c",
                     vital_value >= 32 & vital_value <= 44) %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  filter(recorded_dttm >= admission_dttm & recorded_dttm <= t_0) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(temp_t0 = mean(vital_value, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(temp_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(temp = if_else(
    is.na(temp) & cumsum(!is.na(temp)) == 0,
    temp_t0,
    temp
  )) %>%
  ungroup() %>%
  select(-temp_t0)


##MAP
avg_map <- vary_chars %>% 
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "map",
                     vital_value >= 30 & vital_value <= 250) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(recorded_dttm >= block_start & recorded_dttm < block_end) %>% #find the average of all values in the given time block.
  summarise(avg_map = if_else(all(is.na(vital_value)), NA_real_, mean(vital_value, na.rm = TRUE)),
            .groups = "drop") 

vary_chars <- vary_chars %>% 
  left_join(avg_map, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_map, .direction = "down") %>%
  ungroup()

#Borrow MAP from admission to t_0 for remaining leading missing avg_map
map_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "map",
                     vital_value >= 30 & vital_value <= 250) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  filter(recorded_dttm >= admission_dttm & recorded_dttm <= t_0) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_map_t0 = mean(vital_value, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(map_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_map = if_else(
    is.na(avg_map) & cumsum(!is.na(avg_map)) == 0,
    avg_map_t0,
    avg_map
  )) %>%
  ungroup() %>%
  select(-avg_map_t0)


##resp device
resp_support <- vary_chars %>% 
  left_join(clif_respiratory_support %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     device_category) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(recorded_dttm >= block_start & recorded_dttm < block_end) %>% 
  arrange(desc(recorded_dttm), .by_group = TRUE) %>% #Find the time of the location corresponding to the end of the given 2-hour block.
  slice(1) %>% 
  ungroup() %>% 
  select(patient_id, 
         hospital_block_id,
         time_block,
         device_category)

vary_chars <- vary_chars %>% 
  left_join(resp_support, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(device_category, .direction = "down") %>%
  ungroup()

#Borrow resp device from admission to t_0 for remaining leading missing device_category
resp_support_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_respiratory_support %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     device_category) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  filter(recorded_dttm >= admission_dttm & recorded_dttm <= t_0) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  arrange(desc(recorded_dttm), .by_group = TRUE) %>%
  slice(1) %>% 
  ungroup() %>%
  select(patient_id, hospital_block_id, device_category_t0 = device_category)

vary_chars <- vary_chars %>% 
  left_join(resp_support_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(device_category = if_else(
    is.na(device_category) & cumsum(!is.na(device_category)) == 0,
    device_category_t0,
    device_category
  )) %>%
  ungroup() %>%
  select(-device_category_t0)


##heart rate
avg_hr <- vary_chars %>% 
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "heart_rate",
                     vital_value >= 0 & vital_value <= 300) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(recorded_dttm >= block_start & recorded_dttm < block_end) %>% 
  summarise(avg_hr = if_else(all(is.na(vital_value)), NA_real_, mean(vital_value, na.rm = TRUE)),
            .groups = "drop") 

vary_chars <- vary_chars %>% 
  left_join(avg_hr, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_hr, .direction = "down") %>%
  ungroup() 

#Borrow heart rate from admission to t_0 for remaining leading missing avg_hr
hr_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_vitals %>% 
              collect() %>% 
              select(patient_id,
                     recorded_dttm, 
                     vital_category,
                     vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "heart_rate",
                     vital_value >= 0 & vital_value <= 300) %>%
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  filter(recorded_dttm >= admission_dttm & recorded_dttm <= t_0) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_hr_t0 = mean(vital_value, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(hr_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_hr = if_else(
    is.na(avg_hr) & cumsum(!is.na(avg_hr)) == 0,
    avg_hr_t0,
    avg_hr
  )) %>%
  ungroup() %>%
  select(-avg_hr_t0)


#Sodium
SODIUM_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "sodium"]
SODIUM_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "sodium"]

sodium <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "sodium") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    sodium = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         sodium >= SODIUM_MIN & sodium <= SODIUM_MAX) %>% 
  summarise(avg_sodium = if_else(all(is.na(sodium)), NA_real_, mean(sodium, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(sodium, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_sodium, .direction = "down") %>%
  ungroup() 

#Borrow sodium from admission to t_0 for remaining missing avg_sodium
sodium_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "sodium") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    sodium = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         sodium >= SODIUM_MIN & sodium <= SODIUM_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_sodium_t0 = mean(sodium, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(sodium_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_sodium = if_else(
    is.na(avg_sodium) & cumsum(!is.na(avg_sodium)) == 0,
    avg_sodium_t0,
    avg_sodium
  )) %>%
  ungroup() %>%
  select(-avg_sodium_t0)


##potassium
POTASSIUM_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "potassium"]
POTASSIUM_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "potassium"]

potassium <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "potassium") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    potassium = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         potassium >= POTASSIUM_MIN & potassium <= POTASSIUM_MAX) %>% 
  summarise(avg_potassium = if_else(all(is.na(potassium)), NA_real_, mean(potassium, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(potassium, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_potassium, .direction = "down") %>%
  ungroup() 

#Borrow potassium from admission to t_0 for remaining missing avg_potassium
potassium_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "potassium") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    potassium = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         potassium >= POTASSIUM_MIN & potassium <= POTASSIUM_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_potassium_t0 = mean(potassium, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(potassium_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_potassium = if_else(
    is.na(avg_potassium) & cumsum(!is.na(avg_potassium)) == 0,
    avg_potassium_t0,
    avg_potassium
  )) %>%
  ungroup() %>%
  select(-avg_potassium_t0)


##WBC
WBC_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "wbc"]
WBC_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "wbc"]

wbc <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "wbc") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    wbc = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         wbc >= WBC_MIN & wbc <= WBC_MAX) %>% 
  summarise(avg_wbc = if_else(all(is.na(wbc)), NA_real_, mean(wbc, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(wbc, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_wbc, .direction = "down") %>%
  ungroup()

#Borrow WBC from admission to t_0 for remaining missing avg_wbc
wbc_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "wbc") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    wbc = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         wbc >= WBC_MIN & wbc <= WBC_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_wbc_t0 = mean(wbc, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(wbc_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_wbc = if_else(
    is.na(avg_wbc) & cumsum(!is.na(avg_wbc)) == 0,
    avg_wbc_t0,
    avg_wbc
  )) %>%
  ungroup() %>%
  select(-avg_wbc_t0)

##Bicarb
BICARB_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "bicarb"]
BICARB_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "bicarb"]

bicarb <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "bicarbonate") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    bicarb = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         bicarb >= BICARB_MIN & bicarb <= BICARB_MAX) %>% 
  summarise(avg_bicarb = if_else(all(is.na(bicarb)), NA_real_, mean(bicarb, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(bicarb, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_bicarb, .direction = "down") %>%
  ungroup()

#Borrow bicarb from admission to t_0 for remaining missing avg_bicarb
bicarb_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "bicarbonate") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    bicarb = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         bicarb >= BICARB_MIN & bicarb <= BICARB_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_bicarb_t0 = mean(bicarb, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(bicarb_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_bicarb = if_else(
    is.na(avg_bicarb) & cumsum(!is.na(avg_bicarb)) == 0,
    avg_bicarb_t0,
    avg_bicarb
  )) %>%
  ungroup() %>%
  select(-avg_bicarb_t0)

##Lactate
LACTATE_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "lactate"]
LACTATE_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "lactate"]

lactate <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "lactate") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    lactate = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         lactate >= LACTATE_MIN & lactate <= LACTATE_MAX) %>% 
  summarise(avg_lactate = if_else(all(is.na(lactate)), NA_real_, mean(lactate, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(lactate, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_lactate, .direction = "down") %>%
  ungroup()

#Borrow lactate from admission to t_0 for remaining missing avg_lactate
lactate_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "lactate") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    lactate = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         lactate >= LACTATE_MIN & lactate <= LACTATE_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_lactate_t0 = mean(lactate, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(lactate_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_lactate = if_else(
    is.na(avg_lactate) & cumsum(!is.na(avg_lactate)) == 0,
    avg_lactate_t0,
    avg_lactate
  )) %>%
  ungroup() %>%
  select(-avg_lactate_t0)


##Sofa_brain
GCS_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "gcs_total"]
GCS_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "gcs_total"]

gcs_total <- vary_chars %>% 
  left_join(clif_patient_assessments %>% 
              collect() %>% 
              select(patient_id, recorded_dttm, assessment_category, numerical_value) %>%
              filter(str_to_lower(trimws(assessment_category)) == "gcs_total") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    gcs_total = if_else(
      !is.na(recorded_dttm),
      numerical_value,
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(recorded_dttm >= block_start & recorded_dttm < block_end,
         gcs_total >= GCS_MIN & gcs_total <= GCS_MAX) %>% 
  summarise(avg_gcs = if_else(all(is.na(gcs_total)), NA_real_, mean(gcs_total, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(gcs_total, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_gcs, .direction = "down") %>%
  mutate(
    sofa_brain = case_when(
      is.na(avg_gcs) ~ NA_character_, #categorize all sofa scores
      avg_gcs >= 15 ~ "0",
      avg_gcs >= 13 & avg_gcs < 15 ~ "1",
      avg_gcs >= 10 & avg_gcs < 13 ~ "2",
      avg_gcs >= 6  & avg_gcs < 10 ~ "3",
      avg_gcs < 6 ~ "4"
    ),
    sofa_brain = if_else(
      rep(all(is.na(sofa_brain)), n()),
      "not measured",
      sofa_brain
    )
  ) %>%
  ungroup() %>% 
  select(-avg_gcs)

#Borrow GCS from admission to t_0 for remaining missing sofa_brain
gcs_total_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_patient_assessments %>% 
              collect() %>% 
              select(patient_id, recorded_dttm, assessment_category, numerical_value) %>%
              filter(str_to_lower(trimws(assessment_category)) == "gcs_total") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    gcs_total = if_else(
      !is.na(recorded_dttm),
      numerical_value,
      NA_real_
    )
  ) %>% 
  filter(recorded_dttm >= admission_dttm & recorded_dttm <= t_0,
         gcs_total >= GCS_MIN & gcs_total <= GCS_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_gcs_t0 = mean(gcs_total, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    sofa_brain_t0 = case_when(
      is.na(avg_gcs_t0) ~ NA_character_,
      avg_gcs_t0 >= 15 ~ "0",
      avg_gcs_t0 >= 13 & avg_gcs_t0 < 15 ~ "1",
      avg_gcs_t0 >= 10 & avg_gcs_t0 < 13 ~ "2",
      avg_gcs_t0 >= 6  & avg_gcs_t0 < 10 ~ "3",
      avg_gcs_t0 < 6 ~ "4"
    )
  ) %>%
  select(patient_id, hospital_block_id, sofa_brain_t0)

vary_chars <- vary_chars %>% 
  left_join(gcs_total_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(sofa_brain = if_else(
    is.na(sofa_brain) & cumsum(!is.na(sofa_brain)) == 0,
    sofa_brain_t0,
    sofa_brain
  )) %>%
  ungroup() %>%
  select(-sofa_brain_t0)


#Sofa_liver
BILIRUBIN_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "bilirubin_total"]
BILIRUBIN_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "bilirubin_total"]

bilirubin_total <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "bilirubin_total") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    bilirubin_total = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         bilirubin_total >= BILIRUBIN_MIN & bilirubin_total <= BILIRUBIN_MAX) %>% 
  summarise(bilirubin_total = if_else(all(is.na(bilirubin_total)), NA_real_, mean(bilirubin_total, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(bilirubin_total, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(bilirubin_total, .direction = "down") %>%
  mutate(
    sofa_liver = case_when(
      is.na(bilirubin_total) ~ NA_character_,
      bilirubin_total < 1.2 ~ "0",
      bilirubin_total >= 1.2 & bilirubin_total < 2.0 ~ "1",
      bilirubin_total >= 2.0 & bilirubin_total < 6.0 ~ "2",
      bilirubin_total >= 6.0 & bilirubin_total < 12.0 ~ "3",
      bilirubin_total >= 12.0 ~ "4"
    ),
    sofa_liver = if_else(
      rep(all(is.na(sofa_liver)), n()),
      "not measured",
      sofa_liver
    )
  ) %>%
  ungroup() %>% 
  select(-bilirubin_total)

#Borrow bilirubin from admission to t_0 for remaining missing sofa_liver
bilirubin_total_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "bilirubin_total") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    bilirubin_total = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         bilirubin_total >= BILIRUBIN_MIN & bilirubin_total <= BILIRUBIN_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(bilirubin_total_t0 = mean(bilirubin_total, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    sofa_liver_t0 = case_when(
      is.na(bilirubin_total_t0) ~ NA_character_,
      bilirubin_total_t0 < 1.2 ~ "0",
      bilirubin_total_t0 >= 1.2 & bilirubin_total_t0 < 2.0 ~ "1",
      bilirubin_total_t0 >= 2.0 & bilirubin_total_t0 < 6.0 ~ "2",
      bilirubin_total_t0 >= 6.0 & bilirubin_total_t0 < 12.0 ~ "3",
      bilirubin_total_t0 >= 12.0 ~ "4"
    )
  ) %>%
  select(patient_id, hospital_block_id, sofa_liver_t0)

vary_chars <- vary_chars %>% 
  left_join(bilirubin_total_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(sofa_liver = if_else(
    is.na(sofa_liver) & cumsum(!is.na(sofa_liver)) == 0,
    sofa_liver_t0,
    sofa_liver
  )) %>%
  ungroup() %>%
  select(-sofa_liver_t0)


#Sofa_coagulation
PLATELET_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "platelet_count"]
PLATELET_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "platelet_count"]

platelet_count <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "platelet_count") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    platelet_count = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         platelet_count >= PLATELET_MIN & platelet_count <= PLATELET_MAX) %>% 
  summarise(platelet_count = if_else(all(is.na(platelet_count)), NA_real_, mean(platelet_count, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(platelet_count, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(platelet_count, .direction = "down") %>%
  mutate(
    sofa_coagulation = case_when(
      is.na(platelet_count) ~ NA_character_,
      platelet_count >= 150 ~ "0",
      platelet_count >= 100 & platelet_count < 150 ~ "1",
      platelet_count >= 50 & platelet_count < 100 ~ "2",
      platelet_count >= 20 & platelet_count < 50 ~ "3",
      platelet_count < 20 ~ "4"
    ),
    sofa_coagulation = if_else(
      rep(all(is.na(sofa_coagulation)), n()),
      "not measured",
      sofa_coagulation
    )
  ) %>%
  ungroup() %>% 
  select(-platelet_count)

#Borrow platelet count from admission to t_0 for remaining missing sofa_coagulation
platelet_count_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "platelet_count") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    platelet_count = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         platelet_count >= PLATELET_MIN & platelet_count <= PLATELET_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(platelet_count_t0 = mean(platelet_count, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    sofa_coagulation_t0 = case_when(
      is.na(platelet_count_t0) ~ NA_character_,
      platelet_count_t0 >= 150 ~ "0",
      platelet_count_t0 >= 100 & platelet_count_t0 < 150 ~ "1",
      platelet_count_t0 >= 50 & platelet_count_t0 < 100 ~ "2",
      platelet_count_t0 >= 20 & platelet_count_t0 < 50 ~ "3",
      platelet_count_t0 < 20 ~ "4"
    )
  ) %>%
  select(patient_id, hospital_block_id, sofa_coagulation_t0)

vary_chars <- vary_chars %>% 
  left_join(platelet_count_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(sofa_coagulation = if_else(
    is.na(sofa_coagulation) & cumsum(!is.na(sofa_coagulation)) == 0,
    sofa_coagulation_t0,
    sofa_coagulation
  )) %>%
  ungroup() %>%
  select(-sofa_coagulation_t0)


#Sofa_kidney
CREATININE_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "creatinine"]
CREATININE_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "creatinine"]

creatinine <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "creatinine") %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    creatinine = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         creatinine >= CREATININE_MIN & creatinine <= CREATININE_MAX) %>% 
  summarise(avg_creatinine = if_else(all(is.na(creatinine)), NA_real_, mean(creatinine, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(creatinine, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_creatinine, .direction = "down") %>%
  mutate(
    sofa_kidney = case_when(
      is.na(avg_creatinine) ~ NA_character_,
      avg_creatinine < 1.2 ~ "0",
      avg_creatinine >= 1.2 & avg_creatinine < 2.0 ~ "1",
      avg_creatinine >= 2.0 & avg_creatinine < 3.5 ~ "2",
      avg_creatinine >= 3.5 & avg_creatinine < 5.0 ~ "3",
      avg_creatinine >= 5.0 ~ "4"
    ),
    sofa_kidney = if_else(
      rep(all(is.na(sofa_kidney)), n()),
      "not measured",
      sofa_kidney
    )
  ) %>%
  ungroup() %>% 
  select(-avg_creatinine)

#Borrow creatinine from admission to t_0 for remaining missing sofa_kidney
creatinine_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) == "creatinine") %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    creatinine = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         creatinine >= CREATININE_MIN & creatinine <= CREATININE_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_creatinine_t0 = mean(creatinine, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    sofa_kidney_t0 = case_when(
      is.na(avg_creatinine_t0) ~ NA_character_,
      avg_creatinine_t0 < 1.2 ~ "0",
      avg_creatinine_t0 >= 1.2 & avg_creatinine_t0 < 2.0 ~ "1",
      avg_creatinine_t0 >= 2.0 & avg_creatinine_t0 < 3.5 ~ "2",
      avg_creatinine_t0 >= 3.5 & avg_creatinine_t0 < 5.0 ~ "3",
      avg_creatinine_t0 >= 5.0 ~ "4"
    )
  ) %>%
  select(patient_id, hospital_block_id, sofa_kidney_t0)

vary_chars <- vary_chars %>% 
  left_join(creatinine_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(sofa_kidney = if_else(
    is.na(sofa_kidney) & cumsum(!is.na(sofa_kidney)) == 0,
    sofa_kidney_t0,
    sofa_kidney
  )) %>%
  ungroup() %>%
  select(-sofa_kidney_t0)


#PCO2  
PCO2_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "pco2_arterial"]
PCO2_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "pco2_arterial"]

pco2 <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) %in% c("pco2_arterial", "pco2_venous")) %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    pco2 = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         pco2 >= PCO2_MIN & pco2 <= PCO2_MAX) %>% 
  summarise(avg_pco2 = if_else(all(is.na(pco2)), NA_real_, mean(pco2, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(pco2, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_pco2, .direction = "down") %>%
  ungroup()

#Borrow PCO2 from admission to t_0 for remaining missing avg_pco2
pco2_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) %in% c("pco2_arterial", "pco2_venous")) %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    pco2 = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(lab_collect_dttm >= admission_dttm & lab_collect_dttm <= t_0,
         pco2 >= PCO2_MIN & pco2 <= PCO2_MAX) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(avg_pco2_t0 = mean(pco2, na.rm = TRUE),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(pco2_t0,
            by = c("patient_id",
                   "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_pco2 = if_else(
    is.na(avg_pco2) & cumsum(!is.na(avg_pco2)) == 0,
    avg_pco2_t0,
    avg_pco2
  )) %>%
  ungroup() %>%
  select(-avg_pco2_t0)


#PH
PH_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "ph_venous"]
PH_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "ph_venous"]

ph <- vary_chars %>% 
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) %in% c("ph_arterial", "ph_venous")) %>% 
              distinct(),
            by = "patient_id",
            relationship = "many-to-many") %>%
  mutate(
    ph = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  group_by(patient_id, 
           hospital_block_id, 
           time_block) %>% 
  filter(lab_collect_dttm >= block_start & lab_collect_dttm < block_end,
         ph >= PH_MIN & ph <= PH_MAX) %>% 
  summarise(avg_ph = if_else(all(is.na(ph)), NA_real_, mean(ph, na.rm = TRUE)),
            .groups = "drop")

vary_chars <- vary_chars %>% 
  left_join(ph, 
            by = c("patient_id",
                   "hospital_block_id",
                   "time_block")) %>%
  group_by(hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  tidyr::fill(avg_ph, .direction = "down") %>%
  mutate(
    avg_ph = as.character(avg_ph),
    avg_ph = if_else(
      rep(all(is.na(avg_ph)), n()),
      "not measured",
      avg_ph
    )
  ) %>%
  ungroup()

ph_t0 <- vary_chars %>% 
  filter(time_block == 0) %>%
  select(patient_id, hospital_block_id, hospitalization_id, t_0) %>%
  distinct() %>%
  left_join(clif_hospitalization %>%
              select(patient_id, hospitalization_id, admission_dttm) %>%
              collect() %>%
              distinct(),
            by = c("patient_id", "hospitalization_id")) %>%
  left_join(clif_labs %>% 
              collect() %>% 
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              filter(str_to_lower(trimws(lab_category)) %in% c("ph_arterial", "ph_venous")) %>% 
              distinct(),
            by = "patient_id") %>%
  mutate(
    ph = if_else(
      !is.na(lab_collect_dttm),
      coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value))),
      NA_real_
    )
  ) %>% 
  filter(
    lab_collect_dttm >= admission_dttm,
    lab_collect_dttm <= t_0,
    ph >= PH_MIN,
    ph <= PH_MAX
  ) %>% 
  group_by(patient_id, hospital_block_id) %>% 
  summarise(
    avg_ph_t0 = mean(ph, na.rm = TRUE),
    .groups = "drop"
  )

vary_chars <- vary_chars %>% 
  left_join(ph_t0,
            by = c("patient_id", "hospital_block_id")) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(time_block, .by_group = TRUE) %>%
  mutate(avg_ph = if_else(
    is.na(avg_ph) & cumsum(!is.na(avg_ph)) == 0,
    as.character(avg_ph_t0),
    avg_ph
  )) %>%
  ungroup() %>%
  select(-avg_ph_t0)


##SF ratio
SPO2_MIN = outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "spo2"]
SPO2_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "spo2"]

FIO2_MIN <- outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "fio2_set"]
FIO2_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "fio2_set"]

#All respiratory support entries of the cohort
SF1 <- clif_respiratory_support %>%
  select(patient_id, recorded_dttm, fio2_set) %>%
  distinct() %>%
  collect() 

#Intermediate DF2: extract fio2 values, match them to time blocks, and impute withiin the same device period
SF2 <- vary_chars %>%
  select(patient_id, device_category, t_0, time_block, block_start, block_end) %>% 
  #Interval join to filter recorded fio2 values per each time block per patient. 
  left_join(SF1, 
            by = join_by(patient_id, 
                         block_start <= recorded_dttm, 
                         block_end > recorded_dttm)) %>% 
  #Create a cumsum column for changes in the resp. device that the patient uses in the style of hospitalization block
  group_by(patient_id) %>% 
  arrange(block_start, recorded_dttm) %>% 
  mutate(
    device_category_clean = coalesce(device_category, "UNKNOWN"), #All NA values for device_category will be treated as "UNKNOWN"
    new_device_period = case_when(
      #Initiate first record as T
      row_number() == 1 ~ TRUE, 
      #If they have multiple records, compare current record with the previous one, set new device as T (include: NA -> device but not vice versa)
      device_category_clean != lag(device_category_clean) & (device_category_clean != "UNKNOWN") ~ TRUE, 
      TRUE ~ FALSE
    ),
    device_period = cumsum(new_device_period) #a counter variable tracking the cum. sum of total devices used by the patient
  ) %>% 
  ungroup() %>% 
  select(-device_category_clean, -new_device_period) %>% 
  #Impute within the device period forewards and then backwards
  group_by(patient_id, device_period) %>% 
  arrange(recorded_dttm, .by_group = T) %>% 
  tidyr::fill(fio2_set, .direction = "down") %>%
  tidyr::fill(fio2_set, .direction = "up") %>% 
  ungroup() 

#Intermediate DF3: handle edge cases: Although SF ratio for a given time block is calculated as the average of individual SF ratios, there is one and
#only one device per time block. Therefore, We will limit fio2 records within any given block to those from the last device at the end 
#of the block. 
SF3 <- SF2 %>% 
  rename(device_category_transformed = device_category) %>%
  left_join(clif_respiratory_support %>% select(patient_id, recorded_dttm, device_category_og = device_category) %>% collect() %>% distinct(), 
            by = c("patient_id", "recorded_dttm"), 
            relationship = "many-to-many") %>% 
  filter((is.na(device_category_og) & is.na(device_category_transformed)) | device_category_og == device_category_transformed) 

#Intermediate DF4: create a time block for fio2 so that it can be matched to spo2 and filter only permissible values of fio2
SF4 <- SF3 %>% 
  select(-device_category_og) %>% 
  rename(device_category = device_category_transformed,
         fio2_time_start = recorded_dttm) %>% 
  filter(!is.na(fio2_time_start),
         FIO2_MIN <= fio2_set & FIO2_MAX >= fio2_set) %>% 
  group_by(patient_id) %>% 
  arrange(fio2_time_start) %>% 
  #If it is the last fio2 for the given hospitalization id within the 24 hour block, set the end of the fio2 block as 24 hours after t_0. 
  mutate(fio2_time_end = lead(fio2_time_start),
         fio2_time_end = if_else(row_number() == n(), t_0 + hours(24), fio2_time_end)) %>% 
  ungroup() %>% 
  left_join(clif_vitals %>%
              collect() %>% 
              select(patient_id, recorded_dttm, vital_category, vital_value) %>%
              filter(str_to_lower(trimws(vital_category)) == "spo2",
                     vital_value >= SPO2_MIN, 
                     vital_value <= SPO2_MAX) %>% 
              distinct(), 
            by = join_by(patient_id, 
                         fio2_time_start <= recorded_dttm,
                         fio2_time_end > recorded_dttm),
            relationship = "many-to-many") 

#Intermediate DF5: calculate sf on the hospitalization id, device period, and time block level -> there is a unique sf for each unique comb.
#of hosp. id, device category, device period, and time block. 
SF5 <- SF4 %>% 
  mutate(sf = vital_value / fio2_set) %>% 
  group_by(patient_id, device_category, device_period, time_block) %>% 
  summarise(avg_sf = if_else(all(is.na(sf)), NA_real_, mean(sf, na.rm = TRUE))) %>% 
  ungroup() %>% 
  #Check to avoid having multiple devices per time_block -> select the device closest to end of the time block
  group_by(patient_id, time_block) %>% 
  slice_max(order_by = device_period, n = 1, with_ties = FALSE) %>% 
  ungroup() 

#Joining back with vary_chars, create a device block for sf values from the same device, and impute within the block. 
vary_chars_SF5 <- vary_chars %>% 
  left_join(SF5 %>% 
              select(-device_category) %>% 
              collect() %>% 
              distinct(), 
            by = c("patient_id", "time_block")) %>% 
  group_by(patient_id) %>% 
  arrange(time_block, .by_group = TRUE) %>% 
  mutate(device_period_new = 1 + cumsum(row_number() != 1 &
                                          !is.na(device_category) &
                                          (is.na(lag(device_category)) |
                                             device_category != lag(device_category)))) %>% 
  ungroup() %>% 
  group_by(patient_id, device_period_new) %>% 
  tidyr::fill(avg_sf, .direction = "down") %>%
  ungroup() %>% 
  select(-device_period)

#Final DF without auxiliary variables
vary_chars <- vary_chars_SF5 %>% 
  select(-device_period_new)


##PF ratio
PAO2_MIN = outlier_thresholds$lower_limit[outlier_thresholds$variable_name == "pao2"]
PAO2_MAX <- outlier_thresholds$upper_limit[outlier_thresholds$variable_name == "pao2"]

#All respiratory support entries of the cohort
PF1 <- clif_respiratory_support %>%
  select(patient_id, recorded_dttm, fio2_set) %>%
  distinct() %>%
  collect() 

#Intermediate DF2: extract fio2 values, match them to time blocks, and impute withiin the same device period
PF2 <- vary_chars %>%
  select(patient_id, device_category, t_0, time_block, block_start, block_end) %>% 
  #Interval join to filter recorded fio2 values per each time block per patient. 
  left_join(PF1, 
            by = join_by(patient_id, 
                         block_start <= recorded_dttm, 
                         block_end > recorded_dttm)) %>% 
  #Create a cumsum column for changes in the resp. device that the patient uses in the style of hospitalization block
  group_by(patient_id) %>% 
  arrange(block_start, recorded_dttm) %>% 
  mutate(
    device_category_clean = coalesce(device_category, "UNKNOWN"), #All NA values for device_category will be treated as "UNKNOWN"
    new_device_period = case_when(
      row_number() == 1 ~ TRUE, #Initiate the cumsum with the first record: set as T
      #If they have multiple records, compare current record with the previous one, set new device as T (include: NA -> device but not vice versa)
      device_category_clean != lag(device_category_clean) & (device_category_clean != "UNKNOWN") ~ TRUE, 
      TRUE ~ FALSE
    ),
    device_period = cumsum(new_device_period) #a counter variable tracking the cum. sum of total devices used by the patient
  ) %>% 
  ungroup() %>% 
  select(-device_category_clean, -new_device_period) %>% 
  #Impute within the device period forewards and then backwards
  group_by(patient_id, device_period) %>% 
  arrange(recorded_dttm, .by_group = T) %>% 
  tidyr::fill(fio2_set, .direction = "down") %>%
  tidyr::fill(fio2_set, .direction = "up") %>% 
  ungroup() 

#Intermediate DF3:handle edge cases: Although SF ratio for a given time block is calculated as the average of individual SF ratios, there is one and
#only one device per time block. Therefore, We will limit fio2 records within any given block to those from the last device at the end 
#of the block. 
PF3 <- PF2 %>% 
  rename(device_category_transformed = device_category) %>%
  left_join(clif_respiratory_support %>% select(patient_id, recorded_dttm, device_category_og = device_category) %>% collect() %>% distinct(), 
            by = c("patient_id", "recorded_dttm"), 
            relationship = "many-to-many") %>% 
  filter((is.na(device_category_og) & is.na(device_category_transformed)) | device_category_og == device_category_transformed)  

#Intermediate DF4: create a time block for fio2 so that it can be matched to pao2 and filter only permissible values of fio2
PF4 <- PF3 %>% 
  select(-device_category_og) %>% 
  rename(device_category = device_category_transformed,
         fio2_time_start = recorded_dttm) %>% 
  filter(!is.na(fio2_time_start),
         FIO2_MIN <= fio2_set & FIO2_MAX >= fio2_set) %>% 
  group_by(patient_id) %>% 
  arrange(fio2_time_start) %>% 
  #If it is the last fio2 for the given hospitalization id within the 24 hour block, set the end of the fio2 block as 24 hours after t_0. 
  mutate(fio2_time_end = lead(fio2_time_start),
         fio2_time_end = if_else(row_number() == n(), t_0 + hours(24), fio2_time_end)) %>% 
  ungroup() %>% 
  left_join(clif_labs %>% 
              collect() %>%
              select(patient_id, lab_collect_dttm, lab_category, lab_value, lab_value_numeric) %>%
              mutate(pao2 = coalesce(lab_value_numeric, suppressWarnings(as.numeric(lab_value)))) %>% 
              filter(str_to_lower(trimws(lab_category)) == "po2_arterial",
                     pao2 >= PAO2_MIN, 
                     pao2 <= PAO2_MAX) %>% 
              distinct(), 
            by = join_by(patient_id, 
                         fio2_time_start <= lab_collect_dttm,
                         fio2_time_end > lab_collect_dttm),
            relationship = "many-to-many") 

#Intermediate DF5: calculate pf on the hospitalization id, device period, and time block level -> there is a unique sf for each unique comb.
#of hosp. id, device category, device period, and time block. 
PF5 <- PF4 %>% 
  mutate(pf = pao2 / fio2_set) %>% 
  group_by(patient_id, device_category, device_period, time_block) %>% 
  summarise(avg_pf = if_else(all(is.na(pf)), NA_real_, mean(pf, na.rm = TRUE))) %>% 
  ungroup() %>% 
  #Check to avoid having multiple devices per time_block -> select the device closest to end of the time block
  group_by(patient_id, time_block) %>% 
  slice_max(order_by = device_period, n = 1, with_ties = FALSE) %>% 
  ungroup() 

#Joining back with vary_chars
vary_chars_check_PF <- vary_chars %>% 
  left_join(PF5 %>% 
              select(-device_category) %>% 
              collect() %>% 
              distinct(), 
            by = c("patient_id", "time_block")) %>% 
  group_by(patient_id) %>% 
  arrange(time_block, .by_group = TRUE) %>% 
  mutate(device_period_new = 1 + cumsum(row_number() != 1 &
                                          !is.na(device_category) &
                                          (is.na(lag(device_category)) |
                                             device_category != lag(device_category)))) %>% 
  ungroup() %>% 
  group_by(patient_id, device_period_new) %>% 
  tidyr::fill(avg_pf, .direction = "down") %>%
  ungroup() %>% 
  select(-device_period)

#Final DF without auxiliary variables
vary_chars <- vary_chars_check_PF %>% 
  select(-device_period_new)


##Anti hypertensive drip dose
antihypertensive_drip <- vary_chars %>% 
  left_join(clif_medication_admin_continuous %>% 
              collect() %>% 
              select(patient_id, 
                     med_category, 
                     mar_action_category,
                     admin_dttm) %>% 
              filter(str_to_lower(trimws(med_category)) %in% c("nicardipine",
                                                               "nitroglycerin",
                                                               "nitroprusside",
                                                               "labetalol",
                                                               "esmolol"),
                     str_to_lower(trimws(mar_action_category)) %in% c("new bag", 
                                                                      "rate change", 
                                                                      "restarted", 
                                                                      "rate verify")) %>% 
              distinct(), 
            by = join_by(patient_id, 
                         block_start <= admin_dttm, 
                         block_end >= admin_dttm),
            relationship = "many-to-many") %>% 
  group_by(patient_id, block_start) %>% 
  summarise(
    antihypertensive_drip = any(!is.na(med_category))
  ) %>% 
  ungroup() 

vary_chars <- vary_chars %>%
  left_join(antihypertensive_drip, by = c("patient_id", "block_start")) %>%
  mutate(antihypertensive_drip = coalesce(antihypertensive_drip, FALSE))


##Naloxone
naloxone <- vary_chars %>% 
  left_join(clif_medication_admin_continuous %>% 
              collect() %>% 
              select(patient_id, 
                     med_category, 
                     mar_action_category,
                     admin_dttm) %>% 
              filter(med_category %in% c("naloxone"),
                     str_to_lower(trimws(mar_action_category)) %in% c("new bag", 
                                                                      "rate change", 
                                                                      "restarted", 
                                                                      "rate verify")) %>% 
              distinct(), 
            by = join_by(patient_id, 
                         block_start <= admin_dttm, 
                         block_end >= admin_dttm),
            relationship = "many-to-many") %>% 
  group_by(patient_id, block_start) %>% 
  summarise(
    naloxone = any(!is.na(med_category))
  ) %>% 
  ungroup() 

vary_chars <- vary_chars %>%
  left_join(naloxone, by = c("patient_id", "block_start")) %>%
  mutate(naloxone = coalesce(naloxone, FALSE))


##Create a treatment assignment table (first time patient transitions from ed to ICU or stepdown or censored)
trt_assignment_time <- baseline_chars %>%
  select(patient_id, t_0) %>%
  mutate(
    treatment_assignment_window_end = t_0 + hours(24),
    t_negative2 = t_0 - hours(2)
  ) %>%
  left_join(
    clif_adt %>%
      collect() %>%
      select(patient_id, location_category, in_dttm, out_dttm) %>%
      filter(!str_to_lower(trimws(location_category)) %in% c("procedural", "radiology", "dialysis", "other")) %>% 
      mutate(location_category = str_to_lower(trimws(location_category))),
    by = join_by(
      patient_id,
      t_negative2 <= out_dttm,
      treatment_assignment_window_end >= in_dttm
    )
  ) %>%
  arrange(patient_id, in_dttm) %>%
  group_by(patient_id) %>%
  mutate(location_category = coalesce(location_category, "unknown")) %>%
  summarise(
    t_0 = first(t_0),
    treatment_assignment_window_end = first(treatment_assignment_window_end),
    #Maps out the location transition of a patient during t0 and t0+24
    treatment_transition_path = {
      loc <- location_category
      loc <- loc[c(TRUE, loc[-1] != loc[-length(loc)])] #Don't want repetitive location updates, e.g., ED -> Ed -> ED
      paste(loc, collapse = " -> ")
    },
    
    first_location = first(location_category),
    first_location_time = first(in_dttm),
    
    first_non_ed_idx = {
      idx <- which(location_category != "ed")
      if (length(idx) == 0) NA_integer_ else idx[1]
    },
    
    first_non_ed_location = if_else(
      is.na(first_non_ed_idx),
      NA_character_,
      location_category[first_non_ed_idx]
    ),
    
    first_non_ed_time = if_else(
      is.na(first_non_ed_idx),
      as.POSIXct(NA),
      in_dttm[first_non_ed_idx]
    ),
    
    latest_out_dttm = if_else(all(is.na(out_dttm)), 
                              as.POSIXct(NA),
                              max(out_dttm, na.rm = TRUE)),
    
    .groups = "drop"
  ) %>%
  mutate(
    treatment = case_when(
      #we expect everyone to start in the ED, adding as defensive coding
      first_location != "ed" ~ NA_integer_,
      first_non_ed_location == "icu" ~ 1L,
      first_non_ed_location == "stepdown" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    time_to_assignment = case_when(
      #not starting in ED
      first_location != "ed" ~ t_0,
      #assigned ICU/stepdown -> time = when they moved into ICU/stepdown
      first_non_ed_location %in% c("icu", "stepdown") ~ first_non_ed_time,
      #assigned a non ICU/stepdown location that is not NA, e.g., ward. 
      !is.na(first_non_ed_location) ~ first_non_ed_time,
      #only stayed in ED -> censor at earlier of 24h or latest out_dttm (possible multiple ED periods overlapping 24 hour stay, select the latest)
      TRUE ~ pmin(treatment_assignment_window_end,
                  coalesce(latest_out_dttm, treatment_assignment_window_end))
    ),
    
    randomization = as.integer(treatment %in% c(0L, 1L))
  ) %>% 
  select(patient_id, treatment_transition_path, treatment, time_to_assignment, randomization)

#Join back useful cols from above to baseline
baseline_chars <- baseline_chars %>% 
  left_join(trt_assignment_time, 
            by = "patient_id")

#Allows vary_chars to contain time periods up to ICU/Stepdown placement, 
#censoring due to discharge from ED to ward/home/other places, or t_0 + 24 hours (whichever occurs first). 
vary_chars <- vary_chars %>% 
  left_join(trt_assignment_time %>% 
              select(patient_id, time_to_assignment), by = "patient_id") %>% 
  group_by(patient_id) %>%
  arrange(block_start, .by_group = TRUE) %>%
  mutate(
    max_block_end = max(block_end, na.rm = TRUE),
    event_block = !is.na(time_to_assignment) &
      (
        (time_to_assignment >= block_start & time_to_assignment < block_end) | 
          (time_to_assignment == block_end & block_end == max_block_end)
      ),
    stop_row = cumsum(event_block)
  ) %>%
  ungroup() %>%
  filter(stop_row == 0 | event_block) %>%
  select(-event_block, -stop_row, -max_block_end)

#Shifting time block indices to represent end of the time block for transformations needed to construct assembled_df
#NOTE: previously, time block t = x <-> [x hour, x+2 hour). Now,  time block t = x <-> [x-2 hour, x hour).
vary_chars <- vary_chars %>% 
  mutate(time_block = time_block + 2)

#Remove [-2, 0) block. 
vary_chars <- vary_chars %>% 
  filter(time_block != 0)

# -----------------  End defining baseline table and varying characteristics 



# -----------------  Defining outcomes table

outcomes_chars <- baseline_chars %>% 
  select(patient_id, hospital_block_id, t_0) %>% #t_0 is start of follow-up, as defined in SAP. 
  distinct()

##In-hospital death or hospice discharge by day 28 after start of follow-up
#Hospice discharge
hospice_discharge <- outcomes_chars %>% 
  mutate(end_date = t_0 + days(28)) %>% 
  left_join(
    clif_hospitalization %>% 
      collect() %>% 
      select(patient_id, discharge_category, discharge_dttm) %>% 
      mutate(
        discharge_dttm = if_else(
          is.na(discharge_dttm),
          as.POSIXct("9999-12-30", tz = "UTC"),
          discharge_dttm
        )
      ) %>%
      filter(str_to_lower(trimws(discharge_category)) == "hospice") %>% 
      distinct(),
    by = join_by(
      patient_id, 
      t_0 <= discharge_dttm,
      end_date >= discharge_dttm
    )
  ) %>% 
  group_by(patient_id, hospital_block_id) %>%
  summarise(
    hospice_transition = as.integer(any(!is.na(discharge_dttm))),
    hospice_discharge_time_by_28d = if (any(!is.na(discharge_dttm))) {
      min(discharge_dttm, na.rm = TRUE)
    } else {
      as.POSIXct(NA)
    },
    .groups = "drop"
  )

#Mortality
mortality <- outcomes_chars %>% 
  mutate(end_date = t_0 + days(28)) %>% 
  left_join(
    clif_hospitalization %>% 
      collect() %>% 
      select(patient_id, discharge_dttm, discharge_category) %>% 
      mutate(
        discharge_dttm = if_else(
          is.na(discharge_dttm),
          as.POSIXct("9999-12-30", tz = "UTC"),
          discharge_dttm
        )
      ) %>%
      filter(str_to_lower(trimws(discharge_category)) == "expired") %>% 
      distinct(),
    by = join_by(
      patient_id,
      t_0 <= discharge_dttm,
      end_date >= discharge_dttm
    )
  ) %>%
  mutate(
    death_time = if_else(
      !is.na(discharge_dttm),
      discharge_dttm,
      as.POSIXct(NA)
    ),
    death_by_28d = as.integer(!is.na(death_time))
  ) %>%
  group_by(patient_id, hospital_block_id) %>%
  summarise(
    death_by_28d = as.integer(any(death_by_28d == 1)),
    death_time = if (any(!is.na(death_time))) {
      min(death_time, na.rm = TRUE)
    } else {
      as.POSIXct(NA)
    },
    .groups = "drop"
  )


#Alive discharge/censoring from hospital
discharge_censor <- outcomes_chars %>% 
  mutate(end_date = t_0 + days(28)) %>% 
  left_join(clif_hospitalization %>% 
              collect() %>% 
              select(patient_id, hospital_block_id, hospitalization_id,
                     admission_dttm, discharge_dttm, discharge_category) %>% 
              distinct() %>% 
              #This inner join is a guardrail, ensuring among those with multiple hospitalizations, only hospitlaizations 
              #considered a part of the same hospital stay can contribute to discharge date determination in the below code. 
              inner_join(hospital_block_key_obj %>% 
                           collect(),
                         by = c("patient_id", "hospital_block_id", "hospitalization_id"),
                         relationship = "many-to-many"),
            by = c("patient_id", "hospital_block_id")) %>%
  left_join(baseline_chars %>%
              select(patient_id, hospital_block_id, treatment, time_to_assignment) %>%
              distinct(),
            by = c("patient_id", "hospital_block_id")) %>%
  filter(admission_dttm <= end_date,
         is.na(discharge_dttm) | discharge_dttm >= t_0) %>%
  group_by(patient_id, hospital_block_id) %>%
  arrange(discharge_dttm, .by_group = TRUE) %>%
  summarise(
    t_0 = first(t_0),
    end_date = first(end_date),
    treatment = first(treatment),
    time_to_assignment = first(time_to_assignment),
    
    terminal_discharge_dttm = if (all(is.na(discharge_dttm))) {
      as.POSIXct(NA)
    } else {
      max(discharge_dttm, na.rm = TRUE) #Accounting for multiple hospitalizations during the same stay 
    },
    
    terminal_discharge_category = if (all(is.na(discharge_dttm))) {
      NA_character_
    } else {
      discharge_category[which.max(coalesce(discharge_dttm, as.POSIXct("1900-01-01", tz = "UTC")))]
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    terminal_discharge_category_clean = str_to_lower(trimws(terminal_discharge_category)),
    
    
    discharge_censor_dttm_by_28d = case_when(
      #IF any of the following conditions:
      #1. discharge_dttm is unavailable/NA for a patient
      #2. max. available discharge_dttm is not within follow-up
      #3. the discharge_category corresponding to the max. available discharge_dttm is an event is death or hospice
      #4. for placed patients, the max. available discharge_dttm < time to placement (see trt_assignment_time code block)
      #THEN: set discharge_censor_dttm_by_28d to NA
      #ELSE: set discharge_censor_dttm_by_28d to terminal_discharge_dttm
      is.na(terminal_discharge_dttm) | terminal_discharge_dttm < t_0 | terminal_discharge_dttm > end_date ~ as.POSIXct(NA),
      terminal_discharge_category_clean %in% c("expired", "hospice") ~ as.POSIXct(NA),
      treatment %in% c(0L, 1L) &
        !is.na(time_to_assignment) &
        terminal_discharge_dttm < time_to_assignment ~ as.POSIXct(NA),
      TRUE ~ terminal_discharge_dttm
    )
  ) %>%
  select(
    patient_id,
    hospital_block_id,
    discharge_censor_dttm_by_28d
  )

#Combine in-hospital mortality and hospice, and then join back
mortality_hospice <- mortality %>%
  select(patient_id, hospital_block_id, death_by_28d, death_dttm_by_28d = death_time) %>%
  left_join(
    hospice_discharge %>%
      select(
        patient_id,
        hospital_block_id,
        hospice_transition,
        hospice_discharge_dttm_by_28d = hospice_discharge_time_by_28d
      ),
    by = c("patient_id", "hospital_block_id")
  ) %>%
  mutate(
    death_or_hospice_dttm_by_28d = case_when(
      !is.na(death_dttm_by_28d) & !is.na(hospice_discharge_dttm_by_28d) ~ #If multiple hospitalizations resulting in a combination of hospice and death during 28 days, then take the min(event dates)
        pmin(death_dttm_by_28d, hospice_discharge_dttm_by_28d), 
      !is.na(death_dttm_by_28d) ~ death_dttm_by_28d, #If only death, then take the expired event date
      !is.na(hospice_discharge_dttm_by_28d) ~ hospice_discharge_dttm_by_28d, #If only hospice, then take the hospice event date
      TRUE ~ as.POSIXct(NA)
    ),
    death_or_hospice_by_28d = as.integer(!is.na(death_or_hospice_dttm_by_28d)) #Death or hospice indicator
  ) %>% 
  left_join(
    baseline_chars %>%
      select(patient_id, hospital_block_id, treatment) %>%
      distinct(),
    by = c("patient_id", "hospital_block_id")
  ) %>%
  mutate(
    #As specified, outcome events (in-hospital death or hospice or censoring due to discharge from hospital) are only counted among those who were placed to ICU or stepdown.
    #Therefore, we turn event indicators off and specify dates as NA if there was no placement for a patient who died or went into hospice.  
    death_by_28d = if_else(
      is.na(treatment),
      0L,
      death_by_28d
    ),
    death_dttm_by_28d = if_else(
      is.na(treatment),
      as.POSIXct(NA),
      death_dttm_by_28d
    ),
    hospice_transition = if_else(
      is.na(treatment),
      0L,
      hospice_transition
    ),
    hospice_discharge_dttm_by_28d = if_else(
      is.na(treatment),
      as.POSIXct(NA),
      hospice_discharge_dttm_by_28d
    ),
    death_or_hospice_dttm_by_28d = if_else(
      is.na(treatment),
      as.POSIXct(NA),
      death_or_hospice_dttm_by_28d
    ),
    death_or_hospice_by_28d = if_else(
      is.na(treatment),
      0L,
      death_or_hospice_by_28d
    )
  ) %>%
  select(-treatment)


#Combining everything
outcomes_chars <- outcomes_chars %>% 
  left_join(mortality_hospice, 
            by = c("patient_id", "hospital_block_id")) %>%
  left_join(discharge_censor,
            by = c("patient_id", "hospital_block_id")) %>%
  left_join(baseline_chars %>% 
              select(patient_id, hospital_block_id, treatment, time_to_assignment), 
            by = c("patient_id", "hospital_block_id")) %>%
  mutate(end_date = t_0 + days(28),
         #As specified, outcome events are counted only among patients who were placed in the ICU or stepdown unit. Patients without ICU/stepdown placement are not
         #considered at risk for the primary outcome.
         death_or_hospice_dttm_by_28d = case_when(
            treatment %in% c(0L, 1L) ~ death_or_hospice_dttm_by_28d,
            TRUE ~ as.POSIXct(NA)
          ),
         death_dttm_by_28d = case_when(
            treatment %in% c(0L, 1L) ~ death_dttm_by_28d,
            TRUE ~ as.POSIXct(NA)
          ),
         hospice_discharge_dttm_by_28d = case_when(
            treatment %in% c(0L, 1L) ~ hospice_discharge_dttm_by_28d,
            TRUE ~ as.POSIXct(NA)
          ),
         discharge_censor_dttm_by_28d = case_when(
           treatment %in% c(0L, 1L) &
             !is.na(discharge_censor_dttm_by_28d) &
             !is.na(time_to_assignment) &
             discharge_censor_dttm_by_28d < time_to_assignment ~ as.POSIXct(NA),
           TRUE ~ discharge_censor_dttm_by_28d)
         ) %>%
  mutate(
    #As specified, only in-hospital death, hospice discharge, or censoring due to
    #discharge from hospital within 28 days from the start of follow-up contributes to the primary outcome.

    #Note: Although intermediate variables from the mortality_hospice and discharge_censor tables have the suffix _by_28d, 
    #they are formally restricted to the 28-day follow-up window in the below section. Variable naming will be improved in the next iteration.
    
    #Final follow-up time among ICU/stepdown placed patients is the earliest of:
    #1. death or hospice discharge,
    #2. alive discharge from hospital/censoring,
    #3. administrative censoring at day 28.
    discharge_death_hospice_dttm_by_28d = case_when(
      treatment %in% c(0L, 1L) ~ pmin(
        death_or_hospice_dttm_by_28d,
        discharge_censor_dttm_by_28d,
        end_date,
        na.rm = TRUE
      ),
      TRUE ~ as.POSIXct(NA)
    ),
    
    #Composite event indicator: death or hospice counts only if it occurs on or before the final follow-up/censoring time
    death_or_hospice_by_28d = as.integer(
      !is.na(death_or_hospice_dttm_by_28d) &
        !is.na(discharge_death_hospice_dttm_by_28d) &
        death_or_hospice_dttm_by_28d <= discharge_death_hospice_dttm_by_28d
    ),
    
    #Individual event indicators
    death_by_28d = as.integer(
      !is.na(death_dttm_by_28d) &
        !is.na(discharge_death_hospice_dttm_by_28d) &
        death_dttm_by_28d <= discharge_death_hospice_dttm_by_28d
    ),
    
    hospice_transition = as.integer(
      !is.na(hospice_discharge_dttm_by_28d) &
        !is.na(discharge_death_hospice_dttm_by_28d) &
        hospice_discharge_dttm_by_28d <= discharge_death_hospice_dttm_by_28d
    ),
    
    #Set dates to NA when the corresponding event is not counted.
    death_dttm_by_28d = case_when(
      death_by_28d == 1L ~ death_dttm_by_28d,
      TRUE ~ as.POSIXct(NA)
    ),
    
    hospice_discharge_dttm_by_28d = case_when(
      hospice_transition == 1L ~ hospice_discharge_dttm_by_28d,
      TRUE ~ as.POSIXct(NA)
    ),
    
    death_or_hospice_dttm_by_28d = case_when(
      death_or_hospice_by_28d == 1L ~ death_or_hospice_dttm_by_28d,
      TRUE ~ as.POSIXct(NA)
    )
  ) %>%
  select(-treatment, -end_date, -time_to_assignment)


# -----------------  End defining outcomes table



# -----------------  Start constructing assembled_df with baseline characteristics, varying characteristics, and outcome variables for downstream analyses

assembled_df <- vary_chars %>%
  left_join(
    baseline_chars %>%
      select(patient_id, hospital_block_id, year:elixhauser_count),
    by = c("patient_id", "hospital_block_id")
  ) %>%
  bind_rows(
    vary_chars %>%
      left_join(
        baseline_chars %>%
          select(patient_id, hospital_block_id, year:elixhauser_count),
        by = c("patient_id", "hospital_block_id")
      ) %>%
      group_by(patient_id, hospital_block_id) %>%
      summarise(
        t_0 = first(t_0),
        last_time_block = max(time_block, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      #After placement (see time_to_assignment creation), expand the time periods to t_0 + 24 hours in intervals of 2 hours 
      crossing(time_block = seq(2, 24, by = 2)) %>%
      filter(time_block > last_time_block) %>%
      select(patient_id, hospital_block_id, t_0, time_block)
  ) %>% 
  bind_rows(
    vary_chars %>%
      left_join(
        baseline_chars %>%
          select(patient_id, hospital_block_id, year:elixhauser_count),
        by = c("patient_id", "hospital_block_id")
      ) %>%
      group_by(patient_id, hospital_block_id) %>%
      summarise(
        t_0 = first(t_0),
        last_time_block = max(time_block, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      #After t_0 + 24 hours, expand the time periods to t_0 + 28 days of in intervals of one days. 
      crossing(time_block = seq(48, 672, by = 24)) %>%
      select(patient_id, hospital_block_id, t_0, time_block)
  ) %>% 
  arrange(patient_id, hospital_block_id, time_block)


assembled_df <- assembled_df %>% 
  group_by(patient_id, hospital_block_id) %>%
  mutate(
    block_end = t_0 + hours(time_block),
    block_start = lag(block_end, default = first(t_0))
  ) %>%
  #fill in based on the last available value at time of assignment
  tidyr::fill(
    -c(
      time_block,
      block_start,
      block_end,
      hospitalization_id
    ),
    .direction = "down"
  ) %>%
  ungroup() %>%
  left_join(
    baseline_chars %>% 
      select(
        patient_id,
        hospital_block_id,
        treatment,
        treatment_transition_path) %>% 
      distinct(), 
    by = c("patient_id", "hospital_block_id")
  ) %>% 
  left_join(
    outcomes_chars %>%
      select(
        patient_id,
        hospital_block_id,
        death_or_hospice_dttm_by_28d,
        discharge_death_hospice_dttm_by_28d
      ) %>%
      distinct(),
    by = c("patient_id", "hospital_block_id")
  ) %>%
  mutate(
    #time to any event: if never placed, then it is time to discharge to ward, home, or t_0 + 24 hours (whichever occurs first)
    #time to any event: if placed, then it is time to death, hospice, discharge from hospital, or adminsitrative censoring 
    time_to_event = case_when(
      is.na(treatment) ~ time_to_assignment,
      TRUE ~ discharge_death_hospice_dttm_by_28d
    )
  ) %>%
  filter(block_start < time_to_event) %>%
  mutate(
    death_or_hospice_by_28d = if_else(
      !is.na(death_or_hospice_dttm_by_28d) &
        death_or_hospice_dttm_by_28d > block_start &
        death_or_hospice_dttm_by_28d <= block_end,
      1L,
      0L
    ),
    discharge_death_hospice_dttm_by_28d = if_else(
      is.na(treatment),
      as.POSIXct(NA, tz = "UTC"),
      discharge_death_hospice_dttm_by_28d
    ),
    time_to_censor_from_treatment = if_else(
      is.na(treatment),
      time_to_assignment,
      as.POSIXct(NA, tz = "UTC")
    ),
    time_to_treatment = if_else(
      is.na(treatment),
      as.POSIXct(NA, tz = "UTC"),
      time_to_assignment
    )) %>% 
  mutate(treatment_varying = case_when(treatment == 0L ~ treatment,
                                       treatment == 1L & time_to_assignment > block_end ~ 0L,
                                       treatment == 1L & time_to_assignment <= block_end ~ 1L,
                                       TRUE ~ NA_integer_)) %>% 
  rename(treatment_indicator = treatment) 


assembled_df <- assembled_df %>%
  rename(avg_temp = temp) %>% 
  select(
    #Patient identifiers
    patient_id,
    hospital_block_id,
    hospitalization_id,
    
    #Person period identifiers
    t_0,
    time_block,
    block_start,
    block_end,
    
    #Varying characteristics
    unit_location,
    code_status_category,
    device_category,
    avg_temp,
    avg_map,
    avg_hr,
    avg_sodium,
    avg_potassium,
    avg_wbc,
    avg_bicarb,
    avg_lactate,
    avg_pco2,
    avg_ph,
    avg_sf,
    avg_pf,
    sofa_brain,
    sofa_liver,
    sofa_coagulation,
    sofa_kidney,
    antihypertensive_drip,
    naloxone,
    
    #Baseline variables
    year,
    admission_dttm,
    age,
    albumin_baseline,
    pulmonary_vasc,
    renal_failure,
    liver_disease,
    metastatic_cancer,
    elixhauser_index,
    elixhauser_count,
    
    #Treatment variables
    treatment_indicator, #Fixed ICU vs. Stepdown indicator at the patient-level. NA if no ICU or stepdown placement.
    treatment_transition_path, #Fixed treatment transition path at the patient-level (e.g., from ED to Ward or ED to ICU), the number of NA should be 0.
    treatment_varying, #Time-varying treatment initiation variable: among those with eventual ICU or Stepdown placement, 0 for person-time-periods not placed to ICU, 1 for ICU periods. NA if no ICU or stepdown placement.
    time_to_treatment, #Fixed dttm variable for ICU or Stepdown placement, NA if neither.
    time_to_censor_from_treatment, #Fixed dttm variable for patients with No ICU nor Stepdown placement, NA if placed to either ICU or Stepdown.
    
    #Outcome variables
    death_or_hospice_by_28d, #Time-varying outcome variable: among those with eventual ICU or Stepdown placement, 0 = event not yet observed during person-time-period, 1 = event observed. NA if no ICU or stepdown placement.
                             #Patient stop contributing to the risk set once this variable becomes 1 for those who experienced event.
    discharge_death_hospice_dttm_by_28d, #Fixed dttm outcome variable for death, hospice, or lost to follow up (discharged from hospital), NA if no ICU or stepdown placement.
    death_or_hospice_dttm_by_28d #Fixed dttm outcome variable for death or hospice,, NA if no ICU or stepdown placement.
  )
  



# -----------------  Start constructing assembled_df with baseline characteristics, varying characteristics, and outcome variables for downstream analyses



# -----------------  Start exporting data

#Save baseline characteristics
write_csv(baseline_chars, paste0(project_location, "/private_tables/baseline_chars.csv"))

#Save varying characteristics
write_csv(vary_chars, paste0(project_location, "/private_tables/vary_chars.csv"))

#Save outcome variables
write_csv(outcomes_chars, paste0(project_location, "/private_tables/outcomes_chars.csv"))

#Save assembled_df for analysis
write_csv(assembled_df, paste0(project_location, "/private_tables/assembled_df.csv"))


# -----------------  End exporting data


