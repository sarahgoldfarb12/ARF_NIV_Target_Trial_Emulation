###QC functions
check_vary <- function(
    data,
    ids,
    id_col = "patient_id",
    print_n = 100,
    exclude_cols = NULL,
    group = c(id_col, "time_block")
) {
  if (!id_col %in% names(data)) {
    stop(paste0("Column `", id_col, "` not found in data."))
  }
  
  if (missing(ids) || is.null(ids) || length(ids) == 0) {
    stop("Please provide a deterministic vector of IDs via `ids`.")
  }
  
  data %>%
    filter(.data[[id_col]] %in% ids) %>%
    select(-any_of(exclude_cols)) %>%
    arrange(across(all_of(group))) %>%
    print(n = print_n, width = Inf)
}


#Set seed for reproducibility
set.seed(12)


#10 sample ids for QC checks
sample_ids <- vary_chars %>%
  distinct(patient_id) %>%
  slice_sample(n = 20) %>%
  pull(patient_id)


###1.Randomization and censor_from_treatmenting 

exclude_cols <- names(vary_chars)[8:length(names(vary_chars))]

##A. missingness check: all columns except treatment should be 0
trt_assignment_time %>%
  summarise(
    across(
      -patient_id,
      ~ sum(is.na(.)),
      .names = "n_missing_{.col}"
    )
  )


##B. manually compare extracted treatment transition with actual treatment trajectory from clif-adt
check_vary(data = clif_adt %>% select(patient_id, in_dttm, location_category) %>% collect() %>% inner_join(baseline_chars %>% 
                                                                                                             select(patient_id, t_0),
                                                                                                           by = "patient_id"
),
          ids = sample_ids,
          print_n = 100, 
          exclude_cols = NULL, 
          group = c("patient_id")) %>% View()

check_vary(data = trt_assignment_time,
            ids = sample_ids,
            print_n = 100,
            exclude_cols = NULL,
            group = c("patient_id"))


##C. check treatment_transition_path and randomization indicator consistency
qc_t1 <- table(trt_assignment_time$treatment_transition_path, trt_assignment_time$randomization, useNA = "ifany") %>% as.data.frame.matrix()
bind_rows(qc_t1, 
          tibble(`0` = colSums(qc_t1)[1],
                 `1` = colSums(qc_t1)[2]))

##D. check trt transition, trt indicator, and trt assignment time are consistent
table(trt_assignment_time$treatment_transition_path, trt_assignment_time$treatment, useNA = "ifany") 
table(trt_assignment_time$treatment_transition_path, trt_assignment_time$treatment, useNA = "ifany") %>% colSums()
table(trt_assignment_time$treatment_transition_path, trt_assignment_time$treatment, useNA = "ifany") %>% colSums() %>% sum()

##E.Restricting vary_chars to before each individual's time to treatment initiation or censor_from_treatmenting
#Checking dim
dim(vary_chars)

#Check missingness: all columns except treatment should be 0
qc_summary <- bind_rows(
  #Check missingness: all columns except treatment should be 0
  vary_chars %>%
    summarise(
      check = "missingness_except_treatment",
      pass =
        sum(is.na(patient_id)) == 0 &
        sum(is.na(time_to_assignment)) == 0 &
        sum(is.na(censor)) == 0 &
        sum(is.na(at_risk)) == 0
    ),
  
  #Check any row after treatment/censor time: should be 0
  vary_chars %>%
    summarise(
      check = "no_rows_after_event",
      pass = sum(block_start > time_to_assignment, na.rm = TRUE) == 0
    ),
  
  #Check trt and censor should not happen on same row: condition should be T
  vary_chars %>%
    summarise(
      check = "no_treatment_censor_overlap",
      pass = sum(!is.na(treatment) & censor == 1L, na.rm = TRUE) == 0
    ),
  
  #Check difference between number of rows and number of missing treatment should be total number of randomization:
  #2nd to 5th number should be the same + conditions should be TRUE
  vary_chars %>%
    summarise(
      check = "treatment_count_matches_randomization",
      pass =
        sum(!is.na(treatment)) == sum(trt_assignment_time$randomization == 1L, na.rm = TRUE) &
        sum(!is.na(treatment)) == sum(!is.na(trt_assignment_time$treatment))
    ),
  
  #Check treatment distribution: Condition should True
  vary_chars %>%
    summarise(
      check = "censor_count_matches_trt_assignment",
      pass = sum(censor == 1L, na.rm = TRUE) == sum(is.na(trt_assignment_time$treatment))
    ),
  
  #Check event uniqueness: each patient should have should terminal row ie either treatment or censor
  vary_chars %>%
    group_by(patient_id) %>%
    summarise(
      n_terminal_rows = sum(!is.na(max(treatment)) | max(censor) == 1L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    summarise(
      check = "one_terminal_row_per_patient",
      pass =
        n() == n_distinct(baseline_chars$patient_id) &
        sum(n_terminal_rows == 1L) == n_distinct(baseline_chars$patient_id)
    ),
  
  #Check time_to_assignment is within expected window from t_0: condition sohuld be T
  vary_chars %>%
    distinct(patient_id, t_0, time_to_assignment) %>%
    mutate(
      hours_t0_to_event = as.numeric(difftime(time_to_assignment, t_0, units = "hours"))
    ) %>%
    summarise(
      check = "time_to_assignment_within_0_24h",
      pass =
        sum(hours_t0_to_event < 0, na.rm = TRUE) == 0 &
        sum(hours_t0_to_event > 24, na.rm = TRUE) == 0
    ),
  
  #Check Patient count retained matches original varying table
  vary_chars %>%
    summarise(
      check = "patient_count_retained",
      pass =
        n_distinct(patient_id) == n_distinct(baseline_chars$patient_id) &
        n_distinct(patient_id) == n_distinct(trt_assignment_time$patient_id)
    ),
  
  #Checking if there are any time block from after treatment assignment
  vary_chars %>%
    summarise(
      check = "no_time_block_after_assignment",
      pass = sum(!is.na(time_to_assignment) & block_start > time_to_assignment) == 0
    ),
  
  #Check any patient with a missing time block before treatment assignment: condition should be TRUE
  vary_chars %>%
    group_by(patient_id) %>%
    summarise(
      expected_blocks = list(seq(min(time_block), max(time_block), by = 2)),
      observed_blocks = list(sort(unique(time_block))),
      missing_blocks = list(setdiff(expected_blocks[[1]], observed_blocks[[1]])),
      n_missing_blocks = length(missing_blocks[[1]]),
      .groups = "drop"
    ) %>%
    summarise(
      check = "no_missing_time_blocks",
      pass = sum(n_missing_blocks > 0) == 0
    )
)

if (any(!qc_summary$pass)) {
  print(qc_summary)
  stop("QC failed: at least one check is FALSE.")
} else {
  print(qc_summary)
  message("All QC checks passed.")
}

check_vary(
  data = vary_chars ,
  ids = sample_ids,
  print_n = 100,
  exclude_cols = setdiff(names(vary_chars), keep_cols),
  group = c("patient_id")
) %>% 
  View()

check_vary(data = clif_adt %>% select(patient_id, in_dttm, location_category) %>% 
             collect() %>% 
             inner_join(baseline_chars %>% 
                          select(patient_id, t_0), 
                        by = "patient_id"),
           ids = sample_ids,
           print_n = 100, 
           exclude_cols = NULL, 
           group = c("patient_id")) 

  

###Outcome
##A.Primary outcome (hospice and death by 28 days) and first secondary outcome (hospice and death by 60 days)
#Uniqueness of patient id: 
outcomes_chars %>%
  count(patient_id) %>% 
  count(n, name = "patient_count_by_unique_id")

#Total death count:
outcomes_chars %>% 
  summarise(death_or_hospice_28d_n = sum(death_or_hospice_by_28d == 1, na.rm = TRUE),
            death_or_hospice_60d_n = sum(death_or_hospice_by_60d == 1, na.rm = TRUE))

#Temporal logic checks:
outcomes_chars %>% 
  summarise(#if both 28 and 60 day indicators == 1, then death days are the same
            different_death_dttm_with_identical_28_60_day_indicator = sum(death_or_hospice_by_28d == 1 & death_or_hospice_by_60d == 1 & death_dttm_by_60d != death_dttm_by_60d, na.rm = T),
            #if 28 indicator == 1, then death days <= 28 and > 0
            death_dttm_outside_range_with_28_day_indicator = sum(death_or_hospice_by_28d == 1 & (!death_dttm_by_28d > t_0 | !death_dttm_by_28d <= t_0 + days(28)), na.rm = T),
            #if 60 indicator == 1 and 28 day indicator != 1, then death days > 28 and < 60
            death_dttm_outside_range_with_60_day_indicator = sum(death_or_hospice_by_60d == 1 & death_or_hospice_by_28d != 1 & (!death_dttm_by_60d > t_0 + days(28) | !death_dttm_by_60d <= t_0 + days(60)), na.rm = T)
            )

##B.Secondary outcome (resp. free days by day 28)
outcomes_chars %>% 
  count(resp_support_free_days_28 , name = "total patients") %>% 
  arrange(desc(`total patients`)) %>% 
  print(n=100)

##C.Secondary outcome (escalation of care within 7 days)
qa1 <- outcomes_chars %>% 
  inner_join(baseline_chars %>% select(patient_id, randomization), 
             by = "patient_id")

##D.Escalation logic matches randomization
table(qa1$randomization, qa1$escalation_status_7d, useNA = "ifany")



###Missingness
baseline_missingness <- baseline_chars %>%
  summarise(across(
    everything(),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    percent_missing = round(100 * n_missing / nrow(baseline_chars), 1)
  ) %>%
  arrange(desc(percent_missing))

baseline_missingness

vary_missingness <- vary_chars %>%
  summarise(across(
    everything(),
    ~ sum(is.na(.x))
  )) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    percent_missing = round(100 * n_missing / nrow(vary_chars), 1)
  ) %>%
  arrange(desc(percent_missing))

vary_missingness

assembled_df %>%
  semi_join(random_patient_ids, by = "patient_id") %>%
  arrange(patient_id, t_0, time_block, block_start, treatment, treatment_transition_path, time_to_assignment, time_to_censor_from_assignment, death, discharge_death_hospice_dttm_by_28d) %>%
  print(n = Inf)

patient_level_summary <- assembled_df %>%
  group_by(patient_id, hospital_block_id) %>%
  summarise(
    t_0 = first(t_0),
    treatment_indicator = first(treatment_indicator),
    treatment_transition_path = first(treatment_transition_path),
    time_to_censor_from_treatment = first(time_to_censor_from_treatment),
    death_or_hospice_any = as.integer(any(death_or_hospice == 1L, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    treatment_indicator_num = as.integer(as.character(treatment_indicator)),
    treatment_group = case_when(
      treatment_indicator_num == 0L ~ "ICU",
      treatment_indicator_num == 1L ~ "Stepdown",
      is.na(treatment_indicator_num) ~ "Neither"
    ),
    neither_status = case_when(
      treatment_group != "Neither" ~ NA_character_,
      !is.na(time_to_censor_from_treatment) &
        time_to_censor_from_treatment >= t_0 + hours(24) ~ "Reached t_0 + 24h",
      !is.na(time_to_censor_from_treatment) &
        time_to_censor_from_treatment < t_0 + hours(24) ~ "Left before t_0 + 24h",
      TRUE ~ "Missing censor time"
    ),
    assigned_status = case_when(
      treatment_group %in% c("ICU", "Stepdown") & death_or_hospice_any == 1L ~ "Death/hospice",
      treatment_group %in% c("ICU", "Stepdown") & death_or_hospice_any == 0L ~ "No death/hospice",
      TRUE ~ NA_character_
    )
  )

total_by_treatment_group <- patient_level_summary %>%
  count(treatment_group, name = "n_patients")

total_by_treatment_group %>%
  janitor::adorn_totals(where = "row")

neither_by_transition_path <- patient_level_summary %>%
  filter(treatment_group == "Neither") %>%
  count(treatment_transition_path, neither_status, name = "n_patients") %>%
  tidyr::pivot_wider(
    names_from = neither_status,
    values_from = n_patients,
    values_fill = 0
  ) %>%
  arrange(treatment_transition_path)

neither_by_transition_path %>%
  janitor::adorn_totals(where = c("row", "col"))

assigned_by_transition_path <- patient_level_summary %>%
  filter(treatment_group %in% c("ICU", "Stepdown")) %>%
  count(treatment_group, treatment_transition_path, assigned_status, name = "n_patients") %>%
  tidyr::pivot_wider(
    names_from = assigned_status,
    values_from = n_patients,
    values_fill = 0
  ) %>%
  arrange(treatment_group, treatment_transition_path)

assigned_by_transition_path_icu <- assigned_by_transition_path %>%
  filter(treatment_group == "ICU") %>%
  select(-treatment_group) %>%
  janitor::adorn_totals(where = c("row", "col"))

assigned_by_transition_path_icu

assigned_by_transition_path_stepdown <- assigned_by_transition_path %>%
  filter(treatment_group == "Stepdown") %>%
  select(-treatment_group) %>%
  janitor::adorn_totals(where = c("row", "col"))

assigned_by_transition_path_stepdown

## Full QC checks for assembled_df outcome/treatment timing
## Uses:
## - discharge_death_hospice_dttm_by_28d
## - time_to_censor_from_treatment
## - death_or_hospice_dttm_by_28d
## - death_or_hospice
## - treatment_indicator
## Does NOT use time_to_assignment or time_to_event

patient_level_qc <- assembled_df %>%
  group_by(patient_id, hospital_block_id) %>%
  summarise(
    n_rows = n(),
    t_0 = first(t_0),
    
    treatment_indicator = first(treatment_indicator),
    treatment_transition_path = first(treatment_transition_path),
    
    discharge_death_hospice_dttm_by_28d = first(discharge_death_hospice_dttm_by_28d),
    time_to_censor_from_treatment = first(time_to_censor_from_treatment),
    death_or_hospice_dttm_by_28d = first(death_or_hospice_dttm_by_28d),
    
    n_death_or_hospice_rows = sum(death_or_hospice == 1L, na.rm = TRUE),
    any_death_or_hospice = any(death_or_hospice == 1L, na.rm = TRUE),
    
    min_block_start = min(block_start, na.rm = TRUE),
    max_block_end = max(block_end, na.rm = TRUE),
    
    final_time_block = max(time_block, na.rm = TRUE),
    final_block_start = block_start[which.max(time_block)],
    final_block_end = block_end[which.max(time_block)],
    
    .groups = "drop"
  ) %>%
  mutate(
    stop_type = case_when(
      !is.na(discharge_death_hospice_dttm_by_28d) &
        is.na(time_to_censor_from_treatment) ~ "Group 1: assigned/discharge-death-hospice stop",
      
      is.na(discharge_death_hospice_dttm_by_28d) &
        !is.na(time_to_censor_from_treatment) ~ "Group 2: censored from treatment",
      
      !is.na(discharge_death_hospice_dttm_by_28d) &
        !is.na(time_to_censor_from_treatment) ~ "Invalid: both stop times present",
      
      is.na(discharge_death_hospice_dttm_by_28d) &
        is.na(time_to_censor_from_treatment) ~ "Invalid: neither stop time present"
    ),
    stop_time = case_when(
      stop_type == "Group 1: assigned/discharge-death-hospice stop" ~ discharge_death_hospice_dttm_by_28d,
      stop_type == "Group 2: censored from treatment" ~ time_to_censor_from_treatment,
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    )
  )


qc_stop_time_exclusivity <- patient_level_qc %>%
  summarise(
    both_stop_times_present = sum(
      !is.na(discharge_death_hospice_dttm_by_28d) &
        !is.na(time_to_censor_from_treatment)
    ),
    neither_stop_time_present = sum(
      is.na(discharge_death_hospice_dttm_by_28d) &
        is.na(time_to_censor_from_treatment)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_group1 <- patient_level_qc %>%
  filter(stop_type == "Group 1: assigned/discharge-death-hospice stop") %>%
  summarise(
    group1_treatment_indicator_missing = sum(is.na(treatment_indicator)),
    
    group1_treatment_indicator_not_0_or_1 = sum(
      !is.na(treatment_indicator) &
        !treatment_indicator %in% c(0L, 1L)
    ),
    
    group1_final_row_does_not_contain_stop_time = sum(
      discharge_death_hospice_dttm_by_28d <= final_block_start |
        discharge_death_hospice_dttm_by_28d > final_block_end
    ),
    
    group1_event_time_present_but_not_exactly_one_death_row = sum(
      !is.na(death_or_hospice_dttm_by_28d) &
        n_death_or_hospice_rows != 1
    ),
    
    group1_no_event_time_but_has_death_row = sum(
      is.na(death_or_hospice_dttm_by_28d) &
        n_death_or_hospice_rows != 0
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_group1_event_interval <- assembled_df %>%
  filter(!is.na(discharge_death_hospice_dttm_by_28d)) %>%
  summarise(
    check = "group1_death_or_hospice_row_contains_event_time",
    n_fail = sum(
      death_or_hospice == 1L &
        (
          is.na(death_or_hospice_dttm_by_28d) |
            death_or_hospice_dttm_by_28d <= block_start |
            death_or_hospice_dttm_by_28d > block_end
        ),
      na.rm = TRUE
    ),
    pass = n_fail == 0
  )


qc_group1_event_matching_row <- assembled_df %>%
  filter(!is.na(death_or_hospice_dttm_by_28d)) %>%
  group_by(patient_id, hospital_block_id) %>%
  summarise(
    n_matching_event_interval_rows = sum(
      death_or_hospice_dttm_by_28d > block_start &
        death_or_hospice_dttm_by_28d <= block_end,
      na.rm = TRUE
    ),
    n_death_or_hospice_rows = sum(death_or_hospice == 1L, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(
    event_time_not_in_exactly_one_row = sum(n_matching_event_interval_rows != 1),
    event_time_row_indicator_not_exactly_one = sum(n_death_or_hospice_rows != 1)
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_group2 <- patient_level_qc %>%
  filter(stop_type == "Group 2: censored from treatment") %>%
  summarise(
    group2_treatment_indicator_not_missing = sum(!is.na(treatment_indicator)),
    
    group2_has_death_or_hospice_event_time = sum(!is.na(death_or_hospice_dttm_by_28d)),
    
    group2_has_death_or_hospice_row = sum(n_death_or_hospice_rows > 0),
    
    group2_final_row_does_not_contain_censor_time = sum(
      time_to_censor_from_treatment <= final_block_start |
        time_to_censor_from_treatment > final_block_end
    ),
    
    group2_censor_before_t0 = sum(
      time_to_censor_from_treatment < t_0
    ),
    
    group2_censor_after_t0_plus_24h = sum(
      time_to_censor_from_treatment > t_0 + hours(24)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_no_rows_after_stop <- assembled_df %>%
  mutate(
    stop_time = case_when(
      !is.na(discharge_death_hospice_dttm_by_28d) &
        is.na(time_to_censor_from_treatment) ~ discharge_death_hospice_dttm_by_28d,
      
      is.na(discharge_death_hospice_dttm_by_28d) &
        !is.na(time_to_censor_from_treatment) ~ time_to_censor_from_treatment,
      
      TRUE ~ as.POSIXct(NA, tz = "UTC")
    )
  ) %>%
  summarise(
    check = "no_rows_starting_at_or_after_stop_time",
    n_fail = sum(
      !is.na(stop_time) &
        block_start >= stop_time
    ),
    pass = n_fail == 0
  )


qc_block_structure <- assembled_df %>%
  arrange(patient_id, hospital_block_id, time_block) %>%
  group_by(patient_id, hospital_block_id) %>%
  mutate(
    row_num = row_number(),
    prior_block_end = lag(block_end),
    prior_time_block = lag(time_block)
  ) %>%
  summarise(
    duplicate_time_blocks = n() - n_distinct(time_block),
    
    first_block_start_not_t0 = sum(
      row_num == 1 & block_start != t_0,
      na.rm = TRUE
    ),
    
    block_end_not_after_block_start = sum(
      block_end <= block_start,
      na.rm = TRUE
    ),
    
    non_contiguous_blocks = sum(
      row_num > 1 & block_start != prior_block_end,
      na.rm = TRUE
    ),
    
    non_increasing_time_block = sum(
      row_num > 1 & time_block <= prior_time_block,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  summarise(
    duplicate_time_blocks = sum(duplicate_time_blocks),
    first_block_start_not_t0 = sum(first_block_start_not_t0),
    block_end_not_after_block_start = sum(block_end_not_after_block_start),
    non_contiguous_blocks = sum(non_contiguous_blocks),
    non_increasing_time_block = sum(non_increasing_time_block)
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_treatment_path <- patient_level_qc %>%
  mutate(path_lower = str_to_lower(treatment_transition_path)) %>%
  summarise(
    icu_indicator_without_icu_in_path = sum(
      treatment_indicator == 0L &
        !is.na(path_lower) &
        !str_detect(path_lower, "icu")
    ),
    
    stepdown_indicator_without_stepdown_in_path = sum(
      treatment_indicator == 1L &
        !is.na(path_lower) &
        !str_detect(path_lower, "stepdown")
    ),
    
    neither_indicator_with_icu_or_stepdown_first_non_ed = sum(
      is.na(treatment_indicator) &
        str_detect(path_lower, "^ed -> (icu|stepdown)")
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "check",
    values_to = "n_fail"
  ) %>%
  mutate(pass = n_fail == 0)


qc_outcome_count_match <- tibble::tibble(
  source = c("outcomes_chars", "assembled_df"),
  n_death_or_hospice = c(
    outcomes_chars %>%
      filter(death_or_hospice_by_28d == 1L) %>%
      distinct(patient_id, hospital_block_id) %>%
      nrow(),
    
    assembled_df %>%
      group_by(patient_id, hospital_block_id) %>%
      summarise(
        any_death_or_hospice = any(death_or_hospice == 1L, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(any_death_or_hospice) %>%
      nrow()
  )
)


qc_all <- bind_rows(
  qc_stop_time_exclusivity,
  qc_group1,
  qc_group1_event_interval,
  qc_group1_event_matching_row,
  qc_group2,
  qc_no_rows_after_stop,
  qc_block_structure,
  qc_treatment_path
)


qc_failures <- qc_all %>%
  filter(!pass)


stop_type_counts <- patient_level_qc %>%
  count(stop_type, name = "n_patients") %>%
  janitor::adorn_totals(where = "row")


treatment_indicator_counts <- patient_level_qc %>%
  count(treatment_indicator, name = "n_patients") %>%
  janitor::adorn_totals(where = "row")


group1_outcome_counts <- patient_level_qc %>%
  filter(stop_type == "Group 1: assigned/discharge-death-hospice stop") %>%
  mutate(
    outcome_status = if_else(
      n_death_or_hospice_rows == 1L,
      "Death/hospice",
      "No death/hospice"
    )
  ) %>%
  count(treatment_indicator, outcome_status, name = "n_patients") %>%
  tidyr::pivot_wider(
    names_from = outcome_status,
    values_from = n_patients,
    values_fill = 0
  ) %>%
  janitor::adorn_totals(where = c("row", "col"))


group2_censor_counts <- patient_level_qc %>%
  filter(stop_type == "Group 2: censored from treatment") %>%
  mutate(
    censor_status = case_when(
      time_to_censor_from_treatment >= t_0 + hours(24) ~ "Reached t_0 + 24h",
      time_to_censor_from_treatment < t_0 + hours(24) ~ "Left before t_0 + 24h"
    )
  ) %>%
  count(treatment_transition_path, censor_status, name = "n_patients") %>%
  tidyr::pivot_wider(
    names_from = censor_status,
    values_from = n_patients,
    values_fill = 0
  ) %>%
  janitor::adorn_totals(where = c("row", "col"))


## Print main QC outputs
qc_all
qc_failures
qc_outcome_count_match
stop_type_counts
treatment_indicator_counts
group1_outcome_counts
group2_censor_counts

 
