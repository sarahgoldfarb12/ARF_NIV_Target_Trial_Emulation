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


###1.Randomization and censoring 

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

##E.Restricting vary_chars to before each individual's time to treatment initiation or censoring
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
        sum(is.na(time_to_event)) == 0 &
        sum(is.na(censor)) == 0 &
        sum(is.na(at_risk)) == 0
    ),
  
  #Check any row after treatment/censor time: should be 0
  vary_chars %>%
    summarise(
      check = "no_rows_after_event",
      pass = sum(block_start > time_to_event, na.rm = TRUE) == 0
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
      n_terminal_rows = sum(!is.na(treatment) | censor == 1L, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    summarise(
      check = "one_terminal_row_per_patient",
      pass =
        n() == n_distinct(baseline_chars$patient_id) &
        sum(n_terminal_rows == 1L) == n_distinct(baseline_chars$patient_id)
    ),
  
  #Check time_to_event is within expected window from t_0: condition sohuld be T
  vary_chars %>%
    distinct(patient_id, t_0, time_to_event) %>%
    mutate(
      hours_t0_to_event = as.numeric(difftime(time_to_event, t_0, units = "hours"))
    ) %>%
    summarise(
      check = "time_to_event_within_0_24h",
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
      pass = sum(!is.na(time_to_event) & block_start > time_to_event) == 0
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

#Check a sample of the data
keep_cols <- c(
  "patient_id",
  "t_0",
  "time_block",
  "block_start",
  "block_end",
  "treatment",
  "time_to_event",
  "censor",
  "at_risk"
)

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
    c(year, age, albumin_baseline:elixhauser_count),
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
    code_status_category:naloxone,
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


