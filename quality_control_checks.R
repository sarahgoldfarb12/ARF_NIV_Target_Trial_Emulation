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


###1.Randomization 

exclude_cols <- names(vary_chars)[8:length(names(vary_chars))]

##A. manually compare extracted treatment transition with actual treatment trajectory from clif-adt
check_vary(data = clif_adt %>% select(patient_id, in_dttm, location_category) %>% collect() %>% inner_join(baseline_chars %>% 
                                                                                                             select(patient_id, t_0),
                                                                                                           by = "patient_id"
),
          ids = sample_ids,
          print_n = 100, 
          exclude_cols = NULL, 
          group = c("patient_id")) %>% View()

check_vary(data = transition_baseline,
            ids = sample_ids,
            print_n = 100,
            exclude_cols = NULL,
            group = c("patient_id"))


##B. check treatment_transition_path and randomization indicator consistency
qc_t1 <- table(transition_baseline$treatment_transition_path, transition_baseline$randomization, useNA = "ifany") %>% as.data.frame.matrix()
bind_rows(qc_t1, 
          tibble(`0` = colSums(qc_t1)[1],
                 `1` = colSums(qc_t1)[2]))

##C. check trt transition, trt indicator, and trt assignment time are consistent
check_vary(data = clif_adt %>% select(patient_id, in_dttm, location_category) %>% 
             collect() %>% 
             inner_join(baseline_chars %>% 
                          select(patient_id, t_0), 
                                 by = "patient_id"),
                                      ids = sample_ids,
                                      print_n = 100, 
                                      exclude_cols = NULL, 
                                      group = c("patient_id")) %>% View()

check_vary(data = trt_assignment_time %>% 
             left_join(transition_baseline,
                       by = "patient_id"),
           ids = sample_ids,
           print_n = 100,
           exclude_cols = NULL,
           group = c("patient_id"))

table(baseline_chars$treatment_transition_path, baseline_chars$treatment_assignment, useNA = "ifany") 
table(baseline_chars$treatment_transition_path, baseline_chars$treatment_assignment, useNA = "ifany") %>% colSums()
table(baseline_chars$treatment_transition_path, baseline_chars$treatment_assignment, useNA = "ifany") %>% colSums() %>% sum()


###Outcome
##Primary outcome (hospice and death by 28 days) and first secondary outcome (hospice and death by 60 days)
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

##Secondary outcome (resp. free days by day 28)
outcomes_chars %>% 
  count(room_air_days_28, name = "total patients") %>% 
  arrange(desc(`total patients`)) %>% 
  print(n=100)

##Secondary outcome (escalation of care within 7 days)
qa1 <- outcomes_chars %>% 
  inner_join(baseline_chars %>% select(patient_id, randomization), 
             by = "patient_id")

#Escalation logic matches randomization
table(qa1$randomization, qa1$escalation_status_7d, useNA = "ifany")


