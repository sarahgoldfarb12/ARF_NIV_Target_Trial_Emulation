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

exclude_cols <- names(vary_chars)[7:length(names(vary_chars))]
exclude_cols <- setdiff(exclude_cols, c("device_category", "avg_pf", "avg_sf"))
exclude_cols <- c(exclude_cols, "t_0", "hospital_block_id")

set.seed(12)

sample_ids <- vary_chars %>%
  distinct(patient_id) %>%
  slice_sample(n = 10) %>%
  pull(patient_id)

check_vary(
  data = vary_chars,
  ids = sample_ids,
  id_col = "patient_id",
  print_n = 100,
  exclude_cols = exclude_cols
)

###Randomization Indicator
table(baseline_chars $transition_path, baseline_chars$randomization)

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

