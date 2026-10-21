library(dplyr)
library(readr)

# ----------------- Import saved assembled_df
project_root <- rprojroot::find_root(rprojroot::has_dir("config"))
config <- yaml::read_yaml(file.path(project_root, "config", "config.yaml"))

assembled_df <- read_csv(
  file.path(config$project_location, "private_tables", "assembled_df.csv"),
  col_types = cols(.default = col_character())
)

# ----------------- Keep time_block == 2: ONLY checks missingness in block 2
block2 <- assembled_df %>%
  filter(time_block == "2")

#Check Patient-level percentages require one row per patient.
n_unique_pts = length(unique(block2$patient_id))
stopifnot(n_unique_pts == dim(block2)[1])

#Original covariates only
varying_cols <- names(select(block2, unit_location:naloxone))
fixed_cols <- names(select(block2, year:elixhauser_count))

# ----------------- Variable types
binary_cols <- c(
  "antihypertensive_drip", "naloxone",
  "pulmonary_vasc", "renal_failure",
  "liver_disease", "metastatic_cancer"
)

categorical_cols <- c(
  "unit_location", "code_status_category", "device_category",
  "sofa_brain", "sofa_liver", "sofa_coagulation", "sofa_kidney",
  "discharge_category",
  "year" 
)

# ----------------- Identify missing values
# Includes NA, NaN, blanks, and "not measured".
# Zero and FALSE are not missing.
missing_flags <- block2 %>%
  select(all_of(c(fixed_cols, varying_cols))) %>%
  mutate(across(
    everything(),
    ~ is.na(.x) |
      tolower(trimws(.x)) %in% c("", "na", "n/a", "nan", "<na>", "null",
                                 "na_character_", "na_real_", "na_integer_", "na_complex_",
                                 "missing", "not measured", "not recorded",
                                 "not available", "unavailable")
  ))

# ----------------- Missingness by covariate
missingness_block2 <- tibble(
  variable = names(missing_flags),
  
  covariate_type = if_else(
    variable %in% fixed_cols, "fixed", "varying"
  ),
  
  variable_type = case_when(
    variable %in% binary_cols ~ "binary",
    variable %in% categorical_cols ~ "categorical",
    variable %in% c("admission_dttm", "discharge_dttm") ~ "datetime",
    variable == "elixhauser_count" ~ "count",
    variable == "elixhauser_index" ~ "numeric score",
    TRUE ~ "continuous"
  ),
  
  n_patients = nrow(block2),
  missing_n = colSums(missing_flags),
  missing_pct = round(100 * colMeans(missing_flags), 2)
) %>%
  arrange(desc(covariate_type), desc(missing_pct), variable)

print(missingness_block2, n = Inf)

# ----------------- Patient count by number of missing covariates
missingness_per_patient <- tibble(
  n_missing_vars = rowSums(missing_flags)
) %>%
  count(n_missing_vars, name = "n_patients") %>%
  tidyr::complete(n_missing_vars = 0:10, fill = list(n_patients = 0)) %>%
  arrange(n_missing_vars) %>% 
  mutate(pct_patients = round(100*n_patients/dim(block2)[1], 2))

print(missingness_per_patient, n = Inf)

# ----------------- Patient count by number of missing covariates, excluding top 5 sources (vars) of missingness
top5_missing_vars <- names(sort(colSums(missing_flags), decreasing = TRUE))[1:5]

missingness_per_patient_excluding_top5_missing_vars <- tibble(
  n_missing_vars = rowSums(
    select(missing_flags, -all_of(top5_missing_vars))
  )
) %>%
  count(n_missing_vars, name = "n_patients") %>%
  tidyr::complete(n_missing_vars = 0:10, fill = list(n_patients = 0)) %>%
  arrange(n_missing_vars) %>% 
  mutate(pct_patients = round(100*n_patients/dim(missing_flags)[1], 2))

print(top5_missing_vars)
print(missingness_per_patient_excluding_top5_missing_vars, n = Inf)
