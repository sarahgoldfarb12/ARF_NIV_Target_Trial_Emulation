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
exclude_cols <- setdiff(exclude_cols, "device_category")
exclude_cols <- c(exclude_cols, "t_0", "hospital_block_id")

set.seed(1)

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