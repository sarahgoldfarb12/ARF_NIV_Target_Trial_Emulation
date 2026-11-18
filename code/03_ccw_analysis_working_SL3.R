#' ---
#' title: "IMC_ICU_ARF_TargetTrialEmulation_Analysis"
#' author: "Chad Hochberg, Lizi Shao, Sarah Goldfarb"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---
#' 
 ## ----Install Packages---------------------------------------------------------
packages <- c(
                  "lubridate", 
                  "tidyverse", 
                  "dplyr",
                  "tableone", 
                  "broom", 
                  "arrow", 
                  "rvest", 
                  "readr", 
                  "data.table", 
                  "collapse", 
                  "yaml",
                  "this.path",
                  "comorbidity",
                  "WeightIt",
                  "tidytable",
                  "future",
                  "ranger"
                  )
                  #"SuperLearner",
                  #"nnet")
    
    install_if_missing <- function(package) {
      if (!require(package, character.only = TRUE)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = TRUE)
      }
    }
    
    sapply(packages, install_if_missing)
    rm(packages, install_if_missing)

  if (!requireNamespace("tlverse", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    } 
    devtools::install_github("tlverse/tlverse")
   }

  library(tlverse)
  library(sl3)
##--- 

## ----Globally Defined Variables-----------------------------------------------
MAX_FU = 24*28 # Maximum Amount of Follow-up in Hours at 28 Days

#' 
#' 
#' 
## ----Read in Config-----------------------------------------------------------
# Find project root
project_root <- normalizePath(getwd())
    
    # Read YAML config
config <- yaml::read_yaml(file.path(project_root, "config", "config.yaml"))
global_config <- yaml::read_yaml(file.path(project_root, "config", "global_config.yaml"))
    
# Assign config values to R variables
tables_location <- config$tables_location
project_location <- config$project_location
site <- config$institution
site_time_zone <- config$time_zone
file_type <- config$file_type

#Create Folders if They Do NOT already exist
if (!dir.exists(paste0(project_location, "/", site, "_project_output/graphs"))) {
      dir.create(paste0(project_location,"/", site, "_project_output/graphs"))
}
#' 
## ----Load in Data-------------------------------------------------------------
filenames <- list.files(
  paste0(project_location, "/private_tables/"),
  pattern = "\\.csv$",
  full.names = FALSE
)

raw_dfs <- list()

for (filename in filenames) {
  data <- read_csv(paste0(project_location, "/private_tables/", filename), show_col_types = FALSE)
  
  name_key <- tools::file_path_sans_ext(filename)
  
  raw_dfs[[name_key]] <- data
}

#' 
#' #Analytic Plan
#' - Primary Approach: Clone Censor Weighting Using Super Learner to Calculate Censoring Weights Based on Pre Triage Trajectory
## ----Format Assembled DF for Clone Censor Weighting---------------------------
final_cohort <- as.data.table(raw_dfs$no_pandemic_one_encounter_per_patient)

censor_df <- as.data.table(raw_dfs$assembled_df) %>%
  #Bring in Needed Data from Final Cohort 
  left_join(
    final_cohort %>% select(-c(patient_id, first_trach_dttm, t_0, first_imv_dttm, first_vasoactive_time, first_vasoactive_name, first_vasoactive_dose,
                               first_vasoactive_unit, first_crrt_mode, first_crrt_dttm))
  ) |>
  mutate(
    out_of_ed = fifelse(unit_location != 'ed' & !is.na(unit_location), 1, 0L),
    is_first_non_ed = as.integer(dplyr::row_number() == which(out_of_ed == 1)[1]),
    last_block_was_ed = fifelse(
      last(unit_location[!is.na(unit_location)]) == "ed", 1L, NA_real_),
    in_ed_at_24hours = fifelse(unit_location=='ed' & time_block>=24, 1L, 0L),
    first_icu_indication=coalesce(first_imv_dttm, first_vasoactive_dttm), #For Purposes of this Analysis Assume This is Just Vasopressors and IMV
    #Identify Patients That Develop an ICU Indication While Still in ED (IMV, Vasopressors)
    icu_indication_censor=fifelse(
      out_of_ed==0L & first_icu_indication<=block_end & 
        first_icu_indication<=time_to_treatment, 1L, NA_real_ #NOte time_to_treatment here is ED triage time
    ),
    temp=fifelse(icu_indication_censor==1L, time_block, NA_real_),
    icu_indication_censor_block=ifelse(is.finite(min(temp, na.rm = TRUE)), min(temp, na.rm = TRUE), NA_real_),
    icu_indication_in_ed=fifelse(!is.na(icu_indication_censor_block), 1L, 0L),
    .by = hospital_block_id
  ) |>
  fill(triage_location, in_ed_at_24hours, icu_indication_in_ed, last_block_was_ed, .direction=c('downup'), .by = hospital_block_id) |>
  mutate(
    #Possible Outcomes from ED: Still in ED at Hour 24, In ICU at end Hour 24, In IMC at Hour 24, Developed an ICU indication in ED, Died in ED
    triage_outcome=fcase(
      in_ed_at_24hours==1, 'In ED at 24 Hours',
      icu_indication_in_ed==1, 'ICU indication in ED',
      last_block_was_ed==1 & max(time_block)<24, 'Left ED/Hospital Within 24 Hours',
      in_ed_at_24hours==0 & triage_location=='icu', 'In ICU at 24 Hours',
      in_ed_at_24hours==0 & triage_location=='stepdown', 'In IMC at 24 Hours',
      in_ed_at_24hours==0 & triage_location=='ward', 'In Ward at 24 Hours'),
    #Identify Which Block This Outcome Occurred In
    triage_block=fcase(
      in_ed_at_24hours==1 & icu_indication_in_ed==0, 24,
      icu_indication_in_ed==1, icu_indication_censor_block,
      is_first_non_ed==1, time_block),
      .by = hospital_block_id
    ) |>
  fill(triage_location, in_ed_at_24hours, triage_outcome, triage_block, .direction=c('downup'), .by = hospital_block_id) |>
  select(-temp) 

clone0 <- copy(censor_df)[, icu_arm:=0]
clone1 <- copy(censor_df)[, icu_arm:=1]

# Check the number of clones in each arm
temp0 <- length(clone0$time_block[which(clone0$time_block==2)])
temp1 <- length(clone1$time_block[which(clone1$time_block==2)])

n_events0 <- clone0 %>%
  summarise(n_events = sum(death_or_hospice_by_28d == 1, na.rm = TRUE), 
            .by = hospital_block_id)
n_events1 <- clone1 %>%
  summarise(n_events = sum(death_or_hospice_by_28d == 1, na.rm = TRUE), 
            .by = hospital_block_id)


#Keep Track of Clone Data for Diagnostics
clone_diagnostics <- data.frame(
  clone = c('clone0', 'clone1'),
  description = c('Not ICU', 'ICU'),
  n_clones_prior_censoring = c(temp0, temp1),
  n_events_per_arm = c(sum(n_events0$n_events, na.rm=T), sum(n_events1$n_events, na.rm=T))
)

#Censor Clones When They Are No Longer Adherent to The Strategies of:
#--clone0=imc by end of 24 hours OR if develops an ICU indication prior to triage censor
clone0 <- clone0 |>
  #Censor if in ICU before hour 24
  mutate(
    first_censored=fcase(
      time_block>=24 & in_ed_at_24hours==1, 1L,
      (icu_arm==0 & triage_location!='stepdown' & is_first_non_ed==1), 1L,
      default = 0L),
    censored_reason=fcase(
      first_censored==1L & in_ed_at_24hours==1, 'Still in ED at 24 Hours',
      first_censored==1L & out_of_ed==1 & triage_location=='icu', 'To ICU',
      first_censored==1L & out_of_ed==1 & triage_location=='ward', 'To Ward'
      ),
    censored=cummax(first_censored), 
    temp_censor=cumsum(censored),
    .by = hospital_block_id
  ) |>
  filter(
    temp_censor<=1, 
    .by = hospital_block_id
  ) |>
  select(-temp_censor) |>
  #Create a Block Status Indicating of Clone is Adherent
  fill(censored_reason, .direction='down', .by=hospital_block_id) |>
  mutate(
    block_status=fcase(
      time_block>=icu_indication_censor_block, 'icu_indication_in_ed',
      time_block>=triage_block & !is.na(censored_reason), censored_reason,
      time_block<=triage_block & is.na(censored_reason), 'following_strategy')) |>
  fill(block_status, .direction='down', .by=hospital_block_id) |>
  #Meeting an ICU indication is not a censoring event
  mutate(
    censored=fifelse(
      block_status=='icu_indication_in_ed', 0L, censored
    )
  )
  
#--clone1=icu by end of 24 hours OR if develops an ICU indication prior to triage censor
clone1 <- clone1 |>
  #Censor if in stepdown before hour 24
  mutate(
    first_censored=fcase(
      time_block>=24 & in_ed_at_24hours==1, 1L,
      (icu_arm==1 & triage_location!='icu' & is_first_non_ed==1), 1L,
      default = 0L),
    censored_reason=fcase(
      first_censored==1L & in_ed_at_24hours==1, 'Still in ED at 24 Hours',
      first_censored==1L & out_of_ed==1 & triage_location=='stepdown', 'To IMC',
      first_censored==1L & out_of_ed==1 & triage_location=='ward', 'To Ward'),
    censored=cummax(first_censored), 
    temp_censor=cumsum(censored),
    .by = hospital_block_id
  ) |>
 filter(
    temp_censor<=1, 
    .by = hospital_block_id
  ) |>
  select(-temp_censor) |>
  #Create a Block Status Indicating of Clone is Adherent
  fill(censored_reason, .direction='down', .by=hospital_block_id) |>
  mutate(
    block_status=fcase(
      time_block>=icu_indication_censor_block, 'icu_indication_in_ed',
      time_block>=triage_block & !is.na(censored_reason), censored_reason,
      time_block<=triage_block & is.na(censored_reason), 'following_strategy')) |>
  fill(block_status, .direction='down', .by=hospital_block_id)  |>
#Meeting an ICU indication is not a censoring event
  mutate(
    censored=fifelse(
      block_status=='icu_indication_in_ed', 0L, censored
    )
  )

#Deal With THose that Leave ED Within 24 Hours
both.arms <- rbind(clone0,clone1) |>
  mutate(
    block_status=fcase(
      !is.na(block_status), block_status,
      last_block_was_ed==1 & time_block==max(time_block, na.rm=T), 'Left ED/Hospital Within 24 Hours', 
      default='following_strategy'),
    censored=fifelse(
      last_block_was_ed==1 & time_block==max(time_block, na.rm=T), 1L, censored),
    .by=hospital_block_id
  ) 

#' 
## ----Plot Censoring-----------------------------------------------------------
table(both.arms$icu_arm, both.arms$block_status, useNA="always")

# Prep data for plotting
plot.counts <- both.arms |>
  dplyr::count(time_block, block_status, icu_arm, name = "total_count")
plot.counts$arm_name <- ifelse(plot.counts$icu_arm==0, "Triage to IMC", "Triage to ICU")

# Create plot
cens_plot <- ggplot(data=plot.counts %>% filter(time_block<=24), aes(x=time_block, y=total_count, fill = as.factor(block_status))) +
  geom_bar(position="stack", stat="identity") +
  xlab("Hours from Time Zero") + 
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)) +
  ylab("Number of individuals") + 
  scale_y_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000)) + 
  scale_fill_manual(values = c("#878787", "#14B5D0", "#E8A33D", "#800000", "royalblue", "purple", "black")) + # set colors
  theme_minimal() +
  theme(axis.text = element_text(size=10), 
        legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())+
  facet_wrap(~arm_name, scales = "free")
cens_plot

ggsave(
  filename=paste0(site,'_clone_censoring_counts.pdf'),
  plot=last_plot(),
  device=cairo_pdf,
  width = 8,
  height = 5,
  units ='in',
  path = paste0(project_location, site, '_project_output/graphs/'))
ggsave(
  filename=paste0(site, '_clone_censoring_counts.svg'),
  plot=last_plot(),
  device=cairo_ps,
  width = 8,
  height = 5,
  units ='in',
  path = paste0(project_location, site, '_project_output/graphs/'))


#----Generate Table and Plot for General Triage Outcomes
#First Create Temporary Data SEt with Time blocks 2-24 with censoring outcomes
triage_outcome_df <- censor_df |>
  mutate(
    block_status=fcase(
      triage_block==time_block, triage_outcome,
      time_block<triage_block, 'In ED Awaiting Triage Outcome',
      time_block=max(time_block) & triage_outcome=='Left ED/Hospital Within 24 Hours', triage_outcome),
    .by=hospital_block_id) |>
  filter(time_block<=triage_block | is.na(triage_block))

triage_outcome_plot <- triage_outcome_df |>
  dplyr::count(time_block, block_status, name = "total_count")

triage_plot <- ggplot(data=triage_outcome_plot %>% filter(time_block<=24), aes(x=time_block, y=total_count, fill =   as.factor(block_status))) +
  geom_bar(position="stack", stat="identity") +
  xlab("Hours from Time Zero") + 
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)) +
  ylab("Number of individuals") + 
  scale_y_continuous(breaks=seq(0, 4000, by=500)) + 
  scale_fill_viridis_d() + # set colors
  theme_minimal() +
  theme(axis.text = element_text(size=10), 
        legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(
  filename=paste0(site,'_triage_outcome_counts.pdf'),
  plot=last_plot(),
  device=cairo_pdf,
  width = 8,
  height = 5,
  units ='in',
  path = paste0(project_location, site, '_project_output/graphs/'))
ggsave(
  filename=paste0(site, '_triage_outcome_counts.svg'),
  plot=last_plot(),
  device=cairo_ps,
  width = 8,
  height = 5,
  units ='in',
  path = paste0(project_location, site, '_project_output/graphs/'))

#' 
#' 
#' 
## ----Censoring Mechanisms as a Table------------------------------------------
#Create Table for Triage Outcomes (this will define study population)
to_tab <- c(
  'age',
  'first_hospital_id',
  #'triage_location_raw',
  'first_niv_device',
  'code_status_t0',
  'triage_block',
  'year',
  'discharge_category',
  'ed_hours'
)
factor_vars <- c(
  'first_hospital_id',
  #'triage_location_raw',
  'first_niv_device',
  'code_status_t0',
  'triage_block',
  'year',
  'discharge_category'
)
nonnorm <- c(
  'age',
  'ed_hours'
)

#Data to Enter are 1 Row Per Patient
tab_dat <- censor_df |>
  filter(dplyr::row_number()==1, .by=hospital_block_id)

tab_triage <- CreateTableOne(vars=to_tab, factorVars = factor_vars, data=tab_dat, strata = 'triage_outcome', addOverall=T)
tab_excel <- print(tab_triage, printToggle = FALSE, nonormal=nonnorm)
write.csv(tab_excel, file=paste0(project_location, '/', site, '_project_output/triage_outcome_table.csv'))

#Triage Locations (Actual Names)
to_tab <- c(
  'triage_location_raw'
)
factor_vars <- c(
  'triage_location_raw'
)

#Data to Enter are 1 Row Per Patient
tab_dat <- censor_df |>
  filter(dplyr::row_number()==1, .by=hospital_block_id)

tab_triage <- CreateTableOne(vars=to_tab, factorVars = factor_vars, data=tab_dat, strata = 'triage_outcome', addOverall=T)
tab_excel <- print(tab_triage, printToggle = FALSE, nonormal=nonnorm)
write.csv(tab_excel, file=paste0(project_location, '/', site, '_project_output/triage_locations_raw_table.csv'))

#' 
#' 
## ----Table One For Paper------------------------------------------------------


#' 
## ----Covariate and Outcome Relationship---------------------------------------


#' 
#' 
#' 
## ----Calculate Censoring Weights - 1st With Ensemble Method-------------------
#Only Need to Calculate Weights Through Triage Block
calc_weights <- both.arms |>
  filter(icu_arm==1,
         time_block<=triage_block) |>
  #These are the 5 Mutually Exclusive Classes at the End of Each Block
  mutate(
    adherence_cat=fcase(
      icu_indication_censor==1, 'icu_indication',
      unit_location=='ed', 'still_in_ed',
      unit_location=='icu', 'to_icu',
      unit_location=='stepdown', 'to_imc',
      unit_location=='ward', 'to_ward'
    )
  )
event_levels <- c(
  "icu_indication", "still_in_ed", "to_icu", "to_imc", "to_ward"
)

calc_weights$adherence_cat <- factor(
  calc_weights$adherence_cat,
  levels = event_levels
)

table(calc_weights$adherence_cat, useNA = "ifany")

#Use a One vs All Multiclass Prediction Using a SuperLearner Ensemble Method
#These are the Covariates Outlined in the DAG
predictor_names <- c(
  "age", "pulmonary_vasc_raw", "renal_failure_raw",
  "liver_disease_raw", "metastatic_cancer_raw",
  "elixhauser_index", "first_hospital_id", "year",
  "avg_temp", "avg_temp_lag1", "avg_temp_lag2",
  "time_block", "code_status_category"
)
#' 
#' 
## ----Calculate Censoring Weights - Using SuperLearner Ensemble Method---------
# 1. Enable parallel processing across available CPU cores
plan(multisession, workers = availableCores() - 1)

# 2. Define Task Predicting Probability of Disposition/Censoring
task <- make_sl3_Task(
  data = calc_weights,
  covariates = predictor_names,
  outcome = "adherence_cat", 
  outcome_type='categorical'
)

# 3. Include a preprocessing learner (Lrnr_one_hot) if you have factor/character predictors
# (Or skip Lrnr_nnet and stick to learners that handle factors natively)
lrnr_glmnet_lasso <- make_learner(Lrnr_glmnet, alpha = 1) #Lasso Regression
lrnr_glmnet_ridge <- make_learner(Lrnr_glmnet, alpha = 0) #Ridge Regression
lrnr_glmnet_enet  <- make_learner(Lrnr_glmnet, alpha = 0.5) #Elastic Net Regression
lrnr_indbinomial  <- make_learner(Lrnr_independent_binomial) #Multinomal Implementation of Logit Model
lrnr_rf           <- make_learner(Lrnr_ranger, num.trees = 500, min.node.size = 10) #Faster Random Forests
lrnr_xgb_shallow  <- make_learner(Lrnr_xgboost, nrounds = 150, max_depth = 3, eta = 0.05, tree_method = "hist") #Xgboost
lrnr_xgb_deep     <- make_learner(Lrnr_xgboost, nrounds = 150, max_depth = 6, eta = 0.05, tree_method = "hist") #Xgboost
lrnr_multinom     <- make_learner(Lrnr_nnet, size = 5,
  decay = 0.01,
  maxit = 500,
  entropy = TRUE,
  trace = FALSE) #Feed-Foward Neuroal Networks
lrnr_rpart        <- make_learner(Lrnr_rpart)
lrnr_mean         <- make_learner(Lrnr_mean) #Mean Model - Other Models Should outperform this

stack <- Stack$new(
  lrnr_glmnet_lasso,
  lrnr_glmnet_ridge,
  lrnr_glmnet_enet,
  lrnr_indbinomial,
  lrnr_rf,
  lrnr_xgb_shallow,
  lrnr_xgb_deep,
  lrnr_multinom,
  lrnr_rpart,
  lrnr_mean
)

# 2. Use Lrnr_optim (uses R's built-in optim() instead of buggy Rsolnp!)
safe_combiner <- function(alpha, X) {
  w <- exp(alpha - max(alpha)) / sum(exp(alpha - max(alpha))) #SOFT Max Function
  sl3:::metalearner_linear_multinomial(w, X)
}

metalearner <- make_learner(
  Lrnr_optim,
  loss_function = loss_loglik_multinomial,
  learner_function = safe_combiner
)

# 3. Train the Super Learner
sl <- make_learner(Lrnr_sl, learners = stack, metalearner = metalearner)
sl_fit <- sl$train(task)

# 4. Inspect balanced weights & predictions
alpha <- data.frame(
  coefficient = sl_fit$coefficients
) |>
  tibble::rownames_to_column("learner")

true_weights <- data.frame(
  true_weight = round(
    exp(alpha$coefficient - max(alpha$coefficient)) /
      sum(exp(alpha$coefficient - max(alpha$coefficient))),
    6
  )
)

sl_fit_coefficients <- dplyr::bind_cols(alpha, true_weights)
write_csv(sl_fit_coefficients, paste0(project_location, site, '_project_output/sl_fit_coefficients.csv' ))
print(sl_fit)

#Get Predictions
probs <- as.data.frame(
  sl3::unpack_predictions(sl_fit$predict(task))
)

#Save for Intermediate Processing
#Merge with IDs
prediction_ids <- calc_weights |>
   select(hospital_block_id, time_block)
stopifnot(nrow(prediction_ids) == nrow(probs))
predictions_with_ids <- dplyr::bind_cols(prediction_ids, probs)

#Save Probability 
write_csv(predictions_with_ids, 
paste0(project_location, site, '_project_output/sl_fit_probabilities.csv' ))

#' 
#' 
## ----Calculation of Non-Eligible Weights for Those that Develop ICU Indication in the ED----
raw_probs <- calc_weights |>
  left_join(predictions_with_ids, by=c('hospital_block_id', 'time_block')) |>
  select(hospital_block_id, time_block, icu_indication_probs=icu_indication, 
         still_in_ed_probs=still_in_ed, to_icu_probs=to_icu,
         to_imc_probs=to_imc, to_ward_probs=to_ward) 

#Outcome Only Multinomial Model for Numerator
numerator_fit <- nnet::multinom(adherence_cat ~ time_block , data=calc_weights)
new_data <- calc_weights |>
  dplyr::distinct(time_block) |>
  dplyr::arrange(time_block)

num_probs <- predict(
  numerator_fit,
  newdata = new_data,
  type = "probs"
)
num_probs <- cbind(new_data, as.data.frame(num_probs)) |>
  pivot_longer(
    cols=c("icu_indication", "still_in_ed", "to_icu", "to_imc", "to_ward" ),
    values_to='numerator_probs',
    names_to='adherence_cat')

##Join with Both Arms to Calculate Weights for Each Clone
both.arms_sl_weights <- both.arms |>
  select(hospital_block_id, time_block, icu_arm,  unit_location, block_status) |>
  left_join(raw_probs, by=c('hospital_block_id', 'time_block')) |>
  left_join(calc_weights %>% select(hospital_block_id, time_block, adherence_cat), 
            by=c('hospital_block_id', 'time_block')) |>
  mutate(
    block_weight_denom=fcase(
      icu_arm==0 & time_block<24, to_imc_probs + still_in_ed_probs + icu_indication_probs,
      icu_arm==0 & time_block==24, to_imc_probs + icu_indication_probs,
      icu_arm==1 & time_block<24, to_icu_probs + still_in_ed_probs + icu_indication_probs,
      icu_arm==1 & time_block==24, to_icu_probs + icu_indication_probs
    )
  ) |>
  left_join(num_probs)

write_csv(both.arms_sl_weights, paste0(project_location, '/private_tables/both.arms_sl_weights.csv'))

#' 
## ----Calculation of Non-Eligible Weights for Those that Remain in ED----------


#' 
## ----Calculate the Weights to be Used in the Outcome Model--------------------
#These are the cumulative product of inverse probability of censoring weights and weights for those that are removed because they are ineligible 


#' 
