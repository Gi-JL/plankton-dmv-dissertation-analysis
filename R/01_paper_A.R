############################
# Analysis code for Paper A of the cumulative dissertation
# Conditional and mixed logit models and related heterogeneity analyses
# Input-data sources and download instructions: see data/README.md
############################

# Run this script from the repository root:
# source("R/01_paper_A.R")
#
# Alternatively, from a terminal:
# Rscript R/01_paper_A.R
#
# Relative paths such as "./data/..." assume that the repository root
# is the current working directory.

############## Overview ###############
# 1. Conditional logit model in WTP-space
# 2. Conditional logit model in WTP-space: interaction between stable plankton composition and carbon sequestration
# 3. Mixed logit in WTP-space (including individual betas and linear regression on MPA-preference)
# 4. Conditional logit models in WTP-space: inland vs. coast
# 5. Conditional logit model in WTP-space: interaction with distance from the coast
# 6. Conditional logit models in WTP-space: country comparison
# 7. Pooled conditional logit model in WTP-space: country comparison
# 8. Conditional logit models in WTP-space: central south comparison
# 9. Conditional logit model in WTP-space: interaction with income
# 10. Documentation of starting value search (Conditional logit model in WTP-space)
# 11. Documentation of starting value search (Mixed logit model in WTP-space)
#######################################

############################
# Apollo version
############################

# All logit models reported in this paper were estimated using
# Apollo version 0.3.5. The script therefore installs this version
# automatically if another version is currently installed.

required_apollo_version <- "0.3.5"
required_rsghb_version <- "1.2.2"

installed_packages <- utils::installed.packages()

apollo_installed <-
  "apollo" %in% rownames(installed_packages)

installed_apollo_version <-
  if (apollo_installed) {
    installed_packages["apollo", "Version"]
  } else {
    NA_character_
  }

# Install the required Apollo version if necessary.
if (
  !apollo_installed ||
  installed_apollo_version != required_apollo_version
) {

  # An already loaded Apollo namespace cannot safely be replaced.
  if ("apollo" %in% loadedNamespaces()) {
    stop(
      "Apollo ", installed_apollo_version,
      " is already loaded in the current R session, but Apollo ",
      required_apollo_version,
      " is required. Restart R and rerun the script."
    )
  }

  message(
    "Apollo ",
    required_apollo_version,
    " is required. Installing the required version..."
  )

  # remotes is used because Apollo 0.3.5 and its RSGHB dependency
  # are archived CRAN package versions.
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # RSGHB 1.2.2 is an archived dependency of Apollo 0.3.5.
  installed_packages <- utils::installed.packages()

  if (
    !"RSGHB" %in% rownames(installed_packages) ||
    installed_packages["RSGHB", "Version"] != required_rsghb_version
  ) {
    remotes::install_version(
      "RSGHB",
      version = required_rsghb_version,
      repos = "https://cloud.r-project.org",
      dependencies = NA,
      upgrade = "never",
      force = TRUE
    )
  }

  # Install the Apollo version used for the reported analyses.
  remotes::install_version(
    "apollo",
    version = required_apollo_version,
    repos = "https://cloud.r-project.org",
    dependencies = NA,
    upgrade = "never",
    force = TRUE
  )
}

# Verify the installed version before continuing.
installed_packages <- utils::installed.packages()

if (
  !"apollo" %in% rownames(installed_packages) ||
  installed_packages["apollo", "Version"] != required_apollo_version
) {
  stop(
    "Installation of Apollo ",
    required_apollo_version,
    " was not successful. The analysis cannot continue."
  )
}

message(
  "Using Apollo ",
  installed_packages["apollo", "Version"]
)

############################
# Output directories
############################

dir.create(
  "./outputs/paper_A/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "./outputs/paper_A/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


######################
# 1. Conditional logit model in WTP-space
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(dce_answers, "RID")

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(pref1 = str_replace(pref1, "A", ""))

dce_data <- merge(choi, design, by = "choiceset", all = TRUE)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom = ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom = ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom = ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

# Sorting all columns behind column 5 (alphabetically by column names)
cols_to_sort <- names(dce_data)[5:ncol(dce_data)]
sorted_cols <- sort(cols_to_sort)

# Rearranging the dataframe
dce_data <- dce_data[, c(names(dce_data)[1:4], sorted_cols)]

database <- merge(dce_data, other_info, by = "RID", all = TRUE)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database"])

### Initialise Apollo and core settings
apollo_initialise()

apollo_control = list(
  modelName = "Clogit",
  modelDescr = "Conditional Logit in WTP-space",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = 0,
  b10 = 0,
  b11 = 0,
  b2 = 0,
  b30 = 0,
  b31 = 0,
  b4 = 0
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_output <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/01_baseline_conditional_logit.rds"
)

model_output$Parameter <- rownames(model_output)
rownames(model_output) <- NULL
model_output <- model_output[
  , c("Parameter", setdiff(names(model_output), "Parameter"))
]

write.csv(
  model_output,
  "./outputs/paper_A/tables/01_baseline_conditional_logit.csv",
  row.names = FALSE
)

rm(list=ls()[!ls() %in% c("model", "database")])


######################
# 2. Conditional logit model in WTP-space:
#    interaction between stable plankton composition and carbon sequestration
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(dce_answers, "RID")

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(pref1 = str_replace(pref1, "A", ""))

dce_data <- merge(choi, design, by = "choiceset", all = TRUE)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom = ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom = ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom = ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <- names(dce_data)[5:ncol(dce_data)]
sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[, c(names(dce_data)[1:4], sorted_cols)]

database <- merge(dce_data, other_info, by = "RID", all = TRUE)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database"])

apollo_initialise()

apollo_control = list(
  modelName = "Clogit",
  modelDescr = "Conditional Logit in WTP-space",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = 0,
  b10 = 0,
  b11 = 0,
  b2 = 0,
  b30 = 0,
  b31 = 0,
  b4 = 0,
  bi = 0
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly +
        bi * alt1.plankton_stable * alt1.carbon -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly +
        bi * alt2.plankton_stable * alt2.carbon -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly +
        bi * alt3.plankton_stable * alt3.carbon -
        alt3.cost
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_output <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/02_plankton_carbon_interaction.rds"
)

model_output$Parameter <- rownames(model_output)
rownames(model_output) <- NULL
model_output <- model_output[
  , c("Parameter", setdiff(names(model_output), "Parameter"))
]

write.csv(
  model_output,
  "./outputs/paper_A/tables/02_plankton_carbon_interaction.csv",
  row.names = FALSE
)

rm(list=ls()[!ls() %in% c("model", "database")])


######################
# 3. Mixed logit in WTP-space
#    (including individual betas and linear regression on MPA-preference)
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010
### 3) "PCODE_PT_2020_3035.shp" and associated shapefile components
### --> Eurostat GISCO
### 4) "Europe_coastline/Europe_coastline.shp" and associated shapefile components
### --> EEA coastline for analysis (line), version 2.0, September 2015

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(dce_answers, "RID")

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(pref1 = str_replace(pref1, "A", ""))

dce_data <- merge(choi, design, by = "choiceset", all = TRUE)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom = ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom = ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom = ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <- names(dce_data)[5:ncol(dce_data)]
sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[, c(names(dce_data)[1:4], sorted_cols)]

database <- merge(dce_data, other_info, by = "RID", all = TRUE)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database"])

apollo_initialise()

apollo_control = list(
  modelName = "Mixed Logit",
  modelDescr = "Mixed Logit in WTP-space",
  indivID = "RID",
  mixing = TRUE,
  nCores = 5
)

apollo_beta = c(
  mu_log_b0 = -88.0667,
  sigma_log_b0 = 54.8022,
  mu_log_b10 = 48.5357,
  sigma_log_b10 = 46.7755,
  mu_log_b11 = 15.2579,
  sigma_log_b11 = -45.4263,
  mu_log_b2 = 0.2546,
  sigma_log_b2 = 0.6713,
  mu_log_b30 = 35.6045,
  sigma_log_b30 = -29.5274,
  mu_log_b31 = 31.1781,
  sigma_log_b31 = -16.6907,
  mu_log_b4 = -3.6243,
  sigma_log_b4 = -0.3534
)

apollo_fixed = c()

# Fix the random seed for reproducibility of the MLHS draws in the repository analysis.
# No random seed was fixed in the original analysis reported in Paper A.
set.seed(12345)

apollo_draws = list(
  interDrawsType = "mlhs",
  interNDraws = 5000,
  interNormDraws = c(
    "draws_0",
    "draws_1",
    "draws_2",
    "draws_3",
    "draws_4",
    "draws_5",
    "draws_6"
  )
)

apollo_randCoeff = function(apollo_beta, apollo_inputs) {

  randcoeff = list()

  randcoeff[["b0"]] =
    mu_log_b0 + sigma_log_b0 * draws_0

  randcoeff[["b10"]] =
    mu_log_b10 + sigma_log_b10 * draws_1

  randcoeff[["b11"]] =
    mu_log_b11 + sigma_log_b11 * draws_2

  randcoeff[["b2"]] =
    mu_log_b2 + sigma_log_b2 * draws_3

  randcoeff[["b30"]] =
    mu_log_b30 + sigma_log_b30 * draws_4

  randcoeff[["b31"]] =
    mu_log_b31 + sigma_log_b31 * draws_5

  randcoeff[["b4"]] =
    -exp(mu_log_b4 + sigma_log_b4 * draws_6)

  return(randcoeff)
}

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_output <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/03_mixed_logit.rds"
)

model_output$Parameter <- rownames(model_output)
rownames(model_output) <- NULL
model_output <- model_output[
  , c("Parameter", setdiff(names(model_output), "Parameter"))
]

write.csv(
  model_output,
  "./outputs/paper_A/tables/03_mixed_logit.csv",
  row.names = FALSE
)

# Compute conditional individual betas
individual_betas <- apollo_conditionals(
  model,
  apollo_probabilities,
  apollo_inputs
)

individual_betas$b0$ID

length(unique(individual_betas$b0$ID))

setdiff(
  1:length(unique(individual_betas$b0$ID)),
  unique(individual_betas$b0$ID)
)

individual_betas_table <- cbind(
  individual_betas$b0,
  individual_betas$b10[,2:3],
  individual_betas$b11[,2:3],
  individual_betas$b2[,2:3],
  individual_betas$b30[,2:3],
  individual_betas$b31[,2:3],
  individual_betas$b4[,2:3]
)

colnames(individual_betas_table) <- c(
  "RID",
  "b0_est", "b0_sd",
  "b10_est", "b10_sd",
  "b11_est", "b11_sd",
  "b2_est", "b2_sd",
  "b30_est", "b30_sd",
  "b31_est", "b31_sd",
  "b4_est", "b4_sd"
)

rm(list=ls()[!ls() %in% c("individual_betas_table", "model")])

library(sf)

dce_data <- read.csv2("./data/DMV data.csv")

dce_data$Postcode[dce_data$WS_ID == 7] <-
  paste0("0", dce_data$Postcode[dce_data$WS_ID == 7])

eu_pc <- st_read("./data/PCODE_PT_2020_3035.shp")

dce_data <- dce_data %>%
  rename(
    RID = ID
  )

dce_data <- dce_data %>%
  drop_na(Q1)

dce_data$Postcode[dce_data$Country == "Poland"] <-
  paste0(
    "PL_",
    dce_data$Postcode[dce_data$Country == "Poland"]
  )

dce_data$Postcode[dce_data$Country == "Poland"] <-
  sub(
    "^(.{5})",
    "\\1-",
    dce_data$Postcode[dce_data$Country == "Poland"]
  )

dce_data$Postcode[dce_data$Country == "Italy"] <-
  paste0(
    "IT_",
    dce_data$Postcode[dce_data$Country == "Italy"]
  )

dce_data$Postcode[dce_data$Country == "Basque"] <-
  paste0(
    "ES_",
    dce_data$Postcode[dce_data$Country == "Basque"]
  )

dce_data$Postcode[dce_data$Country == "Germany"] <-
  paste0(
    "DE_",
    dce_data$Postcode[dce_data$Country == "Germany"]
  )

dce_data$Postcode[dce_data$Country == "France"] <-
  paste0(
    "FR_",
    dce_data$Postcode[dce_data$Country == "France"]
  )

dce_data <- dce_data %>%
  rename(
    PC_CNTR = Postcode
  )

coords <- eu_pc %>%
  dplyr::select(PC_CNTR) %>%
  as.data.frame()

dce_data <- merge(
  dce_data,
  coords,
  by = "PC_CNTR"
)

dce_data <- dce_data %>%
  relocate(PC_CNTR, .after = ISCED)

dce_data <- dce_data %>%
  relocate(geometry, .after = PC_CNTR)

dce_data_points <- sf::st_as_sf(dce_data)

coastline <- st_read(
  "./data/Europe_coastline/Europe_coastline.shp"
)

st_crs(coastline)
st_crs(dce_data_points)

dce_data_points <- st_transform(
  dce_data_points,
  crs = st_crs(coastline)
)

dce_data_points$dist2coast <-
  apply(
    st_distance(dce_data_points, coastline),
    1,
    min
  )

point_df <- as.data.frame(dce_data_points)

dist2coast <- point_df[, c("RID", "dist2coast")]

dce_data <- merge(
  dce_data,
  dist2coast,
  by = "RID"
)

individual_betas_table$FoH <-
  individual_betas_table$b30_est -
  individual_betas_table$b31_est

dce_data$Age <-
  2023 - dce_data$`Year.of.birth`

individual_betas_table <- individual_betas_table %>%
  mutate(
    Age = dce_data$Age[match(RID, dce_data$RID)],
    dist2coast = dce_data$dist2coast[match(RID, dce_data$RID)],
    Gender = dce_data$Gender[match(RID, dce_data$RID)],
    ISCED = dce_data$ISCED[match(RID, dce_data$RID)],
    `Net income group` =
      dce_data$`Net.income.group`[match(RID, dce_data$RID)]
  )

individual_betas_table <- individual_betas_table %>%
  mutate(
    age_group = ifelse(
      Age <= 34,
      1,
      ifelse(Age <= 54, 2, 3)
    )
  )

individual_betas_table <- individual_betas_table %>%
  mutate(
    `Net income group_n` = ifelse(
      `Net income group` <= 4,
      1,
      ifelse(`Net income group` <= 8, 2, 3)
    )
  )

individual_betas_table$age_group <-
  as.factor(individual_betas_table$age_group)

individual_betas_table$`Net income group_n` <-
  as.factor(individual_betas_table$`Net income group_n`)

individual_betas_table <- individual_betas_table %>%
  mutate(
    ISCED_group = ifelse(ISCED >= 6, 2, 1)
  )

individual_betas_table$ISCED_group <-
  as.factor(individual_betas_table$ISCED_group)

individual_betas_table$dist2coast_km <-
  individual_betas_table$dist2coast / 1000

write.csv(
  individual_betas_table,
  "./outputs/paper_A/tables/03_individual_specific_WTP.csv",
  row.names = FALSE
)

lm_model <- lm(
  FoH ~
    age_group +
    dist2coast_km +
    Gender +
    ISCED_group +
    `Net income group_n`,
  data = individual_betas_table
)

summary(lm_model)

saveRDS(
  lm_model,
  "./outputs/paper_A/models/03_mpa_linear_regression.rds"
)

lm_coefficients <- as.data.frame(
  coef(summary(lm_model))
)

lm_coefficients$Parameter <- rownames(lm_coefficients)
rownames(lm_coefficients) <- NULL

lm_coefficients <- lm_coefficients[
  , c("Parameter", setdiff(names(lm_coefficients), "Parameter"))
]

write.csv(
  lm_coefficients,
  "./outputs/paper_A/tables/03_mpa_linear_regression_coefficients.csv",
  row.names = FALSE
)

lm_fit <- data.frame(
  Observations = nobs(lm_model),
  R2 = summary(lm_model)$r.squared,
  Adjusted_R2 = summary(lm_model)$adj.r.squared,
  AIC = AIC(lm_model),
  BIC = BIC(lm_model)
)

write.csv(
  lm_fit,
  "./outputs/paper_A/tables/03_mpa_linear_regression_fit.csv",
  row.names = FALSE
)

# Diagnostic check of non-price mean estimates (Sarrias 2020)
diagnostic_check <- data.frame(
  Parameter = c(
    "ASC",
    "Stable plankton composition",
    "Control of bloom events",
    "Carbon sequestration",
    "Fully protected MPA",
    "Highly protected MPA"
  ),
  Relative_deviation = c(
    (mean(individual_betas_table$b0_est) - model$estimate[1]) /
      model$estimate[1],
    (mean(individual_betas_table$b10_est) - model$estimate[3]) /
      model$estimate[3],
    (mean(individual_betas_table$b11_est) - model$estimate[5]) /
      model$estimate[5],
    (mean(individual_betas_table$b2_est) - model$estimate[7]) /
      model$estimate[7],
    (mean(individual_betas_table$b30_est) - model$estimate[9]) /
      model$estimate[9],
    (mean(individual_betas_table$b31_est) - model$estimate[11]) /
      model$estimate[11]
  )
)

print(diagnostic_check)

write.csv(
  diagnostic_check,
  "./outputs/paper_A/tables/03_individual_WTP_diagnostic_check.csv",
  row.names = FALSE
)


######################
# 4. Conditional logit models in WTP-space: inland vs. coast
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(dce_answers, "RID")

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(pref1 = str_replace(pref1, "A", ""))

dce_data <- merge(
  choi,
  design,
  by = "choiceset",
  all = TRUE
)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom = ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom = ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom = ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <- names(dce_data)[5:ncol(dce_data)]
sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[, c(names(dce_data)[1:4], sorted_cols)]

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database_all <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database_all"])

### Inland model
apollo_initialise()

apollo_control = list(
  modelName = "Clogit inland vs coast",
  modelDescr = "Conditional Logit in WTP-space",
  indivID = "RID",
  mixing = FALSE
)

database <- database_all %>%
  filter(Coastal == 0)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_inland <- model

est_inland <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_inland[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_inland[["Estimate"]] / est_inland[["Rob.s.e."]]))

saveRDS(
  model_inland,
  "./outputs/paper_A/models/04_inland_conditional_logit.rds"
)

est_inland_output <- est_inland
est_inland_output$Parameter <- rownames(est_inland_output)
rownames(est_inland_output) <- NULL
est_inland_output <- est_inland_output[
  , c("Parameter", setdiff(names(est_inland_output), "Parameter"))
]

write.csv(
  est_inland_output,
  "./outputs/paper_A/tables/04_inland_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_inland[["Rob.p(1-sided)"]]
std_errors <- est_inland[["Rob.s.e."]]

results_df_inland <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_inland$CI_Lower <-
  results_df_inland$Estimate -
  1.96 * results_df_inland$Std_Error

results_df_inland$CI_Upper <-
  results_df_inland$Estimate +
  1.96 * results_df_inland$Std_Error

results_df_inland$Parameter <-
  rownames(results_df_inland)

results_df_inland <- results_df_inland[
  , c(
    "Parameter",
    setdiff(names(results_df_inland), "Parameter")
  )
]

rownames(results_df_inland) <- NULL

print(results_df_inland)

rm(
  list = setdiff(
    ls(),
    c(
      "model_inland",
      "results_df_inland",
      "est_inland",
      "database_all"
    )
  )
)

### Coast model
apollo_initialise()

apollo_control = list(
  modelName = "Clogit inland vs coast",
  modelDescr = "Conditional Logit in WTP-space",
  indivID = "RID",
  mixing = FALSE
)

database <- database_all %>%
  filter(Coastal == 1)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_coast <- model

est_coast <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_coast[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_coast[["Estimate"]] / est_coast[["Rob.s.e."]]))

saveRDS(
  model_coast,
  "./outputs/paper_A/models/04_coast_conditional_logit.rds"
)

est_coast_output <- est_coast
est_coast_output$Parameter <- rownames(est_coast_output)
rownames(est_coast_output) <- NULL
est_coast_output <- est_coast_output[
  , c("Parameter", setdiff(names(est_coast_output), "Parameter"))
]

write.csv(
  est_coast_output,
  "./outputs/paper_A/tables/04_coast_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_coast[["Rob.p(1-sided)"]]
std_errors <- est_coast[["Rob.s.e."]]

results_df_coast <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_coast$CI_Lower <-
  results_df_coast$Estimate -
  1.96 * results_df_coast$Std_Error

results_df_coast$CI_Upper <-
  results_df_coast$Estimate +
  1.96 * results_df_coast$Std_Error

results_df_coast$Parameter <-
  rownames(results_df_coast)

results_df_coast <- results_df_coast[
  , c(
    "Parameter",
    setdiff(names(results_df_coast), "Parameter")
  )
]

rownames(results_df_coast) <- NULL

print(results_df_coast)

combined_results <- rbind(
  transform(
    results_df_coast,
    Group = "coast"
  ),
  transform(
    results_df_inland,
    Group = "inland"
  )
)

str(combined_results)

write.csv(
  combined_results,
  "./outputs/paper_A/tables/04_coast_inland_comparison.csv",
  row.names = FALSE
)

rm(
  list = ls()[
    !ls() %in%
      c(
        "combined_results",
        "model_inland",
        "model_coast",
        "database"
      )
  ]
)


######################
# 5. Conditional logit model in WTP-space:
#    interaction with distance from the coast
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010
### 3) "PCODE_PT_2020_3035.shp" and associated shapefile components
### --> Eurostat GISCO
### 4) "Europe_coastline/Europe_coastline.shp" and associated shapefile components
### --> EEA coastline for analysis (line), version 2.0, September 2015

rm(list=ls())
library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)
library(sf)

dce_data <- read.csv2("./data/DMV data.csv")

dce_data$Postcode[dce_data$WS_ID == 7] <-
  paste0(
    "0",
    dce_data$Postcode[dce_data$WS_ID == 7]
  )

eu_pc <- st_read(
  "./data/PCODE_PT_2020_3035.shp"
)

dce_data <- dce_data %>%
  rename(
    RID = ID
  )

dce_data <- dce_data %>%
  drop_na(Q1)

dce_data$Postcode[dce_data$Country == "Poland"] <-
  paste0(
    "PL_",
    dce_data$Postcode[dce_data$Country == "Poland"]
  )

dce_data$Postcode[dce_data$Country == "Poland"] <-
  sub(
    "^(.{5})",
    "\\1-",
    dce_data$Postcode[dce_data$Country == "Poland"]
  )

dce_data$Postcode[dce_data$Country == "Italy"] <-
  paste0(
    "IT_",
    dce_data$Postcode[dce_data$Country == "Italy"]
  )

dce_data$Postcode[dce_data$Country == "Basque"] <-
  paste0(
    "ES_",
    dce_data$Postcode[dce_data$Country == "Basque"]
  )

dce_data$Postcode[dce_data$Country == "Germany"] <-
  paste0(
    "DE_",
    dce_data$Postcode[dce_data$Country == "Germany"]
  )

dce_data$Postcode[dce_data$Country == "France"] <-
  paste0(
    "FR_",
    dce_data$Postcode[dce_data$Country == "France"]
  )

dce_data <- dce_data %>%
  rename(
    PC_CNTR = Postcode
  )

coords <- eu_pc %>%
  dplyr::select(PC_CNTR) %>%
  as.data.frame()

dce_data <- merge(
  dce_data,
  coords,
  by = "PC_CNTR"
)

dce_data <- dce_data %>%
  relocate(PC_CNTR, .after = ISCED)

dce_data <- dce_data %>%
  relocate(geometry, .after = PC_CNTR)

dce_data_points <- sf::st_as_sf(dce_data)

coastline <- st_read(
  "./data/Europe_coastline/Europe_coastline.shp"
)

st_crs(coastline)
st_crs(dce_data_points)

dce_data_points <- st_transform(
  dce_data_points,
  crs = st_crs(coastline)
)

dce_data_points$dist2coast <-
  apply(
    st_distance(
      dce_data_points,
      coastline
    ),
    1,
    min
  )

point_df <- as.data.frame(dce_data_points)

dist2coast <- point_df[
  , c("RID", "dist2coast")
]

dce_data <- merge(
  dce_data,
  dist2coast,
  by = "RID"
)

dce_data$dist2coast_km <-
  dce_data$dist2coast / 1000

dce_data$dist2coast <- NULL

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,c(1, 45:60)]

other_info <- dce_data[
  ,-c(12:13, 45:60)
]

colnames(dce_answers)[2:17] <- 1:16

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(
    pref1 = str_replace(
      pref1,
      "A",
      ""
    )
  )

dce_data <- merge(
  choi,
  design,
  by = "choiceset",
  all = TRUE
)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom =
      ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom =
      ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[5:ncol(dce_data)]

sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[
  , c(names(dce_data)[1:4], sorted_cols)
]

database <- merge(
  dce_data,
  other_info,
  by = "RID"
)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84,
    RID != 129
  )

rm(list=ls()[ls() != "database"])

apollo_initialise()

apollo_control = list(
  modelName = "Clogit interaction",
  modelDescr =
    "Conditional Logit in WTP-space - interaction dist2coast",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844,
  bi10 = 0,
  bi11 = 0,
  bi2 = 0,
  bi30 = 0,
  bi31 = 0,
  bi4 = 0,
  bi0 = 0
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -(b4 + bi4 * dist2coast_km) * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost +
        bi10 * alt1.plankton_stable * dist2coast_km +
        bi11 * alt1.plankton_nobloom * dist2coast_km +
        bi2 * dist2coast_km * alt1.carbon +
        bi30 * alt1.mpa_fully * dist2coast_km +
        bi31 * alt1.mpa_highly * dist2coast_km
    )

  V[['alt2']] =
    -(b4 + bi4 * dist2coast_km) * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost +
        bi10 * alt2.plankton_stable * dist2coast_km +
        bi11 * alt2.plankton_nobloom * dist2coast_km +
        bi2 * dist2coast_km * alt2.carbon +
        bi30 * alt2.mpa_fully * dist2coast_km +
        bi31 * alt2.mpa_highly * dist2coast_km
    )

  V[['alt3']] =
    -(b4 + bi4 * dist2coast_km) * (
      b0 +
        bi0 * dist2coast_km +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost +
        bi10 * alt3.plankton_stable * dist2coast_km +
        bi11 * alt3.plankton_nobloom * dist2coast_km +
        bi2 * dist2coast_km * alt3.carbon +
        bi30 * alt3.mpa_fully * dist2coast_km +
        bi31 * alt3.mpa_highly * dist2coast_km
    )

  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] =
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P = apollo_panelProd(
    P,
    apollo_inputs,
    functionality
  )

  P = apollo_prepareProb(
    P,
    apollo_inputs,
    functionality
  )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_output <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/05_distance_interaction.rds"
)

model_output$Parameter <- rownames(model_output)
rownames(model_output) <- NULL
model_output <- model_output[
  , c("Parameter", setdiff(names(model_output), "Parameter"))
]

write.csv(
  model_output,
  "./outputs/paper_A/tables/05_distance_interaction.csv",
  row.names = FALSE
)

rm(list=ls()[!ls() %in% c("model", "database")])


######################
# 6. Conditional logit models in WTP-space: country comparison
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(
  dce_answers,
  "RID"
)

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <- as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(
    pref1 = str_replace(pref1, "A", "")
  )

dce_data <- merge(
  choi,
  design,
  by = "choiceset",
  all = TRUE
)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom =
      ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom =
      ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[5:ncol(dce_data)]

sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[
  , c(names(dce_data)[1:4], sorted_cols)
]

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database"])

apollo_initialise()

apollo_control = list(
  modelName = "Clogit",
  modelDescr = "Conditional Logit in WTP space",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

### Poland
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Country == "Poland"
  )

  P[['model']] =
    apollo_mnl(mnl_settings, functionality)

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_poland <- model

est_poland <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_poland[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_poland[["Estimate"]] / est_poland[["Rob.s.e."]]))

saveRDS(
  model_poland,
  "./outputs/paper_A/models/06_poland_conditional_logit.rds"
)

est_poland_output <- est_poland
est_poland_output$Parameter <- rownames(est_poland_output)
rownames(est_poland_output) <- NULL
est_poland_output <- est_poland_output[
  , c("Parameter", setdiff(names(est_poland_output), "Parameter"))
]

write.csv(
  est_poland_output,
  "./outputs/paper_A/tables/06_poland_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_poland[["Rob.p(1-sided)"]]
std_errors <- est_poland[["Rob.s.e."]]

results_df_poland <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_poland$CI_Lower <-
  results_df_poland$Estimate -
  1.96 * results_df_poland$Std_Error

results_df_poland$CI_Upper <-
  results_df_poland$Estimate +
  1.96 * results_df_poland$Std_Error

results_df_poland$Parameter <-
  rownames(results_df_poland)

results_df_poland <- results_df_poland[
  , c(
    "Parameter",
    setdiff(names(results_df_poland), "Parameter")
  )
]

rownames(results_df_poland) <- NULL

print(results_df_poland)


### Italy
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Country == "Italy"
  )

  P[['model']] =
    apollo_mnl(mnl_settings, functionality)

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_italy <- model

est_italy <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_italy[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_italy[["Estimate"]] / est_italy[["Rob.s.e."]]))

saveRDS(
  model_italy,
  "./outputs/paper_A/models/06_italy_conditional_logit.rds"
)

est_italy_output <- est_italy
est_italy_output$Parameter <- rownames(est_italy_output)
rownames(est_italy_output) <- NULL
est_italy_output <- est_italy_output[
  , c("Parameter", setdiff(names(est_italy_output), "Parameter"))
]

write.csv(
  est_italy_output,
  "./outputs/paper_A/tables/06_italy_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_italy[["Rob.p(1-sided)"]]
std_errors <- est_italy[["Rob.s.e."]]

results_df_italy <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_italy$CI_Lower <-
  results_df_italy$Estimate -
  1.96 * results_df_italy$Std_Error

results_df_italy$CI_Upper <-
  results_df_italy$Estimate +
  1.96 * results_df_italy$Std_Error

results_df_italy$Parameter <-
  rownames(results_df_italy)

results_df_italy <- results_df_italy[
  , c(
    "Parameter",
    setdiff(names(results_df_italy), "Parameter")
  )
]

rownames(results_df_italy) <- NULL

print(results_df_italy)


### Basque Country
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Country == "Basque"
  )

  P[['model']] =
    apollo_mnl(mnl_settings, functionality)

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_basque <- model

est_basque <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_basque[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_basque[["Estimate"]] / est_basque[["Rob.s.e."]]))

saveRDS(
  model_basque,
  "./outputs/paper_A/models/06_basque_conditional_logit.rds"
)

est_basque_output <- est_basque
est_basque_output$Parameter <- rownames(est_basque_output)
rownames(est_basque_output) <- NULL
est_basque_output <- est_basque_output[
  , c("Parameter", setdiff(names(est_basque_output), "Parameter"))
]

write.csv(
  est_basque_output,
  "./outputs/paper_A/tables/06_basque_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_basque[["Rob.p(1-sided)"]]
std_errors <- est_basque[["Rob.s.e."]]

results_df_basque <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_basque$CI_Lower <-
  results_df_basque$Estimate -
  1.96 * results_df_basque$Std_Error

results_df_basque$CI_Upper <-
  results_df_basque$Estimate +
  1.96 * results_df_basque$Std_Error

results_df_basque$Parameter <-
  rownames(results_df_basque)

results_df_basque <- results_df_basque[
  , c(
    "Parameter",
    setdiff(names(results_df_basque), "Parameter")
  )
]

rownames(results_df_basque) <- NULL

print(results_df_basque)


### Germany
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Country == "Germany"
  )

  P[['model']] =
    apollo_mnl(mnl_settings, functionality)

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_germany <- model

est_germany <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_germany[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_germany[["Estimate"]] / est_germany[["Rob.s.e."]]))

saveRDS(
  model_germany,
  "./outputs/paper_A/models/06_germany_conditional_logit.rds"
)

est_germany_output <- est_germany
est_germany_output$Parameter <- rownames(est_germany_output)
rownames(est_germany_output) <- NULL
est_germany_output <- est_germany_output[
  , c("Parameter", setdiff(names(est_germany_output), "Parameter"))
]

write.csv(
  est_germany_output,
  "./outputs/paper_A/tables/06_germany_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_germany[["Rob.p(1-sided)"]]
std_errors <- est_germany[["Rob.s.e."]]

results_df_germany <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_germany$CI_Lower <-
  results_df_germany$Estimate -
  1.96 * results_df_germany$Std_Error

results_df_germany$CI_Upper <-
  results_df_germany$Estimate +
  1.96 * results_df_germany$Std_Error

results_df_germany$Parameter <-
  rownames(results_df_germany)

results_df_germany <- results_df_germany[
  , c(
    "Parameter",
    setdiff(names(results_df_germany), "Parameter")
  )
]

rownames(results_df_germany) <- NULL

print(results_df_germany)


### France
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(alt1 = 1, alt2 = 2, alt3 = 3),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Country == "France"
  )

  P[['model']] =
    apollo_mnl(mnl_settings, functionality)

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_france <- model

est_france <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings = list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
est_france[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_france[["Estimate"]] / est_france[["Rob.s.e."]]))

saveRDS(
  model_france,
  "./outputs/paper_A/models/06_france_conditional_logit.rds"
)

est_france_output <- est_france
est_france_output$Parameter <- rownames(est_france_output)
rownames(est_france_output) <- NULL
est_france_output <- est_france_output[
  , c("Parameter", setdiff(names(est_france_output), "Parameter"))
]

write.csv(
  est_france_output,
  "./outputs/paper_A/tables/06_france_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <- model$estimate
p_values <- est_france[["Rob.p(1-sided)"]]
std_errors <- est_france[["Rob.s.e."]]

results_df_france <- data.frame(
  Estimate = coef_estimates,
  Std_Error = std_errors,
  p_value = p_values
)

results_df_france$CI_Lower <-
  results_df_france$Estimate -
  1.96 * results_df_france$Std_Error

results_df_france$CI_Upper <-
  results_df_france$Estimate +
  1.96 * results_df_france$Std_Error

results_df_france$Parameter <-
  rownames(results_df_france)

results_df_france <- results_df_france[
  , c(
    "Parameter",
    setdiff(names(results_df_france), "Parameter")
  )
]

rownames(results_df_france) <- NULL

print(results_df_france)

combined_results <- rbind(
  transform(results_df_italy, Group = "Italy"),
  transform(results_df_poland, Group = "Poland"),
  transform(results_df_basque, Group = "Basque"),
  transform(results_df_germany, Group = "Germany"),
  transform(results_df_france, Group = "France")
)

write.csv(
  combined_results,
  "./outputs/paper_A/tables/06_country_comparison.csv",
  row.names = FALSE
)

rm(
  list = ls()[
    !ls() %in%
      c(
        "combined_results",
        "model_poland",
        "model_italy",
        "model_basque",
        "model_germany",
        "model_france",
        "database"
      )
  ]
)


######################
# 7. Pooled conditional logit model in WTP-space: country comparison
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(
  dce_answers,
  "RID"
)

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <-
  as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(
    pref1 = str_replace(pref1, "A", "")
  )

dce_data <- merge(
  choi,
  design,
  by = "choiceset",
  all = TRUE
)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable = alt1.plankton == 1,
    alt2.plankton_stable = alt2.plankton == 1,
    alt3.plankton_stable = alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(alt1.plankton == 3, FALSE, TRUE),
    alt2.plankton_nobloom =
      ifelse(alt2.plankton == 3, FALSE, TRUE),
    alt3.plankton_nobloom =
      ifelse(alt3.plankton == 3, FALSE, TRUE),
    alt1.mpa_fully = alt1.mpa == 1,
    alt2.mpa_fully = alt2.mpa == 1,
    alt3.mpa_fully = alt3.mpa == 1,
    alt1.mpa_highly = alt1.mpa == 2,
    alt2.mpa_highly = alt2.mpa == 2,
    alt3.mpa_highly = alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[5:ncol(dce_data)]

sorted_cols <- sort(cols_to_sort)

dce_data <- dce_data[
  , c(names(dce_data)[1:4], sorted_cols)
]

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

unique(database$Country)

database <- database %>%
  mutate(
    Country_Poland = Country == "Poland",
    Country_Italy = Country == "Italy",
    Country_Basque = Country == "Basque",
    Country_Germany = Country == "Germany",
    Country_France = Country == "France"
  )

rm(list=ls()[ls() != "database"])

apollo_initialise()

apollo_control = list(
  modelName = "Clogit",
  modelDescr =
    "Conditional Logit in WTP-space - country comparison with pooled model",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844,

  bi_0_fr = 0,
  bi_10_fr = 0,
  bi_11_fr = 0,
  bi_2_fr = 0,
  bi_30_fr = 0,
  bi_31_fr = 0,
  bi_4_fr = 0,

  bi_0_pl = 0,
  bi_10_pl = 0,
  bi_11_pl = 0,
  bi_2_pl = 0,
  bi_30_pl = 0,
  bi_31_pl = 0,
  bi_4_pl = 0,

  bi_0_it = 0,
  bi_10_it = 0,
  bi_11_it = 0,
  bi_2_it = 0,
  bi_30_it = 0,
  bi_31_it = 0,
  bi_4_it = 0,

  bi_0_bas = 0,
  bi_10_bas = 0,
  bi_11_bas = 0,
  bi_2_bas = 0,
  bi_30_bas = 0,
  bi_31_bas = 0,
  bi_4_bas = 0
)

apollo_fixed = c()

apollo_inputs = apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(
    apollo_beta,
    apollo_inputs
  )

  on.exit(
    apollo_detach(
      apollo_beta,
      apollo_inputs
    )
  )

  P = list()
  V = list()

  V[['alt1']] =
    -(
      b4 +
        bi_4_fr * Country_France +
        bi_4_pl * Country_Poland +
        bi_4_it * Country_Italy +
        bi_4_bas * Country_Basque
    ) * (
      b10 * alt1.plankton_stable +
        bi_10_fr * alt1.plankton_stable * Country_France +
        bi_10_pl * alt1.plankton_stable * Country_Poland +
        bi_10_it * alt1.plankton_stable * Country_Italy +
        bi_10_bas * alt1.plankton_stable * Country_Basque +

        b11 * alt1.plankton_nobloom +
        bi_11_fr * alt1.plankton_nobloom * Country_France +
        bi_11_pl * alt1.plankton_nobloom * Country_Poland +
        bi_11_it * alt1.plankton_nobloom * Country_Italy +
        bi_11_bas * alt1.plankton_nobloom * Country_Basque +

        b2 * alt1.carbon +
        bi_2_fr * alt1.carbon * Country_France +
        bi_2_pl * alt1.carbon * Country_Poland +
        bi_2_it * alt1.carbon * Country_Italy +
        bi_2_bas * alt1.carbon * Country_Basque +

        b30 * alt1.mpa_fully +
        bi_30_fr * alt1.mpa_fully * Country_France +
        bi_30_pl * alt1.mpa_fully * Country_Poland +
        bi_30_it * alt1.mpa_fully * Country_Italy +
        bi_30_bas * alt1.mpa_fully * Country_Basque +

        b31 * alt1.mpa_highly +
        bi_31_fr * alt1.mpa_highly * Country_France +
        bi_31_pl * alt1.mpa_highly * Country_Poland +
        bi_31_it * alt1.mpa_highly * Country_Italy +
        bi_31_bas * alt1.mpa_highly * Country_Basque -

        alt1.cost
    )

  V[['alt2']] =
    -(
      b4 +
        bi_4_fr * Country_France +
        bi_4_pl * Country_Poland +
        bi_4_it * Country_Italy +
        bi_4_bas * Country_Basque
    ) * (
      b10 * alt2.plankton_stable +
        bi_10_fr * alt2.plankton_stable * Country_France +
        bi_10_pl * alt2.plankton_stable * Country_Poland +
        bi_10_it * alt2.plankton_stable * Country_Italy +
        bi_10_bas * alt2.plankton_stable * Country_Basque +

        b11 * alt2.plankton_nobloom +
        bi_11_fr * alt2.plankton_nobloom * Country_France +
        bi_11_pl * alt2.plankton_nobloom * Country_Poland +
        bi_11_it * alt2.plankton_nobloom * Country_Italy +
        bi_11_bas * alt2.plankton_nobloom * Country_Basque +

        b2 * alt2.carbon +
        bi_2_fr * alt2.carbon * Country_France +
        bi_2_pl * alt2.carbon * Country_Poland +
        bi_2_it * alt2.carbon * Country_Italy +
        bi_2_bas * alt2.carbon * Country_Basque +

        b30 * alt2.mpa_fully +
        bi_30_fr * alt2.mpa_fully * Country_France +
        bi_30_pl * alt2.mpa_fully * Country_Poland +
        bi_30_it * alt2.mpa_fully * Country_Italy +
        bi_30_bas * alt2.mpa_fully * Country_Basque +

        b31 * alt2.mpa_highly +
        bi_31_fr * alt2.mpa_highly * Country_France +
        bi_31_pl * alt2.mpa_highly * Country_Poland +
        bi_31_it * alt2.mpa_highly * Country_Italy +
        bi_31_bas * alt2.mpa_highly * Country_Basque -

        alt2.cost
    )

  V[['alt3']] =
    -(
      b4 +
        bi_4_fr * Country_France +
        bi_4_pl * Country_Poland +
        bi_4_it * Country_Italy +
        bi_4_bas * Country_Basque
    ) * (
      b0 +
        bi_0_fr * Country_France +
        bi_0_pl * Country_Poland +
        bi_0_it * Country_Italy +
        bi_0_bas * Country_Basque +

        b10 * alt3.plankton_stable +
        bi_10_fr * alt3.plankton_stable * Country_France +
        bi_10_pl * alt3.plankton_stable * Country_Poland +
        bi_10_it * alt3.plankton_stable * Country_Italy +
        bi_10_bas * alt3.plankton_stable * Country_Basque +

        b11 * alt3.plankton_nobloom +
        bi_11_fr * alt3.plankton_nobloom * Country_France +
        bi_11_pl * alt3.plankton_nobloom * Country_Poland +
        bi_11_it * alt3.plankton_nobloom * Country_Italy +
        bi_11_bas * alt3.plankton_nobloom * Country_Basque +

        b2 * alt3.carbon +
        bi_2_fr * alt3.carbon * Country_France +
        bi_2_pl * alt3.carbon * Country_Poland +
        bi_2_it * alt3.carbon * Country_Italy +
        bi_2_bas * alt3.carbon * Country_Basque +

        b30 * alt3.mpa_fully +
        bi_30_fr * alt3.mpa_fully * Country_France +
        bi_30_pl * alt3.mpa_fully * Country_Poland +
        bi_30_it * alt3.mpa_fully * Country_Italy +
        bi_30_bas * alt3.mpa_fully * Country_Basque +

        b31 * alt3.mpa_highly +
        bi_31_fr * alt3.mpa_highly * Country_France +
        bi_31_pl * alt3.mpa_highly * Country_Poland +
        bi_31_it * alt3.mpa_highly * Country_Italy +
        bi_31_bas * alt3.mpa_highly * Country_Basque -

        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(
        alt1 = 1,
        alt2 = 2,
        alt3 = 3
      ),
    avail = 1,
    choiceVar = pref1,
    V = V
  )

  P[['model']] =
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_output <- as.data.frame(
  apollo_modelOutput(
    model,
    modelOutput_settings =
      list(printPVal = 1)
  )
)

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/07_pooled_country_interaction.rds"
)

model_output$Parameter <-
  rownames(model_output)

rownames(model_output) <- NULL

model_output <- model_output[
  , c(
    "Parameter",
    setdiff(
      names(model_output),
      "Parameter"
    )
  )
]

write.csv(
  model_output,
  "./outputs/paper_A/tables/07_pooled_country_interaction.csv",
  row.names = FALSE
)

rm(
  list = ls()[
    !ls() %in%
      c("model", "database")
  ]
)


######################
# 8. Conditional logit models in WTP-space: central south comparison
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <- read.csv2("./data/DMV data.csv")

design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(-Design, -starts_with("...")),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(!is.na(Choice.situation))

nsets <- nrow(design)

design <- design %>%
  rename(
    choiceset = Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <- design %>%
  relocate(Block, .before = choiceset)

dce_answers <- dce_data[,44:59]
other_info <- dce_data[,-(44:59)]

other_info <- other_info %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(2:75),
    condition = ~.x == 99
  )

colnames(dce_answers) <- 1:16

dce_answers <- tibble::rowid_to_column(
  dce_answers,
  "RID"
)

choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <-
  as.numeric(choi$choiceset)

choi <- choi %>%
  mutate(
    pref1 = str_replace(
      pref1,
      "A",
      ""
    )
  )

dce_data <- merge(
  choi,
  design,
  by = "choiceset",
  all = TRUE
)

dce_data <- dce_data %>%
  relocate(RID)

dce_data <- dce_data[
  order(dce_data[,1], dce_data[,2]),
]

row.names(dce_data) <- NULL

dce_data <- dce_data %>%
  mutate(
    alt1.plankton_stable =
      alt1.plankton == 1,
    alt2.plankton_stable =
      alt2.plankton == 1,
    alt3.plankton_stable =
      alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(
        alt1.plankton == 3,
        FALSE,
        TRUE
      ),
    alt2.plankton_nobloom =
      ifelse(
        alt2.plankton == 3,
        FALSE,
        TRUE
      ),
    alt3.plankton_nobloom =
      ifelse(
        alt3.plankton == 3,
        FALSE,
        TRUE
      ),
    alt1.mpa_fully =
      alt1.mpa == 1,
    alt2.mpa_fully =
      alt2.mpa == 1,
    alt3.mpa_fully =
      alt3.mpa == 1,
    alt1.mpa_highly =
      alt1.mpa == 2,
    alt2.mpa_highly =
      alt2.mpa == 2,
    alt3.mpa_highly =
      alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[5:ncol(dce_data)]

sorted_cols <-
  sort(cols_to_sort)

dce_data <-
  dce_data[
    ,
    c(
      names(dce_data)[1:4],
      sorted_cols
    )
  ]

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)

database <- database %>%
  mutate(across(
    ends_with("cost"),
    ~case_when(
      . == 10 & Country == "Poland" ~ 11,
      . == 20 & Country == "Poland" ~ 22,
      . == 40 & Country == "Poland" ~ 41,
      . == 80 & Country == "Poland" ~ 85,
      . == 120 & Country == "Poland" ~ 125,
      . == 180 & Country == "Poland" ~ 188,
      TRUE ~ .
    )
  ))

database <- database %>%
  filter(
    pref1 >= 1 & pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(list=ls()[ls() != "database"])

database <- database %>%
  mutate(
    Group = case_when(
      Country %in%
        c(
          "Poland",
          "France",
          "Germany"
        ) ~ "Central",
      Country %in%
        c(
          "Italy",
          "Basque"
        ) ~ "South",
      TRUE ~ NA_character_
    )
  )

### Central Europe
apollo_initialise()

apollo_control = list(
  modelName = "Clogit central vs south",
  modelDescr =
    "Conditional Logit in WTP-space",
  indivID = "RID",
  mixing = FALSE
)

apollo_beta = c(
  b0 = -77.3914197911875,
  b10 = 56.3067374160443,
  b11 = 12.4180332540675,
  b2 = 0.286408690011223,
  b30 = 32.9563336136676,
  b31 = 38.0346435809518,
  b4 = -0.0150420382306844
)

apollo_fixed = c()

apollo_inputs =
  apollo_validateInputs()

apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(
    apollo_beta,
    apollo_inputs
  )

  on.exit(
    apollo_detach(
      apollo_beta,
      apollo_inputs
    )
  )

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(
        alt1 = 1,
        alt2 = 2,
        alt3 = 3
      ),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Group == "Central"
  )

  P[['model']] =
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model = apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

model_central <- model

est_central <-
  as.data.frame(
    apollo_modelOutput(
      model,
      modelOutput_settings =
        list(printPVal = 1)
    )
  )

# Add one-sided p-values based on robust standard errors.
est_central[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_central[["Estimate"]] / est_central[["Rob.s.e."]]))

saveRDS(
  model_central,
  "./outputs/paper_A/models/08_central_conditional_logit.rds"
)

est_central_output <-
  est_central

est_central_output$Parameter <-
  rownames(est_central_output)

rownames(est_central_output) <-
  NULL

est_central_output <-
  est_central_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(est_central_output),
        "Parameter"
      )
    )
  ]

write.csv(
  est_central_output,
  "./outputs/paper_A/tables/08_central_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <-
  model$estimate

p_values <-
  est_central[["Rob.p(1-sided)"]]

std_errors <-
  est_central[["Rob.s.e."]]

results_df_central <-
  data.frame(
    Estimate = coef_estimates,
    Std_Error = std_errors,
    p_value = p_values
  )

results_df_central$CI_Lower <-
  results_df_central$Estimate -
  1.96 *
  results_df_central$Std_Error

results_df_central$CI_Upper <-
  results_df_central$Estimate +
  1.96 *
  results_df_central$Std_Error

results_df_central$Parameter <-
  rownames(results_df_central)

results_df_central <-
  results_df_central[
    ,
    c(
      "Parameter",
      setdiff(
        names(results_df_central),
        "Parameter"
      )
    )
  ]

rownames(results_df_central) <-
  NULL

print(results_df_central)


### Southern Europe
apollo_probabilities = function(
  apollo_beta,
  apollo_inputs,
  functionality = "estimate"
) {

  apollo_attach(
    apollo_beta,
    apollo_inputs
  )

  on.exit(
    apollo_detach(
      apollo_beta,
      apollo_inputs
    )
  )

  P = list()
  V = list()

  V[['alt1']] =
    -b4 * (
      b10 * alt1.plankton_stable +
        b11 * alt1.plankton_nobloom +
        b2 * alt1.carbon +
        b30 * alt1.mpa_fully +
        b31 * alt1.mpa_highly -
        alt1.cost
    )

  V[['alt2']] =
    -b4 * (
      b10 * alt2.plankton_stable +
        b11 * alt2.plankton_nobloom +
        b2 * alt2.carbon +
        b30 * alt2.mpa_fully +
        b31 * alt2.mpa_highly -
        alt2.cost
    )

  V[['alt3']] =
    -b4 * (
      b0 +
        b10 * alt3.plankton_stable +
        b11 * alt3.plankton_nobloom +
        b2 * alt3.carbon +
        b30 * alt3.mpa_fully +
        b31 * alt3.mpa_highly -
        alt3.cost
    )

  mnl_settings = list(
    alternatives =
      c(
        alt1 = 1,
        alt2 = 2,
        alt3 = 3
      ),
    avail = 1,
    choiceVar = pref1,
    V = V,
    rows = Group == "South"
  )

  P[['model']] =
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P =
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P =
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}

model =
  apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
  )

model_south <- model

est_south <-
  as.data.frame(
    apollo_modelOutput(
      model,
      modelOutput_settings =
        list(printPVal = 1)
    )
  )

# Add one-sided p-values based on robust standard errors.
est_south[["Rob.p(1-sided)"]] <-
  pnorm(-abs(est_south[["Estimate"]] / est_south[["Rob.s.e."]]))

saveRDS(
  model_south,
  "./outputs/paper_A/models/08_south_conditional_logit.rds"
)

est_south_output <-
  est_south

est_south_output$Parameter <-
  rownames(est_south_output)

rownames(est_south_output) <-
  NULL

est_south_output <-
  est_south_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(est_south_output),
        "Parameter"
      )
    )
  ]

write.csv(
  est_south_output,
  "./outputs/paper_A/tables/08_south_conditional_logit.csv",
  row.names = FALSE
)

coef_estimates <-
  model$estimate

p_values <-
  est_south[["Rob.p(1-sided)"]]

std_errors <-
  est_south[["Rob.s.e."]]

results_df_south <-
  data.frame(
    Estimate = coef_estimates,
    Std_Error = std_errors,
    p_value = p_values
  )

results_df_south$CI_Lower <-
  results_df_south$Estimate -
  1.96 *
  results_df_south$Std_Error

results_df_south$CI_Upper <-
  results_df_south$Estimate +
  1.96 *
  results_df_south$Std_Error

results_df_south$Parameter <-
  rownames(results_df_south)

results_df_south <-
  results_df_south[
    ,
    c(
      "Parameter",
      setdiff(
        names(results_df_south),
        "Parameter"
      )
    )
  ]

rownames(results_df_south) <-
  NULL

print(results_df_south)

combined_results <-
  rbind(
    transform(
      results_df_south,
      Group = "South"
    ),
    transform(
      results_df_central,
      Group = "Central"
    )
  )

write.csv(
  combined_results,
  "./outputs/paper_A/tables/08_central_south_comparison.csv",
  row.names = FALSE
)

rm(
  list = ls()[
    !ls() %in%
      c(
        "combined_results",
        "model_central",
        "model_south",
        "database"
      )
  ]
)


######################
# 9. Conditional logit model in WTP-space: interaction with income
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)

dce_data <-
  read.csv2(
    "./data/DMV data.csv"
  )

design <-
  read_delim(
    "./data/final design_baesyian efficient design with interaction.NGD",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_select =
      c(
        -Design,
        -starts_with("...")
      ),
    name_repair =
      "universal",
  n_max = 16
  ) %>%
  filter(
    !is.na(
      Choice.situation
    )
  )

nsets <-
  nrow(design)

design <-
  design %>%
  rename(
    choiceset =
      Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <-
  design %>%
  relocate(
    Block,
    .before =
      choiceset
  )

dce_answers <-
  dce_data[,44:59]

other_info <-
  dce_data[,-(44:59)]

other_info <-
  other_info %>%
  rename(
    RID = ID
  )

other_info <-
  other_info %>%
  replace_with_na_at(
    .vars =
      c(2:75),
    condition =
      ~.x == 99
  )

colnames(dce_answers) <-
  1:16

dce_answers <-
  tibble::rowid_to_column(
    dce_answers,
    "RID"
  )

choi <-
  dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to =
      "choiceset",
    values_to =
      "pref1"
  )

choi$choiceset <-
  as.numeric(
    choi$choiceset
  )

choi <-
  choi %>%
  mutate(
    pref1 =
      str_replace(
        pref1,
        "A",
        ""
      )
  )

dce_data <-
  merge(
    choi,
    design,
    by =
      "choiceset",
    all =
      TRUE
  )

dce_data <-
  dce_data %>%
  relocate(RID)

dce_data <-
  dce_data[
    order(
      dce_data[,1],
      dce_data[,2]
    ),
  ]

row.names(dce_data) <-
  NULL

dce_data <-
  dce_data %>%
  mutate(
    alt1.plankton_stable =
      alt1.plankton == 1,
    alt2.plankton_stable =
      alt2.plankton == 1,
    alt3.plankton_stable =
      alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(
        alt1.plankton == 3,
        FALSE,
        TRUE
      ),
    alt2.plankton_nobloom =
      ifelse(
        alt2.plankton == 3,
        FALSE,
        TRUE
      ),
    alt3.plankton_nobloom =
      ifelse(
        alt3.plankton == 3,
        FALSE,
        TRUE
      ),
    alt1.mpa_fully =
      alt1.mpa == 1,
    alt2.mpa_fully =
      alt2.mpa == 1,
    alt3.mpa_fully =
      alt3.mpa == 1,
    alt1.mpa_highly =
      alt1.mpa == 2,
    alt2.mpa_highly =
      alt2.mpa == 2,
    alt3.mpa_highly =
      alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[
    5:ncol(dce_data)
  ]

sorted_cols <-
  sort(cols_to_sort)

dce_data <-
  dce_data[
    ,
    c(
      names(dce_data)[1:4],
      sorted_cols
    )
  ]

database <-
  merge(
    dce_data,
    other_info,
    by = "RID"
  )

database <-
  database %>%
  mutate(
    across(
      ends_with("cost"),
      ~case_when(
        . == 10 &
          Country == "Poland" ~ 11,
        . == 20 &
          Country == "Poland" ~ 22,
        . == 40 &
          Country == "Poland" ~ 41,
        . == 80 &
          Country == "Poland" ~ 85,
        . == 120 &
          Country == "Poland" ~ 125,
        . == 180 &
          Country == "Poland" ~ 188,
        TRUE ~ .
      )
    )
  )

database <-
  database %>%
  filter(
    pref1 >= 1 &
      pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84,
    RID != 129
  )

rm(
  list =
    ls()[
      ls() != "database"
    ]
)

apollo_initialise()

apollo_control =
  list(
    modelName =
      "Clogit interaction",
    modelDescr =
      "Conditional Logit in WTP-space - interaction income",
    indivID =
      "RID",
    mixing =
      FALSE
  )

apollo_beta =
  c(
    b0 =
      -77.3914197911875,
    b10 =
      56.3067374160443,
    b11 =
      12.4180332540675,
    b2 =
      0.286408690011223,
    b30 =
      32.9563336136676,
    b31 =
      38.0346435809518,
    b4 =
      -0.0150420382306844,
    bi10 = 0,
    bi11 = 0,
    bi2 = 0,
    bi30 = 0,
    bi31 = 0,
    bi4 = 0,
    bi0 = 0
  )

apollo_fixed =
  c()

apollo_inputs =
  apollo_validateInputs()

apollo_probabilities =
  function(
    apollo_beta,
    apollo_inputs,
    functionality =
      "estimate"
  ) {

    apollo_attach(
      apollo_beta,
      apollo_inputs
    )

    on.exit(
      apollo_detach(
        apollo_beta,
        apollo_inputs
      )
    )

    P =
      list()

    V =
      list()

    V[['alt1']] =
      -(
        b4 +
          bi4 *
          `Net.income.group`
      ) * (
        b10 *
          alt1.plankton_stable +
          b11 *
          alt1.plankton_nobloom +
          b2 *
          alt1.carbon +
          b30 *
          alt1.mpa_fully +
          b31 *
          alt1.mpa_highly -
          alt1.cost +
          bi10 *
          alt1.plankton_stable *
          `Net.income.group` +
          bi11 *
          alt1.plankton_nobloom *
          `Net.income.group` +
          bi2 *
          `Net.income.group` *
          alt1.carbon +
          bi30 *
          alt1.mpa_fully *
          `Net.income.group` +
          bi31 *
          alt1.mpa_highly *
          `Net.income.group`
      )

    V[['alt2']] =
      -(
        b4 +
          bi4 *
          `Net.income.group`
      ) * (
        b10 *
          alt2.plankton_stable +
          b11 *
          alt2.plankton_nobloom +
          b2 *
          alt2.carbon +
          b30 *
          alt2.mpa_fully +
          b31 *
          alt2.mpa_highly -
          alt2.cost +
          bi10 *
          alt2.plankton_stable *
          `Net.income.group` +
          bi11 *
          alt2.plankton_nobloom *
          `Net.income.group` +
          bi2 *
          `Net.income.group` *
          alt2.carbon +
          bi30 *
          alt2.mpa_fully *
          `Net.income.group` +
          bi31 *
          alt2.mpa_highly *
          `Net.income.group`
      )

    V[['alt3']] =
      -(
        b4 +
          bi4 *
          `Net.income.group`
      ) * (
        b0 +
          bi0 *
          `Net.income.group` +
          b10 *
          alt3.plankton_stable +
          b11 *
          alt3.plankton_nobloom +
          b2 *
          alt3.carbon +
          b30 *
          alt3.mpa_fully +
          b31 *
          alt3.mpa_highly -
          alt3.cost +
          bi10 *
          alt3.plankton_stable *
          `Net.income.group` +
          bi11 *
          alt3.plankton_nobloom *
          `Net.income.group` +
          bi2 *
          `Net.income.group` *
          alt3.carbon +
          bi30 *
          alt3.mpa_fully *
          `Net.income.group` +
          bi31 *
          alt3.mpa_highly *
          `Net.income.group`
      )

    mnl_settings =
      list(
        alternatives =
          c(
            alt1 = 1,
            alt2 = 2,
            alt3 = 3
          ),
        avail = 1,
        choiceVar =
          pref1,
        V = V
      )

    P[['model']] =
      apollo_mnl(
        mnl_settings,
        functionality
      )

    P =
      apollo_panelProd(
        P,
        apollo_inputs,
        functionality
      )

    P =
      apollo_prepareProb(
        P,
        apollo_inputs,
        functionality
      )

    return(P)
  }

model =
  apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
  )

model_output <-
  as.data.frame(
    apollo_modelOutput(
      model,
      modelOutput_settings =
        list(printPVal = 1)
    )
  )

# Add one-sided p-values based on robust standard errors.
model_output[["Rob.p(1-sided)"]] <-
  pnorm(-abs(model_output[["Estimate"]] / model_output[["Rob.s.e."]]))

saveRDS(
  model,
  "./outputs/paper_A/models/09_income_interaction.rds"
)

model_output$Parameter <-
  rownames(model_output)

rownames(model_output) <-
  NULL

model_output <-
  model_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(model_output),
        "Parameter"
      )
    )
  ]

write.csv(
  model_output,
  "./outputs/paper_A/tables/09_income_interaction.csv",
  row.names = FALSE
)

rm(
  list =
    ls()[
      !ls() %in%
        c(
          "model",
          "database"
        )
    ]
)


############################
# Optional starting-value searches
############################

# These searches document how the starting values used above were obtained.
# They are computationally intensive and are therefore skipped by default.
# Set to TRUE to reproduce the starting-value searches.

run_starting_value_search <- FALSE

if (run_starting_value_search) {


######################
# 10. Documentation of starting value search
#     (Conditional logit model in WTP-space)
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)

dce_data <-
  read.csv2(
    "./data/DMV data.csv"
  )

design <-
  read_delim(
    "./data/final design_baesyian efficient design with interaction.NGD",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_select =
      c(
        -Design,
        -starts_with("...")
      ),
    name_repair =
      "universal",
  n_max = 16
  ) %>%
  filter(
    !is.na(
      Choice.situation
    )
  )

nsets <-
  nrow(design)

design <-
  design %>%
  rename(
    choiceset =
      Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <-
  design %>%
  relocate(
    Block,
    .before = choiceset
  )

dce_answers <-
  dce_data[,44:59]

other_info <-
  dce_data[,-(44:59)]

other_info <-
  other_info %>%
  rename(
    RID = ID
  )

other_info <-
  other_info %>%
  replace_with_na_at(
    .vars =
      c(2:75),
    condition =
      ~.x == 99
  )

colnames(dce_answers) <-
  1:16

dce_answers <-
  tibble::rowid_to_column(
    dce_answers,
    "RID"
  )

choi <-
  dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to =
      "choiceset",
    values_to =
      "pref1"
  )

choi$choiceset <-
  as.numeric(
    choi$choiceset
  )

choi <-
  choi %>%
  mutate(
    pref1 =
      str_replace(
        pref1,
        "A",
        ""
      )
  )

dce_data <-
  merge(
    choi,
    design,
    by =
      "choiceset",
    all =
      TRUE
  )

dce_data <-
  dce_data %>%
  relocate(RID)

dce_data <-
  dce_data[
    order(
      dce_data[,1],
      dce_data[,2]
    ),
  ]

row.names(dce_data) <-
  NULL

dce_data <-
  dce_data %>%
  mutate(
    alt1.plankton_stable =
      alt1.plankton == 1,
    alt2.plankton_stable =
      alt2.plankton == 1,
    alt3.plankton_stable =
      alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(
        alt1.plankton == 3,
        FALSE,
        TRUE
      ),
    alt2.plankton_nobloom =
      ifelse(
        alt2.plankton == 3,
        FALSE,
        TRUE
      ),
    alt3.plankton_nobloom =
      ifelse(
        alt3.plankton == 3,
        FALSE,
        TRUE
      ),
    alt1.mpa_fully =
      alt1.mpa == 1,
    alt2.mpa_fully =
      alt2.mpa == 1,
    alt3.mpa_fully =
      alt3.mpa == 1,
    alt1.mpa_highly =
      alt1.mpa == 2,
    alt2.mpa_highly =
      alt2.mpa == 2,
    alt3.mpa_highly =
      alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[
    5:ncol(dce_data)
  ]

sorted_cols <-
  sort(
    cols_to_sort
  )

dce_data <-
  dce_data[
    ,
    c(
      names(dce_data)[1:4],
      sorted_cols
    )
  ]

database <-
  merge(
    dce_data,
    other_info,
    by =
      "RID",
    all =
      TRUE
  )

database <-
  database %>%
  mutate(
    across(
      ends_with("cost"),
      ~case_when(
        . == 10 &
          Country == "Poland" ~ 11,
        . == 20 &
          Country == "Poland" ~ 22,
        . == 40 &
          Country == "Poland" ~ 41,
        . == 80 &
          Country == "Poland" ~ 85,
        . == 120 &
          Country == "Poland" ~ 125,
        . == 180 &
          Country == "Poland" ~ 188,
        TRUE ~ .
      )
    )
  )

database <-
  database %>%
  filter(
    pref1 >= 1 &
      pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(
  list =
    ls()[
      ls() != "database"
    ]
)

apollo_initialise()

apollo_control =
  list(
    modelName =
      "Clogit PLN",
    modelDescr =
      "Conditional Logit in WTP-space",
    indivID =
      "RID",
    mixing =
      FALSE
  )

apollo_beta =
  c(
    b0 = 0,
    b10 = 0,
    b11 = 0,
    b2 = 0,
    b30 = 0,
    b31 = 0,
    b4 = 0
  )

apollo_fixed =
  c()

apollo_inputs =
  apollo_validateInputs()

apollo_probabilities =
  function(
    apollo_beta,
    apollo_inputs,
    functionality =
      "estimate"
  ) {

    apollo_attach(
      apollo_beta,
      apollo_inputs
    )

    on.exit(
      apollo_detach(
        apollo_beta,
        apollo_inputs
      )
    )

    P =
      list()

    V =
      list()

    V[['alt1']] =
      -b4 * (
        b10 *
          alt1.plankton_stable +
          b11 *
          alt1.plankton_nobloom +
          b2 *
          alt1.carbon +
          b30 *
          alt1.mpa_fully +
          b31 *
          alt1.mpa_highly -
          alt1.cost
      )

    V[['alt2']] =
      -b4 * (
        b10 *
          alt2.plankton_stable +
          b11 *
          alt2.plankton_nobloom +
          b2 *
          alt2.carbon +
          b30 *
          alt2.mpa_fully +
          b31 *
          alt2.mpa_highly -
          alt2.cost
      )

    V[['alt3']] =
      -b4 * (
        b0 +
          b10 *
          alt3.plankton_stable +
          b11 *
          alt3.plankton_nobloom +
          b2 *
          alt3.carbon +
          b30 *
          alt3.mpa_fully +
          b31 *
          alt3.mpa_highly -
          alt3.cost
      )

    mnl_settings =
      list(
        alternatives =
          c(
            alt1 = 1,
            alt2 = 2,
            alt3 = 3
          ),
        avail = 1,
        choiceVar =
          pref1,
        V = V
      )

    P[['model']] =
      apollo_mnl(
        mnl_settings,
        functionality
      )

    P =
      apollo_panelProd(
        P,
        apollo_inputs,
        functionality
      )

    P =
      apollo_prepareProb(
        P,
        apollo_inputs,
        functionality
      )

    return(P)
  }

model =
  apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
  )

searchStart_settings =
  list(
    searchRoutine =
      "gss",
    nSearches =
      10,
    nSearchDraws =
      1000,
    searchTol =
      1e-6,
    searchMaxEval =
      1000,
    searchMaxIter =
      100,
    searchParallel =
      FALSE
  )

apollo_beta =
  apollo_searchStart(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs,
    searchStart_settings
  )

starting_values <-
  data.frame(
    Parameter =
      names(apollo_beta),
    Starting_value =
      unname(apollo_beta)
  )

write.csv(
  starting_values,
  "./outputs/paper_A/tables/10_conditional_logit_starting_values.csv",
  row.names = FALSE
)

# Best candidate (LL=-2229.2)
#       Value
# b0  -77.3839
# b10  56.3098
# b11  12.4279
# b2    0.2863
# b30  32.9715
# b31  38.0468
# b4   -0.0150
#
# Estimating two different conditional logit models with the above
# starting values and with initialising all parameters at zero yielded
# the same results.


######################
# 11. Documentation of starting value search
#     (Mixed logit model in WTP-space)
######################

### Place the following files in the data/ directory:
### 1) "final design_baesyian efficient design with interaction.NGD"
### 2) "DMV data.csv"
### --> Source: https://zenodo.org/records/12638010

rm(list=ls())
library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)

dce_data <-
  read.csv2(
    "./data/DMV data.csv"
  )

design <-
  read_delim(
    "./data/final design_baesyian efficient design with interaction.NGD",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_select =
      c(
        -Design,
        -starts_with("...")
      ),
    name_repair =
      "universal",
  n_max = 16
  ) %>%
  filter(
    !is.na(
      Choice.situation
    )
  )

nsets <-
  nrow(design)

design <-
  design %>%
  rename(
    choiceset =
      Choice.situation
  )

design$alt3.plankton = 3
design$alt3.carbon = 0
design$alt3.mpa = 3
design$alt3.cost = 0

design <-
  design %>%
  relocate(
    Block,
    .before =
      choiceset
  )

dce_answers <-
  dce_data[,44:59]

other_info <-
  dce_data[,-(44:59)]

other_info <-
  other_info %>%
  rename(
    RID = ID
  )

other_info <-
  other_info %>%
  replace_with_na_at(
    .vars =
      c(2:75),
    condition =
      ~.x == 99
  )

colnames(dce_answers) <-
  1:16

dce_answers <-
  tibble::rowid_to_column(
    dce_answers,
    "RID"
  )

choi <-
  dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to =
      "choiceset",
    values_to =
      "pref1"
  )

choi$choiceset <-
  as.numeric(
    choi$choiceset
  )

choi <-
  choi %>%
  mutate(
    pref1 =
      str_replace(
        pref1,
        "A",
        ""
      )
  )

dce_data <-
  merge(
    choi,
    design,
    by =
      "choiceset",
    all =
      TRUE
  )

dce_data <-
  dce_data %>%
  relocate(RID)

dce_data <-
  dce_data[
    order(
      dce_data[,1],
      dce_data[,2]
    ),
  ]

row.names(dce_data) <-
  NULL

dce_data <-
  dce_data %>%
  mutate(
    alt1.plankton_stable =
      alt1.plankton == 1,
    alt2.plankton_stable =
      alt2.plankton == 1,
    alt3.plankton_stable =
      alt3.plankton == 1,
    alt1.plankton_nobloom =
      ifelse(
        alt1.plankton == 3,
        FALSE,
        TRUE
      ),
    alt2.plankton_nobloom =
      ifelse(
        alt2.plankton == 3,
        FALSE,
        TRUE
      ),
    alt3.plankton_nobloom =
      ifelse(
        alt3.plankton == 3,
        FALSE,
        TRUE
      ),
    alt1.mpa_fully =
      alt1.mpa == 1,
    alt2.mpa_fully =
      alt2.mpa == 1,
    alt3.mpa_fully =
      alt3.mpa == 1,
    alt1.mpa_highly =
      alt1.mpa == 2,
    alt2.mpa_highly =
      alt2.mpa == 2,
    alt3.mpa_highly =
      alt3.mpa == 2
  )

cols_to_sort <-
  names(dce_data)[
    5:ncol(dce_data)
  ]

sorted_cols <-
  sort(cols_to_sort)

dce_data <-
  dce_data[
    ,
    c(
      names(dce_data)[1:4],
      sorted_cols
    )
  ]

database <-
  merge(
    dce_data,
    other_info,
    by =
      "RID",
    all =
      TRUE
  )

database <-
  database %>%
  mutate(
    across(
      ends_with("cost"),
      ~case_when(
        . == 10 &
          Country == "Poland" ~ 11,
        . == 20 &
          Country == "Poland" ~ 22,
        . == 40 &
          Country == "Poland" ~ 41,
        . == 80 &
          Country == "Poland" ~ 85,
        . == 120 &
          Country == "Poland" ~ 125,
        . == 180 &
          Country == "Poland" ~ 188,
        TRUE ~ .
      )
    )
  )

database <-
  database %>%
  filter(
    pref1 >= 1 &
      pref1 <= 3,
    RID != 12,
    RID != 61,
    RID != 62,
    RID != 84
  )

rm(
  list =
    ls()[
      ls() != "database"
    ]
)

apollo_initialise()

apollo_control =
  list(
    modelName =
      "Mixed Logit",
    modelDescr =
      "Mixed Logit in WTP-space",
    indivID =
      "RID",
    mixing =
      TRUE,
    nCores =
      5
  )

apollo_beta =
  c(
    mu_log_b0 =
      -87.93,
    sigma_log_b0 =
      67.90,
    mu_log_b10 =
      51.28,
    sigma_log_b10 =
      48.62,
    mu_log_b11 =
      15.58,
    sigma_log_b11 =
      41.96,
    mu_log_b2 =
      0.30,
    sigma_log_b2 =
      0.63,
    mu_log_b30 =
      34.04,
    sigma_log_b30 =
      27.18,
    mu_log_b31 =
      31.98,
    sigma_log_b31 =
      10.84,
    mu_log_b4 =
      -0.5,
    sigma_log_b4 =
      0.2
  )

# The starting values for the non-price attributes were retrieved
# from a starting value search in preference-space, while the
# cost-parameter values were chosen through a trial-and-error approach.

apollo_fixed =
  c()

apollo_draws =
  list(
    interDrawsType =
      "mlhs",
    interNDraws =
      500,
    interNormDraws =
      c(
        "draws_0",
        "draws_1",
        "draws_2",
        "draws_3",
        "draws_4",
        "draws_5",
        "draws_6"
      )
  )

apollo_randCoeff =
  function(
    apollo_beta,
    apollo_inputs
  ) {

    randcoeff =
      list()

    randcoeff[["b0"]] =
      mu_log_b0 +
      sigma_log_b0 *
      draws_0

    randcoeff[["b10"]] =
      mu_log_b10 +
      sigma_log_b10 *
      draws_1

    randcoeff[["b11"]] =
      mu_log_b11 +
      sigma_log_b11 *
      draws_2

    randcoeff[["b2"]] =
      mu_log_b2 +
      sigma_log_b2 *
      draws_3

    randcoeff[["b30"]] =
      mu_log_b30 +
      sigma_log_b30 *
      draws_4

    randcoeff[["b31"]] =
      mu_log_b31 +
      sigma_log_b31 *
      draws_5

    randcoeff[["b4"]] =
      -exp(
        mu_log_b4 +
        sigma_log_b4 *
        draws_6
      )

    return(randcoeff)
  }

apollo_inputs =
  apollo_validateInputs()

apollo_probabilities =
  function(
    apollo_beta,
    apollo_inputs,
    functionality =
      "estimate"
  ) {

    apollo_attach(
      apollo_beta,
      apollo_inputs
    )

    on.exit(
      apollo_detach(
        apollo_beta,
        apollo_inputs
      )
    )

    P =
      list()

    V =
      list()

    V[['alt1']] =
      -b4 * (
        b10 *
          alt1.plankton_stable +
          b11 *
          alt1.plankton_nobloom +
          b2 *
          alt1.carbon +
          b30 *
          alt1.mpa_fully +
          b31 *
          alt1.mpa_highly -
          alt1.cost
      )

    V[['alt2']] =
      -b4 * (
        b10 *
          alt2.plankton_stable +
          b11 *
          alt2.plankton_nobloom +
          b2 *
          alt2.carbon +
          b30 *
          alt2.mpa_fully +
          b31 *
          alt2.mpa_highly -
          alt2.cost
      )

    V[['alt3']] =
      -b4 * (
        b0 +
          b10 *
          alt3.plankton_stable +
          b11 *
          alt3.plankton_nobloom +
          b2 *
          alt3.carbon +
          b30 *
          alt3.mpa_fully +
          b31 *
          alt3.mpa_highly -
          alt3.cost
      )

    mnl_settings =
      list(
        alternatives =
          c(
            alt1 = 1,
            alt2 = 2,
            alt3 = 3
          ),
        avail = 1,
        choiceVar =
          pref1,
        V = V
      )

    P[['model']] =
      apollo_mnl(
        mnl_settings,
        functionality
      )

    P =
      apollo_panelProd(
        P,
        apollo_inputs,
        functionality
      )

    P =
      apollo_avgInterDraws(
        P,
        apollo_inputs,
        functionality
      )

    P =
      apollo_prepareProb(
        P,
        apollo_inputs,
        functionality
      )

    return(P)
  }

model =
  apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
  )

searchStart_settings =
  list(
    searchRoutine =
      "gss",
    nSearches =
      10,
    nSearchDraws =
      1000,
    searchTol =
      1e-6,
    searchMaxEval =
      1000,
    searchMaxIter =
      100,
    searchParallel =
      FALSE
  )

apollo_beta =
  apollo_searchStart(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs,
    searchStart_settings
  )

starting_values <-
  data.frame(
    Parameter =
      names(apollo_beta),
    Starting_value =
      unname(apollo_beta)
  )

write.csv(
  starting_values,
  "./outputs/paper_A/tables/11_mixed_logit_starting_values.csv",
  row.names = FALSE
)

} # End optional starting-value searches


# Best candidate (LL=-1883.476)
#                 Value
# mu_log_b0       -88.0667
# sigma_log_b0     54.8022
# mu_log_b10       48.5357
# sigma_log_b10    46.7755
# mu_log_b11       15.2579
# sigma_log_b11   -45.4263
# mu_log_b2         0.2546
# sigma_log_b2      0.6713
# mu_log_b30       35.6045
# sigma_log_b30   -29.5274
# mu_log_b31       31.1781
# sigma_log_b31   -16.6907
# mu_log_b4        -3.6243
# sigma_log_b4     -0.3534
#
# This candidate was determined to be the most suitable following
# a comparison of the results obtained from an LL comparison with
# other candidates from different searches that employed the same
# code as above, but with varying start values. The first search
# started with all beta values set to zero and assuming all
# parameters (including cost) to be normally distributed.
# Subsequently, the best candidate was chosen as a basis for the
# next search, and so on.
