############################
# Analysis code for Paper C of the cumulative dissertation
# Quasi-binomial, conditional logit, and ordered logistic regression models
# Input-data sources and download instructions: see data/README.md
############################

# Run this script from the repository root:
# source("R/03_paper_C.R")
#
# Alternatively, from a terminal:
# Rscript R/03_paper_C.R
#
# Relative paths such as "./data/..." assume that the repository root
# is the current working directory.

############## Overview ###############
# 1. Quasi-binomial logit regression models of policy support
# 2. Conditional logit model estimates in willingness-to-pay (WTP) space
# 3. Ordered logistic regression results – estimating perceived influence of other workshop participants
# 4. Documentation of starting-value search for conditional logit models in WTP-space
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
  "./outputs/paper_C/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "./outputs/paper_C/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# 1. Quasi-binomial logit regression models of policy support
############################

# Required files in the data/ directory:
# 1) DMV data_EN.csv
#    Source: https://zenodo.org/records/17710200
# 2) Pro_contra_statements_respondents_csv_UTF_8.csv
#    Source: supplementary materials

rm(list = ls())

library(dplyr)
library(tidyr)


# Helper function for extracting GLM coefficient tables.
extract_glm_coefficients <- function(models) {

  bind_rows(
    lapply(
      seq_along(models),
      function(i) {

        model_summary <- summary(models[[i]])
        coef_matrix <- coef(model_summary)

        data.frame(
          Model = names(models)[i],
          Parameter = rownames(coef_matrix),
          Estimate = coef_matrix[,1],
          Std_Error = coef_matrix[,2],
          Statistic = coef_matrix[,3],
          p_value = coef_matrix[,4],
          row.names = NULL
        )
      }
    )
  )
}


# Read main survey and DCE data.
dce_data <- read.csv2(
  "./data/DMV data_EN.csv"
) %>%
  rename(
    RID = ID
  )


# Calculate E-PVQ value scores and grand-mean-centred value scores.
dce_data <- dce_data %>%
  mutate(
    egoistic = round(
      rowMeans(
        across(
          c(
            Q4.1,
            Q4.3,
            Q4.6,
            Q4.10,
            Q4.15
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    egoistic_c =
      egoistic -
      mean(
        egoistic,
        na.rm = TRUE
      ),

    hedonic = round(
      rowMeans(
        across(
          c(
            Q4.2,
            Q4.5,
            Q4.11
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    hedonic_c =
      hedonic -
      mean(
        hedonic,
        na.rm = TRUE
      ),

    altruistic = round(
      rowMeans(
        across(
          c(
            Q4.4,
            Q4.7,
            Q4.8,
            Q4.13,
            Q4.17
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    altruistic_c =
      altruistic -
      mean(
        altruistic,
        na.rm = TRUE
      ),

    biospheric = round(
      rowMeans(
        across(
          c(
            Q4.9,
            Q4.12,
            Q4.14,
            Q4.16
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    biospheric_c =
      biospheric -
      mean(
        biospheric,
        na.rm = TRUE
      )
  )


# Read respondent-level coded pro/contra statements and aggregate them by respondent.
ind_stat <- read.csv2(
  "./data/Pro_contra_statements_respondents_csv_UTF_8.csv",
  fileEncoding = "UTF-8-BOM"
) %>%
  filter(!is.na(RID))

ind_stat_merged <- ind_stat %>%
  group_by(RID) %>%
  summarise(
    Pro = sum(
      Pro,
      na.rm = TRUE
    ),
    Contra = sum(
      Contra,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join statement totals to the main data and calculate the leave-one-out
# pro-policy tendency of each respondent's workshop discussion environment.
dce_data <- dce_data %>%
  left_join(
    ind_stat_merged,
    by = "RID"
  ) %>%
  mutate(
    Pro = replace_na(
      Pro,
      0
    ),
    Contra = replace_na(
      Contra,
      0
    ),
    net_i =
      Pro -
      Contra
  ) %>%
  group_by(WS_ID) %>%
  mutate(
    WS_policy_tend =
      sum(
        Pro,
        na.rm = TRUE
      ) -
      sum(
        Contra,
        na.rm = TRUE
      ),

    LOO_tend =
      WS_policy_tend -
      net_i
  ) %>%
  ungroup()


# Create income group variable and policy-support outcome.
# BAU_not_chosen counts the number of choice situations in which the respondent
# selected a policy alternative rather than the business-as-usual alternative.
choice_vars <- paste0(
  "CS",
  1:16
)

stopifnot(
  all(
    choice_vars %in%
      names(dce_data)
  )
)

dce_data <- dce_data %>%
  mutate(
    `Net income group_n` = case_when(
      is.na(`Net.income.group`) ~ NA_real_,
      `Net.income.group` <= 4 ~ 1,
      `Net.income.group` <= 8 ~ 2,
      `Net.income.group` > 8 ~ 3
    ),

    `Net income group_n` =
      factor(
        `Net income group_n`,
        levels = c(
          1,
          2,
          3
        )
      ),

    BAU_not_chosen =
      rowSums(
        across(
          all_of(choice_vars),
          ~ .x < 3
        ),
        na.rm = TRUE
      )
  )

stopifnot(
  all(
    dce_data$BAU_not_chosen >= 0 &
      dce_data$BAU_not_chosen <= 16
  )
)


# Define model formulas.
model_formulas <- list(

  model1 =
    cbind(
      BAU_not_chosen,
      16 - BAU_not_chosen
    ) ~
    LOO_tend +
    biospheric_c +
    `Net income group_n`,

  model2 =
    cbind(
      BAU_not_chosen,
      16 - BAU_not_chosen
    ) ~
    LOO_tend +
    egoistic_c +
    `Net income group_n`,

  model3 =
    cbind(
      BAU_not_chosen,
      16 - BAU_not_chosen
    ) ~
    LOO_tend +
    hedonic_c +
    `Net income group_n`,

  model4 =
    cbind(
      BAU_not_chosen,
      16 - BAU_not_chosen
    ) ~
    LOO_tend +
    altruistic_c +
    `Net income group_n`,

  model5 =
    cbind(
      BAU_not_chosen,
      16 - BAU_not_chosen
    ) ~
    LOO_tend +
    biospheric_c +
    egoistic_c +
    hedonic_c +
    altruistic_c +
    `Net income group_n`
)


# Fit initial binomial models.
binom_models <- lapply(
  model_formulas,
  glm,
  data = dce_data,
  family = binomial(
    link = "logit"
  )
)

lapply(
  binom_models,
  summary
)


# Save binomial models used for the overdispersion assessment.
saveRDS(
  binom_models,
  "./outputs/paper_C/models/01_binomial_diagnostic_models.rds"
)


# Save binomial coefficient tables.
binom_coefficients <-
  extract_glm_coefficients(
    binom_models
  )

write.csv(
  binom_coefficients,
  "./outputs/paper_C/tables/01_binomial_diagnostic_coefficients.csv",
  row.names = FALSE
)


# Assess overdispersion in the binomial models.
#
# In a correctly specified binomial GLM, the Pearson chi-square statistic
# should be approximately equal to the residual degrees of freedom. The
# dispersion ratio therefore compares observed residual variation with the
# residual variation expected under the binomial model.
#
# The function returns:
# - chisq: Pearson chi-square statistic, calculated as the sum of squared
#          Pearson residuals.
# - ratio: Pearson chi-square statistic divided by residual degrees of freedom.
#          Values close to 1 suggest no strong overdispersion; values clearly
#          above 1 indicate extra-binomial variation.
# - df: residual degrees of freedom.
# - p: p-value from comparing the Pearson chi-square statistic with a
#      chi-square distribution with the corresponding residual degrees of
#      freedom. A small p-value suggests more residual variation than expected
#      under the standard binomial model.

overdisp_fun <- function(model) {

  # Extract Pearson residuals from the fitted binomial model.
  rp <- residuals(
    model,
    type = "pearson"
  )

  # Calculate the Pearson chi-square statistic.
  pearson_chisq <-
    sum(
      rp^2
    )

  # Calculate the dispersion ratio.
  ratio <-
    pearson_chisq /
    model$df.residual

  # Calculate the p-value for excess residual variation.
  p_value <- pchisq(
    pearson_chisq,
    df = model$df.residual,
    lower.tail = FALSE
  )

  # Return the diagnostic results as a named vector.
  c(
    chisq = pearson_chisq,
    ratio = ratio,
    df = model$df.residual,
    p = p_value
  )
}


overdisp_results <- lapply(
  binom_models,
  overdisp_fun
)

overdisp_results


# Convert overdispersion results into a persistent table.
overdisp_table <-
  as.data.frame(
    do.call(
      rbind,
      overdisp_results
    )
  )

overdisp_table$Model <-
  rownames(
    overdisp_table
  )

rownames(
  overdisp_table
) <- NULL

overdisp_table <-
  overdisp_table[
    ,
    c(
      "Model",
      "chisq",
      "ratio",
      "df",
      "p"
    )
  ]

write.csv(
  overdisp_table,
  "./outputs/paper_C/tables/01_overdispersion_diagnostics.csv",
  row.names = FALSE
)


# Fit final quasi-binomial models.
quasi_models <- lapply(
  model_formulas,
  glm,
  data = dce_data,
  family = quasibinomial(
    link = "logit"
  )
)

lapply(
  quasi_models,
  summary
)


# Save complete final quasi-binomial models.
saveRDS(
  quasi_models,
  "./outputs/paper_C/models/01_quasibinomial_policy_support_models.rds"
)


# Save coefficient tables for final quasi-binomial models.
quasi_coefficients <-
  extract_glm_coefficients(
    quasi_models
  )

write.csv(
  quasi_coefficients,
  "./outputs/paper_C/tables/01_quasibinomial_policy_support_coefficients.csv",
  row.names = FALSE
)


# Save model-fit and dispersion information.
quasi_fit <-
  bind_rows(
    lapply(
      seq_along(
        quasi_models
      ),
      function(i) {

        model <-
          quasi_models[[i]]

        data.frame(
          Model =
            names(
              quasi_models
            )[i],

          Observations =
            nobs(model),

          Residual_deviance =
            deviance(model),

          Residual_df =
            model$df.residual,

          Dispersion =
            summary(model)$dispersion
        )
      }
    )
  )

write.csv(
  quasi_fit,
  "./outputs/paper_C/tables/01_quasibinomial_policy_support_fit.csv",
  row.names = FALSE
)


############################
# 2. Conditional logit model estimates in willingness-to-pay (WTP) space
############################

# Required files in the data/ directory:
# 1) final design_baesyian efficient design with interaction.NGD
# 2) DMV data_EN.csv
#    Source: https://zenodo.org/records/17710200
# 3) Pro_contra_statements_respondents_csv_UTF_8.csv
#    Source: supplementary materials

rm(list = ls())

library(apollo)
library(dplyr)
library(readr)
library(naniar)
library(tidyr)
library(stringr)


# Read main survey and DCE data.
dce_data <- read.csv2(
  "./data/DMV data_EN.csv"
) %>%
  rename(
    RID = ID
  )


# Calculate E-PVQ value scores and grand-mean-centred value scores.
dce_data <- dce_data %>%
  mutate(
    egoistic = round(
      rowMeans(
        across(
          c(
            Q4.1,
            Q4.3,
            Q4.6,
            Q4.10,
            Q4.15
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    egoistic_c =
      egoistic -
      mean(
        egoistic,
        na.rm = TRUE
      ),

    hedonic = round(
      rowMeans(
        across(
          c(
            Q4.2,
            Q4.5,
            Q4.11
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    hedonic_c =
      hedonic -
      mean(
        hedonic,
        na.rm = TRUE
      ),

    altruistic = round(
      rowMeans(
        across(
          c(
            Q4.4,
            Q4.7,
            Q4.8,
            Q4.13,
            Q4.17
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    altruistic_c =
      altruistic -
      mean(
        altruistic,
        na.rm = TRUE
      ),

    biospheric = round(
      rowMeans(
        across(
          c(
            Q4.9,
            Q4.12,
            Q4.14,
            Q4.16
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    biospheric_c =
      biospheric -
      mean(
        biospheric,
        na.rm = TRUE
      )
  )


# Read respondent-level coded pro/contra statements and aggregate them by respondent.
ind_stat <- read.csv2(
  "./data/Pro_contra_statements_respondents_csv_UTF_8.csv",
  fileEncoding = "UTF-8-BOM"
) %>%
  filter(!is.na(RID))

ind_stat_merged <- ind_stat %>%
  group_by(RID) %>%
  summarise(
    Pro = sum(
      Pro,
      na.rm = TRUE
    ),
    Contra = sum(
      Contra,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join statement totals and calculate leave-one-out pro-policy tendency.
dce_data <- dce_data %>%
  left_join(
    ind_stat_merged,
    by = "RID"
  ) %>%
  mutate(
    Pro = replace_na(
      Pro,
      0
    ),
    Contra = replace_na(
      Contra,
      0
    ),
    net_i =
      Pro -
      Contra
  ) %>%
  group_by(WS_ID) %>%
  mutate(
    WS_policy_tend =
      sum(
        Pro,
        na.rm = TRUE
      ) -
      sum(
        Contra,
        na.rm = TRUE
      ),

    LOO_tend =
      WS_policy_tend -
      net_i
  ) %>%
  ungroup()


# Read experimental design and reshape DCE responses to long format.
design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(
    -Design,
    -starts_with("...")
  ),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(
    !is.na(
      Choice.situation
    )
  )

nsets <- nrow(
  design
)

design <- design %>%
  rename(
    choiceset =
      Choice.situation
  )

design$alt3.plankton <- 3
design$alt3.carbon <- 0
design$alt3.mpa <- 3
design$alt3.cost <- 0

design <- design %>%
  relocate(
    Block,
    .before = choiceset
  )


choice_vars <- paste0(
  "CS",
  1:16
)

stopifnot(
  all(
    choice_vars %in%
      names(dce_data)
  )
)

dce_answers <-
  dce_data[
    choice_vars
  ]

other_info <-
  dce_data[
    ,
    !(
      names(dce_data) %in%
        choice_vars
    )
  ]


# Recode 99 values in survey variables as missing values.
other_info <- other_info %>%
  replace_with_na_at(
    .vars = c(
      2:75
    ),
    condition =
      ~ .x == 99
  )

colnames(
  dce_answers
) <- 1:16

dce_answers <-
  tibble::rowid_to_column(
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
  as.numeric(
    choi$choiceset
  )

choi <- choi %>%
  mutate(
    pref1 =
      str_replace(
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
) %>%
  relocate(RID)

dce_data <-
  dce_data[
    order(
      dce_data[,1],
      dce_data[,2]
    ),
  ]

row.names(
  dce_data
) <- NULL


# Create dummy-coded attribute indicators for the conditional logit models.
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


# Sort attribute columns alphabetically after the first four columns.
cols_to_sort <-
  names(
    dce_data
  )[
    5:ncol(
      dce_data
    )
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

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)


# Convert Polish cost levels to their euro equivalents.
database <- database %>%
  mutate(
    across(
      ends_with("cost"),
      ~ case_when(
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


# Keep valid DCE choices and exclude the respondent with an incomplete biospheric score.
database <- database %>%
  filter(
    pref1 >= 1 &
      pref1 <= 3,
    RID != 35
  )

rm(
  list =
    ls()[
      ls() != "database"
    ]
)


# -------------------------------------------------------------------------
# Model 1: BAU ASC interactions with LOO_tend and biospheric values
# -------------------------------------------------------------------------

apollo_initialise()

apollo_control <- list(
  modelName =
    "Clogit_biospheric",

  modelDescr =
    "Conditional logit in WTP-space: biospheric values",

  indivID =
    "RID",

  mixing =
    FALSE
)

apollo_beta <- c(
  b0 = -77.3839,
  b0_LOO = 0,
  b0_BV = 0,
  b10 = 56.3098,
  b11 = 12.4279,
  b2 = 0.2863,
  b30 = 32.9715,
  b31 = 38.0468,
  b4 = -0.0150
)

apollo_fixed <- c()

apollo_inputs <-
  apollo_validateInputs()


apollo_probabilities <- function(
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

  P <- list()
  V <- list()

  V[["alt1"]] <-
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

  V[["alt2"]] <-
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

  V[["alt3"]] <-
    -b4 * (
      b0 +
        b0_LOO *
        LOO_tend +
        b0_BV *
        biospheric_c +
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

  mnl_settings <- list(
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

  P[["model"]] <-
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P <-
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P <-
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}


model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings =
    list(
      hessianRoutine =
        "maxLik"
    )
)

model1 <- model

model1_output <-
  as.data.frame(
    apollo_modelOutput(
      model1,
      modelOutput_settings =
        list(
          printPVal = 1
        )
    )
  )

# Add one-sided p-values based on the robust standard errors used in the dissertation.
model1_output[["Rob.p(1-sided)"]] <-
  pnorm(
    -abs(
      model1_output[["Estimate"]] /
        model1_output[["Rob.s.e."]]
    )
  )

saveRDS(
  model1,
  "./outputs/paper_C/models/02_clogit_biospheric.rds"
)

model1_output$Parameter <-
  rownames(
    model1_output
  )

rownames(
  model1_output
) <- NULL

model1_output <-
  model1_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(model1_output),
        "Parameter"
      )
    )
  ]

write.csv(
  model1_output,
  "./outputs/paper_C/tables/02_clogit_biospheric.csv",
  row.names = FALSE
)

rm(
  list =
    ls()[
      !ls() %in%
        c(
          "model1",
          "database"
        )
    ]
)


# -------------------------------------------------------------------------
# Model 2: BAU ASC interactions with LOO_tend and altruistic values
# -------------------------------------------------------------------------

apollo_initialise()

apollo_control <- list(
  modelName =
    "Clogit_altruistic",

  modelDescr =
    "Conditional logit in WTP-space: altruistic values",

  indivID =
    "RID",

  mixing =
    FALSE
)

apollo_beta <- c(
  b0 = -77.3839,
  b0_LOO = 0,
  b0_AV = 0,
  b10 = 56.3098,
  b11 = 12.4279,
  b2 = 0.2863,
  b30 = 32.9715,
  b31 = 38.0468,
  b4 = -0.0150
)

apollo_fixed <- c()

apollo_inputs <-
  apollo_validateInputs()


apollo_probabilities <- function(
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

  P <- list()
  V <- list()

  V[["alt1"]] <-
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

  V[["alt2"]] <-
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

  V[["alt3"]] <-
    -b4 * (
      b0 +
        b0_LOO *
        LOO_tend +
        b0_AV *
        altruistic_c +
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

  mnl_settings <- list(
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

  P[["model"]] <-
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P <-
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P <-
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}


model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings =
    list(
      hessianRoutine =
        "maxLik"
    )
)

model2 <- model

model2_output <-
  as.data.frame(
    apollo_modelOutput(
      model2,
      modelOutput_settings =
        list(
          printPVal = 1
        )
    )
  )

# Add one-sided p-values based on the robust standard errors used in the dissertation.
model2_output[["Rob.p(1-sided)"]] <-
  pnorm(
    -abs(
      model2_output[["Estimate"]] /
        model2_output[["Rob.s.e."]]
    )
  )

saveRDS(
  model2,
  "./outputs/paper_C/models/02_clogit_altruistic.rds"
)

model2_output$Parameter <-
  rownames(
    model2_output
  )

rownames(
  model2_output
) <- NULL

model2_output <-
  model2_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(model2_output),
        "Parameter"
      )
    )
  ]

write.csv(
  model2_output,
  "./outputs/paper_C/tables/02_clogit_altruistic.csv",
  row.names = FALSE
)

rm(
  list =
    ls()[
      !ls() %in%
        c(
          "model1",
          "model2",
          "database"
        )
    ]
)


# -------------------------------------------------------------------------
# Model 3: BAU ASC interactions with LOO_tend and hedonic values
# -------------------------------------------------------------------------

apollo_initialise()

apollo_control <- list(
  modelName =
    "Clogit_hedonic",

  modelDescr =
    "Conditional logit in WTP-space: hedonic values",

  indivID =
    "RID",

  mixing =
    FALSE
)

apollo_beta <- c(
  b0 = -77.3839,
  b0_LOO = 0,
  b0_HV = 0,
  b10 = 56.3098,
  b11 = 12.4279,
  b2 = 0.2863,
  b30 = 32.9715,
  b31 = 38.0468,
  b4 = -0.0150
)

apollo_fixed <- c()

apollo_inputs <-
  apollo_validateInputs()


apollo_probabilities <- function(
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

  P <- list()
  V <- list()

  V[["alt1"]] <-
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

  V[["alt2"]] <-
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

  V[["alt3"]] <-
    -b4 * (
      b0 +
        b0_LOO *
        LOO_tend +
        b0_HV *
        hedonic_c +
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

  mnl_settings <- list(
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

  P[["model"]] <-
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P <-
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P <-
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}


model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings =
    list(
      hessianRoutine =
        "maxLik"
    )
)

model3 <- model

model3_output <-
  as.data.frame(
    apollo_modelOutput(
      model3,
      modelOutput_settings =
        list(
          printPVal = 1
        )
    )
  )

# Add one-sided p-values based on the robust standard errors used in the dissertation.
model3_output[["Rob.p(1-sided)"]] <-
  pnorm(
    -abs(
      model3_output[["Estimate"]] /
        model3_output[["Rob.s.e."]]
    )
  )

saveRDS(
  model3,
  "./outputs/paper_C/models/02_clogit_hedonic.rds"
)

model3_output$Parameter <-
  rownames(
    model3_output
  )

rownames(
  model3_output
) <- NULL

model3_output <-
  model3_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(model3_output),
        "Parameter"
      )
    )
  ]

write.csv(
  model3_output,
  "./outputs/paper_C/tables/02_clogit_hedonic.csv",
  row.names = FALSE
)

rm(
  list =
    ls()[
      !ls() %in%
        c(
          "model1",
          "model2",
          "model3",
          "database"
        )
    ]
)


# -------------------------------------------------------------------------
# Model 4: BAU ASC interactions with LOO_tend and egoistic values
# -------------------------------------------------------------------------

apollo_initialise()

apollo_control <- list(
  modelName =
    "Clogit_egoistic",

  modelDescr =
    "Conditional logit in WTP-space: egoistic values",

  indivID =
    "RID",

  mixing =
    FALSE
)

apollo_beta <- c(
  b0 = -77.3839,
  b0_LOO = 0,
  b0_EV = 0,
  b10 = 56.3098,
  b11 = 12.4279,
  b2 = 0.2863,
  b30 = 32.9715,
  b31 = 38.0468,
  b4 = -0.0150
)

apollo_fixed <- c()

apollo_inputs <-
  apollo_validateInputs()


apollo_probabilities <- function(
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

  P <- list()
  V <- list()

  V[["alt1"]] <-
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

  V[["alt2"]] <-
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

  V[["alt3"]] <-
    -b4 * (
      b0 +
        b0_LOO *
        LOO_tend +
        b0_EV *
        egoistic_c +
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

  mnl_settings <- list(
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

  P[["model"]] <-
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P <-
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P <-
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}


model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings =
    list(
      hessianRoutine =
        "maxLik"
    )
)

model4 <- model

model4_output <-
  as.data.frame(
    apollo_modelOutput(
      model4,
      modelOutput_settings =
        list(
          printPVal = 1
        )
    )
  )

# Add one-sided p-values based on the robust standard errors used in the dissertation.
model4_output[["Rob.p(1-sided)"]] <-
  pnorm(
    -abs(
      model4_output[["Estimate"]] /
        model4_output[["Rob.s.e."]]
    )
  )

saveRDS(
  model4,
  "./outputs/paper_C/models/02_clogit_egoistic.rds"
)

model4_output$Parameter <-
  rownames(
    model4_output
  )

rownames(
  model4_output
) <- NULL

model4_output <-
  model4_output[
    ,
    c(
      "Parameter",
      setdiff(
        names(model4_output),
        "Parameter"
      )
    )
  ]

write.csv(
  model4_output,
  "./outputs/paper_C/tables/02_clogit_egoistic.csv",
  row.names = FALSE
)

rm(
  list =
    ls()[
      !ls() %in%
        c(
          "model1",
          "model2",
          "model3",
          "model4",
          "database"
        )
    ]
)


############################
# 3. Ordered logistic regression results –
#    estimating perceived influence of other workshop participants
############################

# Required files in the data/ directory:
# 1) DMV data_EN.csv
#    Source: https://zenodo.org/records/17710200
# 2) Pro_contra_statements_respondents_csv_UTF_8.csv
#    Source: supplementary materials

rm(list = ls())

library(dplyr)
library(tidyr)
library(MASS)
library(pscl)


# Read main survey and DCE data.
dce_data <- read.csv2(
  "./data/DMV data_EN.csv"
) %>%
  rename(
    RID = ID
  )


# Calculate E-PVQ value scores and grand-mean-centred value scores.
dce_data <- dce_data %>%
  mutate(
    egoistic = round(
      rowMeans(
        across(
          c(
            Q4.1,
            Q4.3,
            Q4.6,
            Q4.10,
            Q4.15
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    egoistic_c =
      egoistic -
      mean(
        egoistic,
        na.rm = TRUE
      ),

    hedonic = round(
      rowMeans(
        across(
          c(
            Q4.2,
            Q4.5,
            Q4.11
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    hedonic_c =
      hedonic -
      mean(
        hedonic,
        na.rm = TRUE
      ),

    altruistic = round(
      rowMeans(
        across(
          c(
            Q4.4,
            Q4.7,
            Q4.8,
            Q4.13,
            Q4.17
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    altruistic_c =
      altruistic -
      mean(
        altruistic,
        na.rm = TRUE
      ),

    biospheric = round(
      rowMeans(
        across(
          c(
            Q4.9,
            Q4.12,
            Q4.14,
            Q4.16
          )
        ),
        na.rm = TRUE
      ),
      2
    ),
    biospheric_c =
      biospheric -
      mean(
        biospheric,
        na.rm = TRUE
      )
  )


# Calculate leave-one-out divergence from the workshop mean for biospheric values.
dce_data <- dce_data %>%
  group_by(WS_ID) %>%
  mutate(
    biospheric_WS_sum =
      sum(
        biospheric,
        na.rm = TRUE
      ),

    biospheric_WS_nobs =
      sum(
        !is.na(
          biospheric
        )
      ),

    biospheric_group_mean_LOO =
      case_when(
        is.na(biospheric) ~ NA_real_,

        biospheric_WS_nobs < 2 ~ NA_real_,

        TRUE ~
          (
            biospheric_WS_sum -
              biospheric
          ) /
          (
            biospheric_WS_nobs -
              1
          )
      ),

    biospheric_in_group_diff_LOO =
      abs(
        biospheric -
          biospheric_group_mean_LOO
      )
  ) %>%
  ungroup() %>%
  dplyr::select(
    -biospheric_WS_sum,
    -biospheric_WS_nobs
  )


# Read respondent-level coded pro/contra statements and aggregate them by respondent.
ind_stat <- read.csv2(
  "./data/Pro_contra_statements_respondents_csv_UTF_8.csv",
  fileEncoding = "UTF-8-BOM"
) %>%
  filter(!is.na(RID))

ind_stat_merged <- ind_stat %>%
  group_by(RID) %>%
  summarise(
    Pro = sum(
      Pro,
      na.rm = TRUE
    ),
    Contra = sum(
      Contra,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Join statement totals and calculate leave-one-out workshop-level relative disagreement.
dce_data <- dce_data %>%
  left_join(
    ind_stat_merged,
    by = "RID"
  ) %>%
  mutate(
    Pro = replace_na(
      Pro,
      0
    ),
    Contra = replace_na(
      Contra,
      0
    )
  ) %>%
  group_by(WS_ID) %>%
  mutate(
    WS_policy_pro =
      sum(
        Pro,
        na.rm = TRUE
      ),

    WS_policy_contra =
      sum(
        Contra,
        na.rm = TRUE
      ),

    WS_policy_pro_LOO =
      WS_policy_pro -
      Pro,

    WS_policy_contra_LOO =
      WS_policy_contra -
      Contra,

    WS_rel_diff_LOO =
      if_else(
        WS_policy_pro_LOO +
          WS_policy_contra_LOO > 0,

        abs(
          WS_policy_pro_LOO -
            WS_policy_contra_LOO
        ) /
          (
            WS_policy_pro_LOO +
              WS_policy_contra_LOO
          ),

        NA_real_
      )
  ) %>%
  ungroup()


# Prepare ordered response variable for ordered logistic regression.
dce_data <- dce_data %>%
  mutate(
    Q17.2 =
      factor(
        Q17.2,
        ordered = TRUE
      )
  )


# Fit ordered logistic regression models.
model1 <- MASS::polr(
  Q17.2 ~
    biospheric_c +
    biospheric_in_group_diff_LOO +
    WS_rel_diff_LOO,
  data = dce_data,
  Hess = TRUE
)

model2 <- MASS::polr(
  Q17.2 ~
    biospheric_c *
    biospheric_in_group_diff_LOO +
    WS_rel_diff_LOO,
  data = dce_data,
  Hess = TRUE
)


# Extract coefficient estimates, standard errors, z-values, and Wald p-values.
polr_results <- function(model) {

  coefs <-
    coef(model)

  vcov_matrix <-
    vcov(model)

  se <-
    sqrt(
      diag(
        vcov_matrix
      )
    )[
      names(coefs)
    ]

  z_values <-
    coefs /
    se

  p_values <-
    2 *
    (
      1 -
        pnorm(
          abs(
            z_values
          )
        )
    )

  data.frame(
    Parameter =
      names(coefs),

    Estimate =
      unname(coefs),

    Std_Error =
      unname(se),

    z_value =
      unname(z_values),

    p_value =
      unname(p_values),

    row.names = NULL
  )
}


# Display model results and pseudo-R2 statistics.
summary(model1)

model1_results <-
  polr_results(
    model1
  )

model1_results

pR2_model1 <-
  pscl::pR2(
    model1
  )

pR2_model1


summary(model2)

model2_results <-
  polr_results(
    model2
  )

model2_results

pR2_model2 <-
  pscl::pR2(
    model2
  )

pR2_model2


# Save complete ordered logistic regression models.
saveRDS(
  list(
    model1 = model1,
    model2 = model2
  ),
  "./outputs/paper_C/models/03_perceived_influence_ordered_logit_models.rds"
)


# Combine and save ordered logistic regression coefficient tables.
model1_results$Model <-
  "Model 1"

model2_results$Model <-
  "Model 2"

ordered_logit_results <-
  bind_rows(
    model1_results,
    model2_results
  )

ordered_logit_results <-
  ordered_logit_results[
    ,
    c(
      "Model",
      "Parameter",
      "Estimate",
      "Std_Error",
      "z_value",
      "p_value"
    )
  ]

write.csv(
  ordered_logit_results,
  "./outputs/paper_C/tables/03_perceived_influence_coefficients.csv",
  row.names = FALSE
)


# Save pseudo-R2 statistics.
pseudo_r2 <-
  as.data.frame(
    rbind(
      model1 =
        pR2_model1,

      model2 =
        pR2_model2
    )
  )

pseudo_r2$Model <-
  rownames(
    pseudo_r2
  )

rownames(
  pseudo_r2
) <- NULL

pseudo_r2 <-
  pseudo_r2[
    ,
    c(
      "Model",
      setdiff(
        names(pseudo_r2),
        "Model"
      )
    )
  ]

write.csv(
  pseudo_r2,
  "./outputs/paper_C/tables/03_perceived_influence_pseudo_R2.csv",
  row.names = FALSE
)


############################
# Optional starting-value search
############################

# This search documents how the starting values used in Section 2 were obtained.
# It is computationally intensive and is therefore skipped by default.
# Set to TRUE to reproduce the starting-value search.
run_starting_value_search <- FALSE

if (run_starting_value_search) {

############################
# 4. Documentation of starting-value search for conditional logit models in WTP-space
############################

# Required files in the data/ directory:
# 1) final design_baesyian efficient design with interaction.NGD
# 2) DMV data_EN.csv
#    Source: https://zenodo.org/records/17710200

rm(list = ls())

library(dplyr)
library(naniar)
library(tidyverse)
library(apollo)


# Read data used for starting-value search.
dce_data <- read.csv2(
  "./data/DMV data_EN.csv"
)


# Read experimental design.
design <- read_delim(
  "./data/final design_baesyian efficient design with interaction.NGD",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  col_select = c(
    -Design,
    -starts_with("...")
  ),
  name_repair = "universal",
  n_max = 16
) %>%
  filter(
    !is.na(
      Choice.situation
    )
  )

nsets <- nrow(
  design
)

design <- design %>%
  rename(
    choiceset =
      Choice.situation
  )

design$alt3.plankton <- 3
design$alt3.carbon <- 0
design$alt3.mpa <- 3
design$alt3.cost <- 0

design <- design %>%
  relocate(
    Block,
    .before = choiceset
  )


# Separate DCE choices from respondent-level variables.
dce_answers <-
  dce_data[
    ,
    44:59
  ]

other_info <-
  dce_data[
    ,
    -(44:59)
  ] %>%
  rename(
    RID = ID
  )

other_info <- other_info %>%
  replace_with_na_at(
    .vars =
      c(
        2:75
      ),
    condition =
      ~ .x == 99
  )

colnames(
  dce_answers
) <- 1:16

dce_answers <-
  tibble::rowid_to_column(
    dce_answers,
    "RID"
  )


# Reshape DCE responses to long format.
choi <- dce_answers %>%
  pivot_longer(
    cols = 2:17,
    names_to = "choiceset",
    values_to = "pref1"
  )

choi$choiceset <-
  as.numeric(
    choi$choiceset
  )

choi <- choi %>%
  mutate(
    pref1 =
      str_replace(
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
) %>%
  relocate(RID)

dce_data <-
  dce_data[
    order(
      dce_data[,1],
      dce_data[,2]
    ),
  ]

row.names(
  dce_data
) <- NULL


# Create dummy-coded attribute indicators for the conditional logit model.
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


# Sort attribute columns alphabetically after the first four columns.
cols_to_sort <-
  names(
    dce_data
  )[
    5:ncol(
      dce_data
    )
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

database <- merge(
  dce_data,
  other_info,
  by = "RID",
  all = TRUE
)


# Convert Polish cost levels to their euro equivalents.
database <- database %>%
  mutate(
    across(
      ends_with("cost"),
      ~ case_when(
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


# Keep valid DCE choices and exclude protest respondents used in the starting-value search.
database <- database %>%
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


# Estimate a baseline conditional logit model for starting-value search.
apollo_initialise()

apollo_control <- list(
  modelName =
    "Clogit_PLN",

  modelDescr =
    "Conditional logit in WTP-space used for starting-value search",

  indivID =
    "RID",

  mixing =
    FALSE
)

apollo_beta <- c(
  b0 = 0,
  b10 = 0,
  b11 = 0,
  b2 = 0,
  b30 = 0,
  b31 = 0,
  b4 = 0
)

apollo_fixed <- c()

apollo_inputs <-
  apollo_validateInputs()


apollo_probabilities <- function(
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

  P <- list()
  V <- list()

  V[["alt1"]] <-
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

  V[["alt2"]] <-
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

  V[["alt3"]] <-
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

  mnl_settings <- list(
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

  P[["model"]] <-
    apollo_mnl(
      mnl_settings,
      functionality
    )

  P <-
    apollo_panelProd(
      P,
      apollo_inputs,
      functionality
    )

  P <-
    apollo_prepareProb(
      P,
      apollo_inputs,
      functionality
    )

  return(P)
}


model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)


searchStart_settings <- list(
  searchRoutine = "gss",
  nSearches = 10,
  nSearchDraws = 1000,
  searchTol = 1e-6,
  searchMaxEval = 1000,
  searchMaxIter = 100,
  searchParallel = FALSE
)


apollo_beta <- apollo_searchStart(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  searchStart_settings
)


# Save resulting starting-value candidate.
starting_values <- data.frame(
  Parameter =
    names(
      apollo_beta
    ),

  Starting_value =
    unname(
      apollo_beta
    )
)

write.csv(
  starting_values,
  "./outputs/paper_C/tables/04_conditional_logit_starting_values.csv",
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

} # End optional starting-value search

# Best candidate from starting-value search (LL = -2229.2):
# b0  = -77.3839
# b10 =  56.3098
# b11 =  12.4279
# b2  =   0.2863
# b30 =  32.9715
# b31 =  38.0468
# b4  =  -0.0150
#
# Estimating the conditional logit models with these starting values and with
# all parameters initially set to zero yielded the same final results.
