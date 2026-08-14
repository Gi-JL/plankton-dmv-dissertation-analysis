############################
# Analysis code for Paper B of the cumulative dissertation
# Wilcoxon and Friedman tests and ordered logistic regression models
# Input-data sources and download instructions: see data/README.md
############################

# Run this script from the repository root:
# source("R/02_paper_B.R")
#
# Alternatively, from a terminal:
# Rscript R/02_paper_B.R
#
# Relative paths such as "./data/..." assume that the repository root
# is the current working directory.

############## Overview ###############
# 1. Data import and preparation
# 2. Wilcoxon rank-sum tests (Coastal vs Inland)
# 3. Friedman tests and pairwise signed-rank comparisons
# 4. Ordered logistic regression results – estimating TEV category importance using personal value scores
# 5. Ordered logistic regression results – estimating ES category importance using personal value scores and socio-demographic data
# 6. Ordered logistic regression results – estimating TEV category importance using personal value scores and socio-demographic data
# 7. Ordered logistic regression results – estimating DCE attribute importance using personal value scores and socio-demographic data
# 8. Ordered logistic regression results – estimating ES category importance using personal value scores and socio-demographic data, using coastal/inland group membership instead of distance from the coast
#######################################


############################
# Output directories
############################

dir.create(
  "./outputs/paper_B/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "./outputs/paper_B/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


######################
# 1. Data import and preparation
######################

rm(list=ls())

library(MASS)
library(dplyr)
library(naniar)
library(tidyr)
library(broom)
library(pscl)
library(sf)
library(stargazer)


############################
# Helper function for saving ordered logistic regression results
############################

# Saves:
# 1) the complete set of fitted polr model objects as an RDS file;
# 2) a structured CSV containing predictor coefficients, standard errors,
#    z-values and two-sided Wald p-values;
# 3) a CSV containing basic model-fit information;
# 4) the original stargazer text output.

save_polr_outputs <- function(models,
                              outcomes,
                              prefix,
                              title) {

  names(models) <- outcomes

  # Save complete fitted model objects
  saveRDS(
    models,
    file.path(
      "./outputs/paper_B/models",
      paste0(prefix, "_models.rds")
    )
  )

  # Extract predictor coefficients from each model
  coefficient_tables <- lapply(
    seq_along(models),
    function(i) {

      model <- models[[i]]

      coef_matrix <- coef(summary(model))

      predictor_names <- names(coef(model))

      coef_matrix <- coef_matrix[
        predictor_names,
        ,
        drop = FALSE
      ]

      z_values <- coef_matrix[,3]

      data.frame(
        Outcome = outcomes[i],
        Parameter = rownames(coef_matrix),
        Estimate = coef_matrix[,1],
        Std_Error = coef_matrix[,2],
        z_value = z_values,
        p_value = 2 * pnorm(
          abs(z_values),
          lower.tail = FALSE
        ),
        row.names = NULL
      )
    }
  )

  coefficient_table <- bind_rows(
    coefficient_tables
  )

  write.csv(
    coefficient_table,
    file.path(
      "./outputs/paper_B/tables",
      paste0(prefix, "_coefficients.csv")
    ),
    row.names = FALSE
  )

  # Save basic model-fit information
  fit_table <- data.frame(
    Outcome = outcomes,
    Observations = sapply(
      models,
      nobs
    ),
    Residual_deviance = sapply(
      models,
      deviance
    ),
    AIC = sapply(
      models,
      AIC
    ),
    row.names = NULL
  )

  write.csv(
    fit_table,
    file.path(
      "./outputs/paper_B/tables",
      paste0(prefix, "_model_fit.csv")
    ),
    row.names = FALSE
  )

  # Preserve the original stargazer representation
  stargazer_output <- capture.output(
    stargazer(
      models,
      type = "text",
      title = title,
      dep.var.labels = outcomes
    )
  )

  # Print to console as in the original script
  cat(
    paste(
      stargazer_output,
      collapse = "\n"
    ),
    "\n"
  )

  # Save text output
  writeLines(
    stargazer_output,
    file.path(
      "./outputs/paper_B/tables",
      paste0(prefix, "_stargazer.txt")
    )
  )

  invisible(
    list(
      models = models,
      coefficients = coefficient_table,
      fit = fit_table
    )
  )
}


all_data <- read.csv2(
  "./data/DMV data.csv"
)

all_data <- all_data %>%
  rename(
    RID = ID
  )

all_data <- all_data %>%
  replace_with_na_at(
    .vars = c(2:91),
    condition = ~.x == 99
  ) %>%
  filter(
    !is.na(Q1)
  )

all_data$Age <-
  2023 -
  all_data$`Year.of.birth`

all_data <- all_data %>%
  mutate(
    egoistic = round(
      rowMeans(
        across(
          c(
            "Q4.1",
            "Q4.3",
            "Q4.6",
            "Q4.10",
            "Q4.15"
          )
        ),
        na.rm = TRUE
      ),
      2
    ),

    hedonic = round(
      rowMeans(
        across(
          c(
            "Q4.2",
            "Q4.5",
            "Q4.11"
          )
        ),
        na.rm = TRUE
      ),
      2
    ),

    altruistic = round(
      rowMeans(
        across(
          c(
            "Q4.4",
            "Q4.7",
            "Q4.8",
            "Q4.13",
            "Q4.17"
          )
        ),
        na.rm = TRUE
      ),
      2
    ),

    biospheric = round(
      rowMeans(
        across(
          c(
            "Q4.9",
            "Q4.12",
            "Q4.14",
            "Q4.16"
          )
        ),
        na.rm = TRUE
      ),
      2
    )
  ) %>%

  mutate(
    egoistic_c =
      egoistic -
      mean(
        egoistic,
        na.rm = TRUE
      ),

    hedonic_c =
      hedonic -
      mean(
        hedonic,
        na.rm = TRUE
      ),

    altruistic_c =
      altruistic -
      mean(
        altruistic,
        na.rm = TRUE
      ),

    biospheric_c =
      biospheric -
      mean(
        biospheric,
        na.rm = TRUE
      )
  )

all_data$Coastal <-
  as.factor(
    all_data$Coastal
  )

all_data$Postcode[
  all_data$WS_ID == 7
] <-
  paste0(
    "0",
    all_data$Postcode[
      all_data$WS_ID == 7
    ]
  )

eu_pc <- st_read(
  "./data/PCODE_PT_2020_3035.shp"
)

all_data$Postcode[
  all_data$Country == "Poland"
] <-
  paste0(
    "PL_",
    all_data$Postcode[
      all_data$Country == "Poland"
    ]
  )

all_data$Postcode[
  all_data$Country == "Poland"
] <-
  sub(
    "^(.{5})",
    "\\1-",
    all_data$Postcode[
      all_data$Country == "Poland"
    ]
  )

all_data$Postcode[
  all_data$Country == "Italy"
] <-
  paste0(
    "IT_",
    all_data$Postcode[
      all_data$Country == "Italy"
    ]
  )

all_data$Postcode[
  all_data$Country == "Basque"
] <-
  paste0(
    "ES_",
    all_data$Postcode[
      all_data$Country == "Basque"
    ]
  )

all_data$Postcode[
  all_data$Country == "Germany"
] <-
  paste0(
    "DE_",
    all_data$Postcode[
      all_data$Country == "Germany"
    ]
  )

all_data$Postcode[
  all_data$Country == "France"
] <-
  paste0(
    "FR_",
    all_data$Postcode[
      all_data$Country == "France"
    ]
  )

all_data <- all_data %>%
  rename(
    PC_CNTR = Postcode
  )

coords <- eu_pc %>%
  dplyr::select(
    PC_CNTR
  ) %>%
  as.data.frame()

all_data <- merge(
  all_data,
  coords,
  by = "PC_CNTR",
  all.x = TRUE
)

all_data <- all_data %>%
  relocate(
    PC_CNTR,
    .after = ISCED
  )

all_data_points <-
  sf::st_as_sf(
    all_data
  )

coastline <- st_read(
  "./data/Europe_coastline/Europe_coastline.shp"
)

st_crs(coastline)
st_crs(all_data_points)

all_data_points <-
  st_transform(
    all_data_points,
    crs = st_crs(coastline)
  )

all_data_points$dist2coast <-
  apply(
    st_distance(
      all_data_points,
      coastline
    ),
    1,
    min
  )

all_data <-
  as.data.frame(
    all_data_points
  )

all_data <- all_data %>%
  rename(
    distkm = dist2coast
  ) %>%
  mutate(
    distkm =
      distkm / 1000
  )

# Some socio-demographic data for respondent 129 are missing.
# This means the respondent cannot be included in all analyses.
all_data_wo129 <-
  all_data %>%
  filter(
    RID != 129
  )


######################
# 2. Wilcoxon rank-sum tests (Coastal vs Inland)
######################

# TEV and non-TEV items
tev_only_vars <- c(
  "Q13.3",
  "Q13.4",
  "Q13.5",
  "Q13.6",
  "Q13.7"
)

non_tev_vars <- c(
  "Q13.8",
  "Q13.9",
  "Q13.10",
  "Q13.11"
)

full_set_vars <- c(
  tev_only_vars,
  non_tev_vars
)

# Keep all TEV and non-TEV items available in the dataset
wilcox_vars_full <-
  full_set_vars[
    full_set_vars %in%
      names(all_data)
  ]

run_wilcox <- function(var) {

  if (
    !all(
      c(
        "Coastal",
        var
      ) %in%
        names(all_data)
    )
  ) {
    return(NULL)
  }

  df <- all_data %>%
    select(
      Coastal,
      all_of(var)
    ) %>%
    drop_na()

  if (
    nrow(df) == 0
  ) {
    return(NULL)
  }

  out <- wilcox.test(
    reformulate(
      "Coastal",
      response = var
    ),
    data = df,
    exact = FALSE
  )

  tibble(
    variable = var,
    p.value = out$p.value,
    W = unname(
      out$statistic
    ),
    n = nrow(df)
  )
}

wilcox_full_results <-
  bind_rows(
    lapply(
      wilcox_vars_full,
      run_wilcox
    )
  )

print(
  wilcox_full_results
)

write.csv(
  wilcox_full_results,
  "./outputs/paper_B/tables/02_wilcoxon_coastal_inland.csv",
  row.names = FALSE
)


######################
# 3. Friedman tests and pairwise signed-rank comparisons
######################

# Select columns corresponding to TEV + non-TEV categories
# TEV (Direct–Option) + non-TEV (fairness + I/We-preference) categories

df <- all_data %>%
  dplyr::select(
    Q13.3,  # Direct use
    Q13.4,  # Indirect use
    Q13.5,  # Bequest
    Q13.6,  # Existence
    Q13.7,  # Option
    Q13.8,  # Intragenerational fairness
    Q13.9,  # Intergenerational fairness
    Q13.10, # I-Preference
    Q13.11  # We-Preference
  )


# Friedman test on TEV-only categories
friedman_tev <-
  friedman.test(
    as.matrix(
      df[,1:5]
    )
  )

print(
  friedman_tev
)


# Friedman test on full set
friedman_full <-
  friedman.test(
    as.matrix(df)
  )

print(
  friedman_full
)


# Save Friedman test results
friedman_results <-
  data.frame(
    Test = c(
      "TEV-only categories",
      "TEV and non-TEV categories"
    ),

    Chi_squared = c(
      unname(
        friedman_tev$statistic
      ),
      unname(
        friedman_full$statistic
      )
    ),

    df = c(
      unname(
        friedman_tev$parameter
      ),
      unname(
        friedman_full$parameter
      )
    ),

    p_value = c(
      friedman_tev$p.value,
      friedman_full$p.value
    )
  )

write.csv(
  friedman_results,
  "./outputs/paper_B/tables/03_friedman_tests.csv",
  row.names = FALSE
)


# Pairwise Wilcoxon signed-rank tests
combinations <-
  combn(
    names(df),
    2,
    simplify = FALSE
  )

p_values <-
  sapply(
    combinations,
    function(cols) {

      test <-
        wilcox.test(
          df[[cols[1]]],
          df[[cols[2]]],
          paired = TRUE
        )

      test$p.value
    }
  )


# Create results table
p_values_table <-
  data.frame(
    comparison =
      sapply(
        combinations,
        function(x) {
          paste(
            x,
            collapse = " vs "
          )
        }
      ),

    p.value =
      p_values
  )


# Add significance markers
p_values_table$significance <-
  ifelse(
    p_values_table$p.value < 0.001,
    "***",

    ifelse(
      p_values_table$p.value < 0.01,
      "**",

      ifelse(
        p_values_table$p.value < 0.05,
        "*",

        ifelse(
          p_values_table$p.value < 0.1,
          ".",
          "ns"
        )
      )
    )
  )

print(
  p_values_table
)

write.csv(
  p_values_table,
  "./outputs/paper_B/tables/03_pairwise_wilcoxon_signed_rank.csv",
  row.names = FALSE
)


######################
# 4. Ordered logistic regression results –
#    estimating TEV category importance using personal value scores
######################

ord_outcomes <- c(
  "Q13.3",
  "Q13.4",
  "Q13.5",
  "Q13.6",
  "Q13.7"
)

for (
  o in ord_outcomes
) {
  all_data[[o]] <-
    factor(
      all_data[[o]],
      ordered = TRUE
    )
}

models_tev_values <-
  lapply(
    ord_outcomes,
    function(o) {

      polr(
        as.formula(
          paste0(
            o,
            " ~ egoistic_c + biospheric_c + altruistic_c + hedonic_c"
          )
        ),
        data = all_data,
        Hess = TRUE
      )
    }
  )

names(models_tev_values) <-
  ord_outcomes

save_polr_outputs(
  models = models_tev_values,
  outcomes = ord_outcomes,
  prefix = "04_tev_personal_values",
  title =
    "Ordered logistic regression results – estimating TEV category importance using personal value scores."
)


######################
# Further data preparation
######################

all_data_wo129 <-
  all_data_wo129 %>%
  mutate(
    age_group =
      ifelse(
        Age <= 24,
        1,

        ifelse(
          Age <= 34,
          2,

          ifelse(
            Age <= 44,
            3,

            ifelse(
              Age <= 54,
              4,

              ifelse(
                Age <= 64,
                5,
                6
              )
            )
          )
        )
      )
  )

all_data_wo129 <-
  all_data_wo129 %>%
  relocate(
    age_group,
    .after = Age
  )

actcol <-
  "Activities.by.the.sea"

all_data_wo129$Swimming_and_bathing <-
  ifelse(
    grepl(
      "^1([^0-9]|$)",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Sun_bathing <-
  ifelse(
    grepl(
      "2",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Surfing_and_sailing <-
  ifelse(
    grepl(
      "3",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Walking <-
  ifelse(
    grepl(
      "4",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Recreational_fishing <-
  ifelse(
    grepl(
      "5",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Diving_and_snorkeling <-
  ifelse(
    grepl(
      "6",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Other_recreational <-
  ifelse(
    grepl(
      "7",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Job_related <-
  ifelse(
    grepl(
      "8",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Other_non_recreational <-
  ifelse(
    grepl(
      "9",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )

all_data_wo129$Living_by_the_sea <-
  ifelse(
    grepl(
      "10",
      all_data_wo129[[actcol]]
    ),
    1,
    0
  )


######################
# 5. Ordered logistic regression results –
#    estimating ES category importance using personal value scores
#    and socio-demographic data
######################

ord_outcomes <- c(
  "Q11.1",
  "Q11.2",
  "Q11.3",
  "Q11.4",
  "Q11.5",
  "Q11.6",
  "Q11.7",
  "Q11.8"
)

for (
  o in ord_outcomes
) {
  all_data_wo129[[o]] <-
    factor(
      all_data_wo129[[o]],
      ordered = TRUE
    )
}

models_es_socio <-
  lapply(
    ord_outcomes,
    function(o) {

      polr(
        as.formula(
          paste0(
            o,
            " ~ egoistic_c + biospheric_c + altruistic_c + hedonic_c + ",
            "`Net.income.group` + Gender + age_group + distkm + ISCED + ",
            "Diving_and_snorkeling + Swimming_and_bathing + Recreational_fishing"
          )
        ),
        data = all_data_wo129,
        Hess = TRUE
      )
    }
  )

names(models_es_socio) <-
  ord_outcomes

save_polr_outputs(
  models = models_es_socio,
  outcomes = ord_outcomes,
  prefix = "05_es_personal_values_sociodemographics",
  title =
    "Ordered logistic regression results – estimating ES category importance using personal value scores and socio-demographic data."
)


######################
# 6. Ordered logistic regression results –
#    estimating TEV category importance using personal value scores
#    and socio-demographic data
######################

ord_outcomes <- c(
  "Q13.3",
  "Q13.4",
  "Q13.5",
  "Q13.6",
  "Q13.7"
)

for (
  o in ord_outcomes
) {
  all_data_wo129[[o]] <-
    factor(
      all_data_wo129[[o]],
      ordered = TRUE
    )
}

models_tev_socio <-
  lapply(
    ord_outcomes,
    function(o) {

      polr(
        as.formula(
          paste0(
            o,
            " ~ egoistic_c + biospheric_c + altruistic_c + hedonic_c + ",
            "`Net.income.group` + Gender + age_group + ISCED + Coastal + ",
            "Diving_and_snorkeling + Swimming_and_bathing + Recreational_fishing"
          )
        ),
        data = all_data_wo129,
        Hess = TRUE
      )
    }
  )

names(models_tev_socio) <-
  ord_outcomes

save_polr_outputs(
  models = models_tev_socio,
  outcomes = ord_outcomes,
  prefix = "06_tev_personal_values_sociodemographics",
  title =
    "Ordered logistic regression results – estimating TEV category importance using personal value scores and socio-demographic data."
)


######################
# 7. Ordered logistic regression results –
#    estimating DCE attribute importance using personal value scores
#    and socio-demographic data
######################

ord_outcomes <- c(
  "Q12.1",
  "Q12.2",
  "Q12.3",
  "Q12.4",
  "Q12.5"
)

for (
  o in ord_outcomes
) {
  all_data_wo129[[o]] <-
    factor(
      all_data_wo129[[o]],
      ordered = TRUE
    )
}

models_dce_socio <-
  lapply(
    ord_outcomes,
    function(o) {

      polr(
        as.formula(
          paste0(
            o,
            " ~ egoistic_c + biospheric_c + altruistic_c + hedonic_c + ",
            "`Net.income.group` + Gender + age_group + ISCED + distkm + ",
            "Diving_and_snorkeling + Swimming_and_bathing + Recreational_fishing"
          )
        ),
        data = all_data_wo129,
        Hess = TRUE
      )
    }
  )

names(models_dce_socio) <-
  ord_outcomes

save_polr_outputs(
  models = models_dce_socio,
  outcomes = ord_outcomes,
  prefix = "07_dce_attributes_personal_values_sociodemographics",
  title =
    "Ordered logistic regression results – estimating DCE attribute importance using personal value scores and socio-demographic data."
)


######################
# 8. Ordered logistic regression results –
#    estimating ES category importance using personal value scores
#    and socio-demographic data, using coastal/inland group membership
#    instead of distance from the coast
######################

ord_outcomes <- c(
  "Q11.1",
  "Q11.2",
  "Q11.3",
  "Q11.4",
  "Q11.5",
  "Q11.6",
  "Q11.7",
  "Q11.8"
)

for (
  o in ord_outcomes
) {
  all_data_wo129[[o]] <-
    factor(
      all_data_wo129[[o]],
      ordered = TRUE
    )
}

models_es_coastal <-
  lapply(
    ord_outcomes,
    function(o) {

      polr(
        as.formula(
          paste0(
            o,
            " ~ egoistic_c + biospheric_c + altruistic_c + hedonic_c + ",
            "`Net.income.group` + Gender + age_group + Coastal + ISCED + ",
            "Diving_and_snorkeling + Swimming_and_bathing + Recreational_fishing"
          )
        ),
        data = all_data_wo129,
        Hess = TRUE
      )
    }
  )

names(models_es_coastal) <-
  ord_outcomes

save_polr_outputs(
  models = models_es_coastal,
  outcomes = ord_outcomes,
  prefix = "08_es_coastal_inland_sociodemographics",
  title =
    "Ordered logistic regression results – estimating ES category importance using personal value scores and socio-demographic data. In this regression coastal or inland group membership was used as independent variable instead of distance of residence from the coast in km."
)
