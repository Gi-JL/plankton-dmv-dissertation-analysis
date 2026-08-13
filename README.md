# Plankton DMV Dissertation Analysis

This repository contains the R code used for the quantitative analyses underlying a cumulative dissertation on the deliberative monetary valuation of marine biodiversity and plankton ecosystem services.

The analyses are based on data collected during 15 deliberative monetary valuation (DMV) workshops conducted in five European study regions between October 2023 and February 2024.

## Repository structure

```text
.
├── README.md
├── run_all.R
├── R/
│   ├── 01_paper_A.R
│   ├── 02_paper_B.R
│   └── 03_paper_C.R
└── data/
    ├── README.md
    └── Pro_contra_statements_respondents_csv_UTF_8.csv
```

The three analysis scripts correspond to the three papers included in the cumulative dissertation:

- `R/01_paper_A.R` – analyses for Paper A
- `R/02_paper_B.R` – analyses for Paper B
- `R/03_paper_C.R` – analyses for Paper C
- `run_all.R` – executes all three paper-specific analysis scripts

## Data

Most input data are not stored directly in this repository. Download instructions and information on the required files are provided in `data/README.md`.

The coded discussion data required for Paper C are included directly in the `data/` directory.

## Software requirements

The analyses are conducted in R and require the following packages:

- `apollo`
- `broom`
- `dplyr`
- `MASS`
- `naniar`
- `pscl`
- `readr`
- `sf`
- `stargazer`
- `stringr`
- `tibble`
- `tidyr`
- `tidyverse`

The following R code checks which required packages are already installed and installs any missing packages:

```r
required_packages <- c(
  "apollo",
  "broom",
  "dplyr",
  "MASS",
  "naniar",
  "pscl",
  "readr",
  "sf",
  "stargazer",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}
```

## Running the analyses

After obtaining the required input data and placing them in the `data/` directory as described in `data/README.md`, individual analyses can be run from the repository root using:

```bash
Rscript R/01_paper_A.R
Rscript R/02_paper_B.R
Rscript R/03_paper_C.R
```

Alternatively, all three scripts can be executed sequentially using:

```bash
Rscript run_all.R
```
