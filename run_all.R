############################
# Run all dissertation analysis scripts
############################

# Run this file from the repository root.
# Each paper-specific script is executed in a separate R process so that
# workspace-clearing commands within one script do not affect the others.

scripts <- c(
  "R/01_paper_A.R",
  "R/02_paper_B.R",
  "R/03_paper_C.R"
)

missing_scripts <- scripts[!file.exists(scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "The following analysis script(s) could not be found:\n",
    paste(missing_scripts, collapse = "\n")
  )
}

rscript <- file.path(R.home("bin"), "Rscript")

if (.Platform$OS.type == "windows") {
  rscript <- paste0(rscript, ".exe")
}

for (script in scripts) {

  message(
    "\n========================================\n",
    "Running: ", script, "\n",
    "========================================"
  )

  status <- system2(
    command = rscript,
    args = shQuote(script)
  )

  if (status != 0) {
    stop("Execution failed for: ", script)
  }
}

message("\nAll analysis scripts completed successfully.")
