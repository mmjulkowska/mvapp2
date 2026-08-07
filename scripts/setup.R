# One-time local setup for MVApp 2.0 development.
# Run from the repo root:  Rscript scripts/setup.R

# 1. Reproducible environment (recommended). Uncomment to adopt renv:
# install.packages("renv"); renv::init()

# 2. Packages needed to run the app + backend in dev mode.
pkgs <- c("shiny", "DT", "ggplot2", "corrplot", "plotly", "RColorBrewer",
          "ggsci", "pheatmap", "dendextend", "multcompView", "cluster",
          "devtools", "testthat")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

# 3. Load + test the backend package.
devtools::load_all("mvapp")
testthat::test_dir("mvapp/tests/testthat")

cat("\nSetup complete. Run the app with:\n",
    "  setwd('app'); shiny::runApp()\n", sep = "")
