# MVApp 2.0

Modernized, package-backed rebuild of **[MVApp](https://github.com/mmjulkowska/MVApp)** —
the multivariate analysis application for plant-phenotyping data
([Julkowska et al., 2019, *Plant Physiology*](https://doi.org/10.1104/pp.19.00509)).

> **This repo does not touch the production app.** The live instance on the KAUST
> server keeps running whatever is deployed there. MVApp 2.0 is developed and tested
> here (and locally) and only reaches production if/when you deliberately deploy it.

## What's different from MVApp 1.x

The science moves out of a single 15k-line `server.R` into a **standalone, testable R
package (`mvapp`)**, with Shiny as a thin UI on top. Every button in the app calls a
named function in `mvapp::` — no statistics live in a server block. That makes the code
reviewable, unit-tested, and reusable in lab RMarkdown notebooks.

## Layout

```
MVApp2/
├── mvapp/                # the analysis backend (an R package)
│   ├── R/                # pure, tested functions
│   │   ├── data-contract.R   # mvapp_dataset: the one structure every analysis consumes
│   │   ├── io-generic.R      # read_long_csv()
│   │   ├── io-phenosight.R    # PS2 / RGB / watering readers        [Roadmap Phase 2]
│   │   ├── toe.R             # Time-Of-Experiment
│   │   ├── tolerance.R       # STI / STI1 / STI2
│   │   ├── ml-importance.R    # RF / LASSO / GBM importance          [Roadmap Phase 4]
│   │   └── path-analysis.R    # lavaan SEM / path analysis           [Roadmap Phase 4]
│   ├── tests/testthat/   # unit tests (run on every push via CI)
│   └── inst/extdata/     # small example dataset
├── app/                  # the Shiny app (thin UI)
│   ├── app.R
│   └── R/                # one module per analysis tab (mod_upload, mod_correlations, …)
├── docker/               # local + deployable image
├── .github/workflows/    # R CMD check + tests
└── scripts/setup.R       # one-time local setup
```

## Run it locally

**Option A — RStudio / R (fastest for development):**

```r
Rscript scripts/setup.R      # installs deps, loads + tests the backend
setwd("app"); shiny::runApp() # opens the app; loads mvapp from ../mvapp in dev mode
```

**Option B — Docker (mirrors production):**

```bash
docker compose -f docker/docker-compose.yml up --build
# then open http://localhost:3838
```

## Run the tests

```r
devtools::test("mvapp")     # or: Rscript -e 'testthat::test_dir("mvapp/tests/testthat")'
```

## The development loop (from the roadmap)

For each of the 12 legacy analysis tabs, one at a time:

1. **Extract** the statistics into a pure function in `mvapp/R/` (input = data frame + params, output = data frame / model / plot).
2. **Test** it in `mvapp/tests/testthat/`.
3. **Wrap** the tab as a Shiny module in `app/R/` that only does I/O + rendering.
4. **Delete** the old block; confirm CI is green; commit.

`mod_correlations` is a worked example of steps 1–3. Grow the app by repeating the loop
and by filling in the `[Roadmap Phase N]` stubs.

## Roadmap

The full phased plan (consolidation → PhenoSight import → PhenoSight analysis → ML /
path analysis) lives in the project roadmap document.

## License

Apache License 2.0 (carried over from MVApp 1.x).
