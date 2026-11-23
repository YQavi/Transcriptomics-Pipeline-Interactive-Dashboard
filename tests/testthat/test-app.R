# Smoke test: the dashboard starts and every output renders without a Shiny error.
# Skipped automatically unless shinytest2 is installed and a data object exists.
library(testthat)

test_that("dashboard loads and renders without output errors", {
  skip_if_not_installed("shinytest2")
  root <- normalizePath(file.path(getwd(), "..", ".."))
  skip_if_not(file.exists(file.path(root, "app", "data", "demo_deseq2_data.RData")) ||
              file.exists(file.path(root, "results", "deseq2", "app_data.RData")),
              "no dashboard data object yet (run step 5)")

  app <- shinytest2::AppDriver$new(app_dir = file.path(root, "app"), name = "smoke", load_timeout = 60000)
  on.exit(app$stop(), add = TRUE)
  app$wait_for_idle()

  n_err <- function() app$get_js("document.querySelectorAll('.shiny-output-error').length")
  expect_equal(n_err(), 0)

  app$set_inputs(pval_thresh = 0.01, lfc_thresh = 1)    # thresholds drive every tab
  app$set_inputs(sig_filter = "Up-regulated")
  app$wait_for_idle()
  expect_equal(n_err(), 0)
})
