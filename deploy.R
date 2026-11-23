# Deploy the dashboard to shinyapps.io. Run from the repository root after the pipeline has
# created app/data/demo_deseq2_data.RData.
#
# One-time account setup (token from your shinyapps.io dashboard):
#   rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
stopifnot(file.exists("app/data/demo_deseq2_data.RData"))
rsconnect::deployApp(
  appDir   = "app",
  appName  = "cervical-caf-rnaseq-dashboard",
  appFiles = c("app.R", "data/demo_deseq2_data.RData")
)
