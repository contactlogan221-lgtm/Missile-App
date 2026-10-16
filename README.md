# Interactive Research Platform for Human Decision-Making Under Uncertainty

Developed by Logan McKinlay as part of research with Prof. Matthew Zefferman, Department of Defense Analysis, Naval Postgraduate School.

**Interactive Application:** https://loganhm.shinyapps.io/missileapp/

---

## About

A browser-based decision-making game built in R Shiny for a behavioural research study. Participants gather probabilistic intelligence about two missile fields, update their beliefs through sampling, and allocate coverage units to neutralise as many missiles as possible. The study examines how different uncertainty visualisations affect the quality of allocation decisions.

## Requirements

The following R packages are required, all available on CRAN:

- shiny
- shinyjs
- httr
- jsonlite

Install them with:

```r
install.packages(c("shiny", "shinyjs", "httr", "jsonlite"))
```

## Running the app locally

```r
shiny::runApp("app_costly_intel_study.R")
```

## Configuration

Before deploying, open the file and fill in your Supabase credentials near the top:

```r
SUPABASE_URL      <- "https://your-project.supabase.co"
SUPABASE_ANON_KEY <- "your-anon-key"
```

If left blank the app runs normally and data is saved to local CSV files in a folder called `study_logs` instead.

## Deployment

The app is designed to be deployed on shinyapps.io. With the `rsconnect` package installed, run:

```r
rsconnect::deployApp("app_costly_intel_study.R")
```

## Admin access

Enter `admin` as the participant code at the start screen to access the full settings panel.
