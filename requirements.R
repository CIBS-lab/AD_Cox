packages <- c("survminer", "survival", "MatchIt", "tibble", "sjPlot", 
              "naniar", "finalfit", "broom", "coxme", "stargazer", 
              "survcomp", "arsenal", "dplyr", "ggplot2", "scater", "ggstatsplot")

installed_packages <- rownames(installed.packages())
missing_packages <- packages[!packages %in% installed_packages]

if(length(missing_packages)) {
  install.packages(missing_packages, dependencies = TRUE)
}

# Load all packages
lapply(packages, library, character.only = TRUE)
