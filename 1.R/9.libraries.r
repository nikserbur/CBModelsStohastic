`%!in%` <- Negate(`%in%`)
libs <- c("qs","data.table","plotly","tidyverse","stats4","fitdistrplus","stabledist","scoringRules","readxl","alphastable")

new.packages <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

purrr::walk(libs, ~library(.x, character.only = TRUE))

