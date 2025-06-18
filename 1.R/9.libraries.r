`%!in%` <- Negate(`%in%`)
libs <- c("qs","data.table")

new.packages <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(lapply(libs, require, character.only = TRUE))

