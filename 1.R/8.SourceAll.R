# Централизованный соурс всех разработанных функций
# Стохастическое моделирование процентных ставок ЦБ
# Последний апдейт 25.12.2024

#### Подключение библиотек и всех функций проекта ####
# Each file is sourced separately. Missing files and errors are reported
# instead of being hidden, and the success message is printed only when
# every file has been loaded.
local({
  # all_fun.R sources 1.R/9.libraries.r itself, so it is not listed here
  # (listing both would install-check and attach every package twice).
  project_files <- c(
    "1.R/1.functions/all_fun.R"
  )

  failed_files <- character(0)
  for (file_path in project_files) {
    if (!file.exists(file_path)) {
      message("[8.SourceAll.R] Missing file: ", file_path)
      failed_files <- c(failed_files, file_path)
      next
    }
    is_loaded <- tryCatch({
      source(file_path)
      TRUE
    }, error = function(e) {
      message("[8.SourceAll.R] Error while sourcing ", file_path, ": ",
              conditionMessage(e))
      FALSE
    })
    if (!is_loaded) failed_files <- c(failed_files, file_path)
  }

  if (length(failed_files) == 0) {
    cat("All libraries and project functions were loaded successfully.\n")
  } else {
    warning("8.SourceAll.R: ", length(failed_files), " of ", length(project_files),
            " file(s) were not loaded: ", paste(failed_files, collapse = ", "),
            call. = FALSE)
  }
})
