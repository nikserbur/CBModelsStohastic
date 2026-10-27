# 9.libraries.r — Единый файл подключения всех библиотек проекта
# Стохастическое моделирование процентных ставок ЦБ

#### Константы ####
`%!in%` <- Negate(`%in%`)

#### Список всех необходимых библиотек ####
libs <- c(
  # Основные библиотеки для работы с данными
  "data.table",        # Быстрая обработка данных
  "tidyverse",         # Коллекция для data science
  "magrittr",          # Pipe оператор %>%
  "RColorBrewer",
  # Визуализация
  "plotly",            # Интерактивные графики
  "visNetwork",        # Интерактивные сетевые диаграммы
  "viridis",           # Colour palettes (5.HiddenMarkingsModels.R)

  # Статистика и моделирование
  "stats",             # Базовая статистика R
  "stats4",            # Дополнительные статистические функции
  "moments",           # Моменты распределений (kurtosis, skewness)
  "MASS",              # Статистические функции и распределения

  # Распределения и стохастические процессы
  "extraDistr",        # Дополнительные распределения
  "fitdistrplus",      # Подгонка распределений
  "stabledist",        # α-stable распределения
  "alphastable",       # Альтернативная библиотека для α-stable
  "actuar",            # dpareto(shape, scale) for fitdist(..., "pareto") in Main.R; must be attached after extraDistr
  "Pareto",            # Pareto random draws (10.DashboardData.R examples)

  # Оценка качества моделей
  "scoringRules",      # CRPS, PIT и другие метрики
  "nortest",           # Anderson-Darling normality test (10.DashboardData.R)

  # Работа с файлами
  # "qs" is not attached here: it was archived on CRAN on 2026-01-17, so
  # install.packages("qs") fails. It is only needed for the optional
  # save_qs_path argument of BuildDistDashboard() in 10.DashboardData.R.
  "readxl",            # Чтение Excel файлов

  # Дополнительные библиотеки для расширенной функциональности
  "zoo",               # Временные ряды
  "forecast",          # Прогнозирование
  "lubridate",         # Работа с датами
  "xts",               # Time series objects (Demo_Volatility_Analysis.R)
  "rugarch",           # GARCH(1,1) benchmark (Demo_Volatility_Analysis.R)
  "highfrequency",     # HAR model (Demo_Volatility_Analysis.R)
  "vars",              # VAR models (PredictTimeSeriesWithSVD() in all_fun.R)
  "Rcpp"               # Compiling 2.Cpp/*.cpp via Rcpp::sourceCpp()
)

#### Автоматическая установка недостающих пакетов ####
cat("Проверка и установка необходимых библиотек...\n")

new.packages <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  cat("Устанавливаем недостающие пакеты:", paste(new.packages, collapse = ", "), "\n")
  install.packages(new.packages, dependencies = TRUE)
} else {
  cat("Все необходимые пакеты уже установлены\n")
}

#### Подключение всех библиотек ####
cat("Подключение библиотек...\n")
purrr::walk(libs, ~{
  suppressPackageStartupMessages(library(.x, character.only = TRUE))
  cat("✓", .x, "\n")
})

#### Дополнительные настройки ####
# Настройки для data.table
setDTthreads(0)  # Использовать все доступные ядра

# Настройки для plotly
options(plotly.browser = TRUE)  # Открывать графики в браузере

cat("Настройки завершены. Система готова к работе.\n")

