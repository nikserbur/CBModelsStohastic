# Main.R — Основной файл комплексного стохастического моделирования ставок ЦБ
# Проект: Стохастическое моделирование процентных ставок с жирными хвостами
# Интеграция: HAR модели, HMM, Gradient Boosting, Taleb R² методики
# Разработчик: Бураков Никита Сергеевич

#' Main.R — Интегрированная стохастическая модель ставок ЦБ РФ
#' Методология: SVD-VAR (детерминистическая траектория) + HAR волатильность + HMM (7 режимов, переходы по Маркову) + расширенное правило Тейлора с экспоненциальным затуханием влияний + α-stable ошибки.
#' Особенности: LM-интерполяция пропусков; смешивание base/shock = 0.7/0.3; φ_AR ≤ 0.85; конвергенция нейтральной к 7% с λ=0.35; режимный сдвиг — смесь средних эмиссий с прогнозными вероятностями HMM; демонстрация «лучший/худший» пути в beautiful_plot.
#' Ожидаемый ориентир: медиана ставки на горизонте 19 мес ~ 7–8% при сохранении хвостовых рисков.
#' Запуск из корня репозитория: Rscript 1.R/Main.R
#### 0. Подключение всех библиотек и функций ####
source("1.R/8.SourceAll.R")

#### 1. Константы и параметры проекта ####
set.seed(123)
N_SIMULATIONS <- 10000
HORIZON <- 19
TARGET_INFLATION <- 4.0
NEUTRAL_RATE <- 7.0

# Параметры правила Тейлора ЦБ РФ (официальные веса)
TAYLOR_INFLATION_WEIGHT <- 1.5   # Реакция на отклонение инфляции от цели
TAYLOR_OUTPUT_WEIGHT <- 0.5      # Реакция на разрыв выпуска
TAYLOR_SMOOTHING <- 0.8          # Параметр сглаживания ставок

performance_profile <- data.table(
  Module = character(),
  Time_sec = numeric(),
  Complexity = character()
)

log_performance <- function(module, start_time, complexity) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  performance_profile <<- rbind(
    performance_profile,
    data.table(Module = module, Time_sec = elapsed, Complexity = complexity)
  )
}

regularisation_log <- data.table(
  Step = character(),
  Constraint = character(),
  LogLikelihood = numeric(),
  Notes = character()
)

append_reg_log <- function(step, constraint, log_likelihood = NA_real_, notes = "") {
  constraint <- enc2utf8(constraint)
  notes <- enc2utf8(notes)
  Encoding(constraint) <- "UTF-8"
  Encoding(notes) <- "UTF-8"
  regularisation_log <<- rbind(
    regularisation_log,
    data.table(
      Step = step,
      Constraint = constraint,
      LogLikelihood = log_likelihood,
      Notes = notes
    )
  )
}

#### 2. Загрузка и подготовка данных ####

# prep_data() (загрузка и объединение файлов из 3.Data) определена в
# 1.R/1.functions/all_fun.R, чтобы её могли использовать и другие скрипты.

input_data <- prep_data()

#### 3. Разделение данных на обучающую и прогнозную выборки ####

# Определение последней даты обучения - берем максимальную дату по всем показателям
last_available_dates <- input_data[, .(
  last_ruonia = max(date[!is.na(ruonia_rate)]),
  last_inflation = max(date[!is.na(cpi_inflation)]),
  last_gdp = max(date[!is.na(gdp_level)]),
  last_volatility = max(date[!is.na(volatility)])
)]

# Используем последнюю доступную дату среди ключевых показателей
last_train_date <- min(last_available_dates$last_ruonia,
                      last_available_dates$last_inflation,
                      last_available_dates$last_gdp)

# Применение SVD-VAR экстраполяции
macro_forecasts <- PredictTimeSeriesWithSVD(
  input_dataset = input_data[date <= last_train_date],#-c("ruonia_rate","key_rate")
  prediction_steps = HORIZON,
  temporal_var = "date"
)

# Создание прогнозных дат и базового прогнозного датасета
forecast_dates <- seq.Date(last_train_date + months(1), last_train_date + months(HORIZON), by = "month")


forecast_data <- macro_forecasts$forecast %>%
  .[,date := forecast_dates] %>%
  .[,-c("V1")]


all_data <- rbindlist(list(input_data[date <= last_train_date],forecast_data),use.names = TRUE) %>%
  setorder(date) %>%
  .[, inflation_gap := cpi_inflation - TARGET_INFLATION] %>%
  .[, policy_tightness := key_rate - cpi_inflation] %>%
  .[, inflation_expectations_gap := expected_inflation - TARGET_INFLATION] %>%
  .[, gdp_growth := (gdp_level / lag(gdp_level, 1) - 1) * 100] %>%
  .[, credit_to_gdp := credit_impulse / gdp_level * 100] %>%
  .[, key_rate_change := key_rate - lag(key_rate, 1)] %>%
  .[, gdp_trend := loess(gdp_level ~ as.numeric(date), span = 0.3)$fitted] %>%
  .[, gdp_gap := (gdp_level - gdp_trend) / gdp_trend * 100] %>%
  .[, inflation_deviation := cpi_inflation - TARGET_INFLATION] %>%
  .[, rate_change := ruonia_rate - lag(ruonia_rate, 1)] %>%
  # Создание волатильности ставок (realized volatility)
  .[, rate_volatility := frollapply(ruonia_rate, 6, function(x) {
    if (sum(!is.na(x)) < 3) return(NA)
    sd(x, na.rm = TRUE)
  }, align = "right")] %>%
  .[, rate_volatility := ifelse(is.na(rate_volatility) | rate_volatility == 0,
                                mean(rate_volatility, na.rm = TRUE), rate_volatility)] %>%
  .[, asset_volatility := volatility] %>%
  .[, volatility_regime := ifelse(rate_volatility > quantile(rate_volatility, 0.75, na.rm = TRUE), "high",
                                  ifelse(rate_volatility < quantile(rate_volatility, 0.25, na.rm = TRUE), "low", "normal"))] %>%
  .[, (names(.)[sapply(., is.numeric) & names(.) != "date"]) := lapply(.SD, function(x) {
                                    if (any(is.na(x))) {
                                      tt <- 1:length(x)
                                      fit <- lm(x ~ tt)
                                      x_hat <- predict(fit, newdata = data.frame(tt = tt))
                                      ix <- which(is.na(x))
                                      x[ix] <- x_hat[ix]
                                      x
                                    } else x
                                  }), .SDcols = names(.)[sapply(., is.numeric) & names(.) != "date"]]

# Прогнозный датасет с интеллектуальными прогнозами
forecast_data <- copy(all_data)[date > last_train_date,-c("ruonia_rate","key_rate")] %>%
  .[, har_volatility := NA_real_] %>%  # Будет заполнено HAR прогнозом
  .[, regime_probs := NA_real_]       # Будут заполнены HMM прогнозами


train_data <-copy(all_data)[date <= last_train_date]

#### 4. ЭТАП 1: Калибровка стохастических распределений по методу максимального правдоподобия ####

# Калибровка распределений
stochastic_params <- CalibrateStochasticDistributions(train_data)

# Проверяем калибровку макропеременных
for (var_name in names(stochastic_params$macro_distributions)) {
  dist_info <- stochastic_params$macro_distributions[[var_name]]
  cat("   ", var_name, " - распределение:", dist_info$distribution,
      ", параметры:", paste(round(unlist(dist_info$params), 4), collapse = ", "), "\n")
}

#### 5. ЭТАП 2: HAR модель волатильности ставок ####

volatility_series <- train_data$rate_volatility[!is.na(train_data$rate_volatility)]
dates <- seq.Date(as.Date("2016-01-01"), by = "month", length.out = length(volatility_series)-1)
volatility <- as.xts(volatility_series[-length(volatility_series)], order.by = dates)

# Параметры HAR
periods <- c(5,22)
horizon <- HORIZON

# Прогноз на 12 месяцев (рекурсивно)
forecast_vals <- numeric(horizon)

for (i in 1:horizon) {
  # Обновляем модель с текущими данными
  har_model <- highfrequency::HARmodel(
    data = volatility,
    periods = periods,
    RVest = c("rCov"),
    type = "HAR",
    h = 1,
    transform = NULL,
    inputType = "RM"
  )

  # Прогноз на 1 шаг вперёд
  pred <- predict(har_model)
  forecast_vals[i] <- as.numeric(pred)

  # Добавляем предсказанное значение к временным рядам
  next_date <- end(volatility) + months(1)
  volatility <- rbind(volatility, xts(pred, order.by = next_date))
}

# Добавление HAR прогноза в прогнозный датасет
forecast_data[, har_volatility := forecast_vals]

#### 6. ЭТАП 3: Улучшенные HMM режимы с марковскими цепями ####

# Обучение HMM на всех доступных данных (комплексный индикатор для режимов)
hmm_training_data <- copy(train_data) %>%
  # Создание композитного индикатора для определения режимов (взвешенная комбинация)
  .[, regime_indicators := {
    important_cols <- c("gdp_gap", "B002", "credit_to_gdp", "budget_impulse", "policy_tightness","ruonia_rate")
    important_weight <- 0.07
    all_numeric_cols <- names(.)[sapply(., is.numeric) & names(.) != "date"]

    other_cols <- setdiff(all_numeric_cols, important_cols)

    remaining_weight <- 1 - length(important_cols) * important_weight
    other_weight <- remaining_weight / length(other_cols)

    important_part <- paste0(important_weight, " * scale(", important_cols, ")[,1]", collapse = " + ")
    other_part <- paste0(other_weight, " * scale(", other_cols, ")[,1]", collapse = " + ")

    full_formula <- paste0(important_part, " + ", other_part)

    eval(parse(text = full_formula))
  }]

# Инициализация параметров HMM
hmm_params <- InitializeHmmParameters(n_states = 7, observations = hmm_training_data$regime_indicators)
# Обучение HMM
timer_hmm_est <- Sys.time()
hmm_model <- EstimateHmmParameters(
  observations = hmm_training_data$regime_indicators,
  initial_params = hmm_params,
  max_iterations = 1000,
  tolerance = 1e-6
)
log_performance("HMM estimation", timer_hmm_est, "O(m^2 * n)")
append_reg_log(
  "HMM estimation",
  "Диагональные элементы >= 0.8, вероятности >= 0",
  if (!is.null(hmm_model$log_likelihood)) hmm_model$log_likelihood else NA_real_,
  paste("Итерации:", hmm_model$iterations)
)

# Декодирование скрытых состояний
hidden_states <- ViterbiDecoding(
  observations = hmm_training_data$regime_indicators,
  hmm_params = hmm_model
)

# Создание обучающего датасета с режимами
train_data_with_regimes <- copy(train_data) %>%
  .[, regime := factor(hidden_states, levels = 1:7, labels = paste0(1:7))] %>%
  .[, regime_1 := as.numeric(regime == "1")] %>%
  .[, regime_2 := as.numeric(regime == "2")] %>%
  .[, regime_3 := as.numeric(regime == "3")] %>%
  .[, regime_4 := as.numeric(regime == "4")] %>%
  .[, regime_5 := as.numeric(regime == "5")] %>%
  .[, regime_6 := as.numeric(regime == "6")] %>%
  .[, regime_7 := as.numeric(regime == "7")] %>%
  # Добавляем прогнозную волатильность HAR для обучающей выборки
  .[, har_volatility := rate_volatility]

# Прогнозирование режимов для прогнозного периода
# hmm_model$transition_matrix[7,7] <- 0.92
# hmm_model$transition_matrix[7,6] <- 0.04
# hmm_model$transition_matrix[7,5] <- 0.02
# hmm_model$transition_matrix[7,4] <- 0.01
# hmm_model$transition_matrix[7,3] <- 0.01
regime_forecasts <- PredictRegimeProbabilities(hmm_model, HORIZON, current_state = 7)

# Добавление прогнозов режимов в прогнозный датасет
forecast_data[, regime_1 := regime_forecasts[, "1"]]
forecast_data[, regime_2 := regime_forecasts[, "2"]]
forecast_data[, regime_3 := regime_forecasts[, "3"]]
forecast_data[, regime_4 := regime_forecasts[, "4"]]
forecast_data[, regime_5 := regime_forecasts[, "5"]]
forecast_data[, regime_6 := regime_forecasts[, "6"]]
forecast_data[, regime_7 := regime_forecasts[, "7"]]

regime_dist <- table(train_data_with_regimes$regime)
for (regime_name in names(regime_dist)) {
  cat("     ", regime_name, ":", regime_dist[regime_name], "наблюдений (",
      round(regime_dist[regime_name]/sum(regime_dist)*100, 1), "%)\n")
}

#### 7. ЭТАП 4: Калибровка параметров правила Тейлора методом максимального правдоподобия ####

# Калибровка параметров правила Тейлора на исторических данных
timer_taylor_calib <- Sys.time()
taylor_calibration <- CalibrateTaylorRuleParameters(train_data)
log_performance("Taylor calibration", timer_taylor_calib, "O(n)")
append_reg_log(
  "Taylor calibration",
  "ρ∈[0.1,0.99], φπ>0, φy∈[-2,2]",
  if (!is.null(taylor_calibration$log_likelihood)) taylor_calibration$log_likelihood else NA_real_,
  paste("Convergence:", taylor_calibration$convergence)
)

# Обновление параметров правила Тейлора на основе калибровки
if (taylor_calibration$convergence) {
  # Используем калиброванные параметры
  TAYLOR_NEUTRAL_RATE <- taylor_calibration$parameters["r_star"]
  TAYLOR_INFLATION_WEIGHT <- taylor_calibration$parameters["phi_pi"]
  TAYLOR_OUTPUT_WEIGHT <- taylor_calibration$parameters["phi_y"]
  # TAYLOR_SMOOTHING <- taylor_calibration$parameters["rho"]
  TAYLOR_ERROR_SIGMA <- taylor_calibration$parameters["sigma"]

  cat("   Используем калиброванные параметры правила Тейлора:\n")
  cat("     r* (нейтральная ставка):", round(TAYLOR_NEUTRAL_RATE, 3), "\n")
  cat("     φ_π (вес инфляции):", round(TAYLOR_INFLATION_WEIGHT, 3), "\n")
  cat("     φ_y (вес выпуска):", round(TAYLOR_OUTPUT_WEIGHT, 3), "\n")
  cat("     ρ (сглаживание):", round(TAYLOR_SMOOTHING, 3), "\n")
  cat("     σ (волатильность ошибок):", round(TAYLOR_ERROR_SIGMA, 3), "\n")

  # Корректировка нейтральной ставки для более реалистичного прогноза
  NEUTRAL_RATE_ADJUSTMENT <- 2.0  # Снижение на 2%
  TAYLOR_NEUTRAL_RATE <- TAYLOR_NEUTRAL_RATE - NEUTRAL_RATE_ADJUSTMENT
  cat("   Скорректированная нейтральная ставка:", round(TAYLOR_NEUTRAL_RATE, 3), "\n")
} else {
  cat("   Используем исходные параметры правила Тейлора (калибровка не удалась):\n")
  TAYLOR_NEUTRAL_RATE <- NEUTRAL_RATE
  TAYLOR_ERROR_SIGMA <- 1.0
}

# Официальное правило Тейлора ЦБ РФ с калиброванными параметрами:
# i_t = ρ * i_{t-1} + (1-ρ) * [r* + φ_π * (π_t - π*) + φ_y * y_t] + ε_t
# где параметры калиброваны методом максимального правдоподобия

# Базовое правило Тейлора ЦБ РФ с калиброванными параметрами
# Используем калиброванные параметры вместо OLS оценки
taylor_cb_model <- list(
  coefficients = c(
    "(Intercept)" = TAYLOR_NEUTRAL_RATE * (1 - TAYLOR_SMOOTHING),
    "lag_rate" = TAYLOR_SMOOTHING,
    "inflation_deviation" = TAYLOR_INFLATION_WEIGHT * (1 - TAYLOR_SMOOTHING),
    "gdp_gap" = TAYLOR_OUTPUT_WEIGHT * (1 - TAYLOR_SMOOTHING)
  ),
  fitted_values = TAYLOR_SMOOTHING * lag(train_data_with_regimes$ruonia_rate, 1) +
                 (1 - TAYLOR_SMOOTHING) * (TAYLOR_NEUTRAL_RATE +
                 TAYLOR_INFLATION_WEIGHT * train_data_with_regimes$inflation_deviation +
                 TAYLOR_OUTPUT_WEIGHT * train_data_with_regimes$gdp_gap)
)

# Расчет остатков для калиброванной модели
taylor_residuals <- train_data_with_regimes$ruonia_rate - taylor_cb_model$fitted_values
taylor_residuals <- taylor_residuals[!is.na(taylor_residuals)]

# Извлечение коэффициентов правила Тейлора (калиброванные)
rho_smoothing <- TAYLOR_SMOOTHING  # Параметр сглаживания
phi_inflation <- TAYLOR_INFLATION_WEIGHT  # Реакция на инфляцию
phi_output <- TAYLOR_OUTPUT_WEIGHT     # Реакция на выпуск

# Корректная оценка α параметра для остатков
alpha_estimate <- EstimateAlphaParameter(abs(taylor_residuals), method = "hill")
alpha_corrected <- min(max(alpha_estimate, 1.1), 4.0)  # Ограничиваем в разумных пределах для финансовых данных

# Параметр α остатков: alpha_corrected

#### 8. ЭТАП 5: Расширенное правило Тейлора с новыми индикаторами ####

# Новое расширенное правило Тейлора с дополнительными индикаторами:
# i_t = ρ*i_{t-1} + (1-ρ)*[r* + φ_π*(π_t-π*) + φ_y*y_t + φ_v*vol_t + Σφ_j*regime_j + φ_πe*(πe_t-π*)] + ε_t

timer_extended_taylor <- Sys.time()
extended_taylor_model <- lm(
  ruonia_rate ~ I(lag(ruonia_rate, 1)) + inflation_deviation + gdp_gap +
                har_volatility + regime_7 +
                expected_inflation + credit_to_gdp,
  data = train_data_with_regimes
)
log_performance("Extended Taylor regression", timer_extended_taylor, "O(n)")

extended_residuals <- residuals(extended_taylor_model)
print(summary(extended_taylor_model))

# Расчет скорректированного R² с учетом параметров распределения
extended_r2_result <- CalculateRobustR2(
  y_actual = train_data_with_regimes$ruonia_rate,
  y_predicted = fitted(extended_taylor_model),
  alpha_tail = alpha_corrected
)

# Shapley анализ важности новых переменных (исключая зависимые переменные)
# Проверяем какие переменные доступны в обучающих данных
available_vars <- names(train_data_with_regimes)
independent_variables <- c("inflation_deviation", "gdp_gap", "har_volatility",
                          paste0("regime_", 1:7), "expected_inflation", "credit_to_gdp")

# Оставляем только доступные переменные
independent_variables <- independent_variables[independent_variables %in% available_vars]

shapley_data <- train_data_with_regimes[, c("ruonia_rate", independent_variables), with = FALSE]
shapley_r2 <- ComputeShapleyR2(
  data = shapley_data,
  target_var = "ruonia_rate",
  alpha_tail = alpha_corrected
)
print(round(shapley_r2, 3))

#### 8.5. Сравнение базового и расширенного правила Тейлора ####

# Расчет метрик для обеих моделей
# Для калиброванной модели используем расчет R² вручную
base_r2 <- 1 - sum(taylor_residuals^2) / sum((train_data_with_regimes$ruonia_rate[!is.na(taylor_cb_model$fitted_values)] - mean(train_data_with_regimes$ruonia_rate[!is.na(taylor_cb_model$fitted_values)]))^2)
extended_r2 <- summary(extended_taylor_model)$r.squared

# Для калиброванной модели AIC/BIC не применимы (не OLS модель)
base_aic <- NA
extended_aic <- AIC(extended_taylor_model)

base_bic <- NA
extended_bic <- BIC(extended_taylor_model)

# Сравнительная таблица
taylor_comparison <- data.table(
  Model = c("Базовое правило Тейлора", "Расширенное правило Тейлора"),
  R_squared = c(base_r2, extended_r2),
  AIC = c(base_aic, extended_aic),
  BIC = c(base_bic, extended_bic),
  Variables = c("inflation_deviation + gdp_gap",
                "inflation_deviation + gdp_gap + har_volatility + regime_1..regime_7 + expected_inflation + credit_to_gdp"),
  Taleb_R2 = c(CalculateRobustR2(train_data_with_regimes$ruonia_rate[!is.na(taylor_cb_model$fitted_values)], taylor_cb_model$fitted_values[!is.na(taylor_cb_model$fitted_values)], alpha_corrected)$r2_taleb_exact,
              extended_r2_result$r2_taleb_exact)
)

print(taylor_comparison[, .(Model, R_squared, Taleb_R2, AIC, BIC)])
#### 9. Интегрированная стохастическая модель ####

# Функция симуляции одного пути simulate_single_path() вложена внутрь IntegratedStochasticCBModel
IntegratedStochasticCBModel <- function(train_data, forecast_data, stochastic_params,
                                       taylor_models, hmm_model, horizon = 12, n_simulations = 10000) {

  # Запуск интегрированной стохастической модели

  # Извлечение параметров
  base_model <- taylor_models$base
  extended_model <- taylor_models$extended
  vol_params <- stochastic_params$volatility
  stable_params <- stochastic_params$stable_errors

  # Начальные условия
  last_rate <- train_data$ruonia_rate[nrow(train_data)]
  last_vol <- train_data$rate_volatility[nrow(train_data)]

  # Предрасчёт коэффициентов расширенной модели и якоря нейтральной ставки
  extended_coef <- coef(extended_model)
  rho_ext <- if ("I(lag(ruonia_rate, 1))" %in% names(extended_coef)) as.numeric(extended_coef[["I(lag(ruonia_rate, 1))"]]) else TAYLOR_SMOOTHING
  rho_ext <- max(0.4, min(0.85, rho_ext))
  b_infl <- if ("inflation_deviation" %in% names(extended_coef)) as.numeric(extended_coef[["inflation_deviation"]]) else 0
  b_gdp <- if ("gdp_gap" %in% names(extended_coef)) as.numeric(extended_coef[["gdp_gap"]]) else 0
  b_vol <- if ("har_volatility" %in% names(extended_coef)) as.numeric(extended_coef[["har_volatility"]]) else 0
  b_exp <- if ("expected_inflation" %in% names(extended_coef)) as.numeric(extended_coef[["expected_inflation"]]) else 0
  b_credit <- if ("credit_to_gdp" %in% names(extended_coef)) as.numeric(extended_coef[["credit_to_gdp"]]) else 0
  regime_means <- hmm_model$emission_params$mean
  NEUTRAL_TARGET <- NEUTRAL_RATE
  lambda_neutral <- 0.35
  # Затухания влияния предикторов по горизонту
  kappa_pi <- 0.25; kappa_y <- 0.15; kappa_vol <- 0.20; kappa_exp <- 0.20; kappa_cr <- 0.10

  # Вложенная функция симуляции одного пути (видит локальные переменные)
  simulate_single_path <- function(sim_id) {
    path <- numeric(horizon)
    vol_path <- numeric(horizon)
    regime_path <- numeric(horizon)

    current_rate <- last_rate
    current_vol <- last_vol
    # Инициализация режимом по последнему историческому состоянию (Витерби)
    current_regime <- as.integer(hidden_states[length(hidden_states)])

    for (t in 1:horizon) {
      # Режимы: марковские переходы от предыдущего состояния
      transition_probs <- hmm_model$transition_matrix[current_regime, ]
      current_regime <- sample(1:7, 1, prob = transition_probs)

      # Режимные веса (прогнозные вероятности) и восстановление уровня ставки по режимам
      regime_cols <- paste0("regime_", 1:length(regime_means))
      regime_probs_t <- if (all(regime_cols %in% names(forecast_data))) as.numeric(forecast_data[t, ..regime_cols]) else rep(0, length(regime_means))
      regime_probs_t <- regime_probs_t / pmax(sum(regime_probs_t), 1e-9)
      regime_mixture_shift <- sum(regime_means[seq_along(regime_probs_t)] * regime_probs_t)

      # Волатильность: HAR + стохастическая компонента (Хестон)
      har_vol_forecast <- forecast_data$har_volatility[t]
      vol_innovation <- rnorm(1, 0, vol_params$sigma_v * sqrt(current_vol))
      stochastic_vol_component <- vol_params$kappa * (vol_params$theta - current_vol) + vol_innovation
      current_vol <- 0.6 * har_vol_forecast + 0.4 * (current_vol + stochastic_vol_component)
      if (!is.finite(current_vol) || is.na(current_vol) || current_vol <= 0) current_vol <- 1e-4

      # Макро: базовые прогнозы + калиброванные шоки (сдвиг в пользу базы)
      base_inflation <- forecast_data$inflation_deviation[t]
      base_gdp <- forecast_data$gdp_gap[t]

      inflation_dist <- stochastic_params$macro_distributions$cpi_inflation
      if (!is.null(inflation_dist) && inflation_dist$distribution == "pareto") {
        raw_shock <- extraDistr::rpareto(1, a = inflation_dist$params[[1]], b = inflation_dist$params[[2]])
        inflation_shock <- raw_shock - inflation_dist$params[[2]]
      } else if (!is.null(inflation_dist) && inflation_dist$distribution == "t_student") {
        inflation_shock <- rt(1, df = inflation_dist$params[[3]]) * inflation_dist$params[[2]] + inflation_dist$params[[1]]
      } else {
        inflation_shock <- rnorm(1, 0, 0.5)
      }

      gdp_dist <- stochastic_params$macro_distributions$gdp_gap
      if (!is.null(gdp_dist) && gdp_dist$distribution == "t_student") {
        gdp_shock <- rt(1, df = gdp_dist$params[[3]]) * gdp_dist$params[[2]] + gdp_dist$params[[1]]
      } else if (!is.null(gdp_dist) && gdp_dist$distribution == "pareto") {
        raw_shock <- extraDistr::rpareto(1, a = gdp_dist$params[[1]], b = gdp_dist$params[[2]])
        gdp_shock <- raw_shock - gdp_dist$params[[2]]
      } else {
        gdp_shock <- rnorm(1, 0, 0.3)
      }

      w_base <- 0.7; w_shock <- 0.3
      current_inflation_dev <- w_base * base_inflation + w_shock * inflation_shock
      current_gdp_gap <- w_base * base_gdp + w_shock * gdp_shock

      # Доп. предикторы для расширенного правила (разрывы/зажигание тренда)
      current_expected_gap <- if ("expected_inflation" %in% names(forecast_data)) (forecast_data$expected_inflation[t] - TARGET_INFLATION) else 0
      current_credit <- if ("credit_to_gdp" %in% names(forecast_data)) forecast_data$credit_to_gdp[t] else 0

      # Нейтральная ставка с конвергенцией к целевому уровню и смесью режимного уровня
      neutral_t <- TAYLOR_NEUTRAL_RATE + (NEUTRAL_TARGET - TAYLOR_NEUTRAL_RATE) * (1 - exp(-lambda_neutral * t)) + regime_mixture_shift

      # Расширенное правило Тейлора со затухающими влияниями предикторов
      decay_pi <- exp(-kappa_pi * t)
      decay_y <- exp(-kappa_y * t)
      decay_vol <- exp(-kappa_vol * t)
      decay_exp <- exp(-kappa_exp * t)
      decay_cr <- exp(-kappa_cr * t)

      linear_component <- neutral_t +
        (b_infl * decay_pi) * current_inflation_dev +
        (b_gdp * decay_y) * current_gdp_gap +
        (b_vol * decay_vol) * current_vol +
        (b_exp * decay_exp) * current_expected_gap +
        (b_credit * decay_cr) * current_credit

      # Сглаживание по AR-коэффициенту расширенной модели
      target_rate <- rho_ext * current_rate + (1 - rho_ext) * linear_component

      # α-stable ошибки (жирные хвосты)
      alpha_for_stable <- min(max(alpha_corrected, 1.01), 1.99)
      gamma_stable <- max(0.5 * sqrt(current_vol), 1e-4)
      if (!is.finite(gamma_stable) || is.na(gamma_stable) || gamma_stable <= 0) gamma_stable <- 1e-4
      stable_error <- stabledist::rstable(1, alpha = alpha_for_stable, beta = 0, gamma = gamma_stable, delta = 0)
      if (!is.finite(stable_error) || is.na(stable_error)) stable_error <- 0

      path[t] <- max(0.1, min(40, target_rate + stable_error))

      current_rate <- path[t]
      vol_path[t] <- current_vol
      regime_path[t] <- current_regime
    }

    list(rates = path, volatility = vol_path, regimes = regime_path)
  }

  # Monte Carlo симуляция
  simulation_results <- replicate(n_simulations, simulate_single_path(), simplify = FALSE)

  # Извлечение результатов
  simulation_paths <- sapply(simulation_results, function(x) x$rates)
  volatility_paths <- sapply(simulation_results, function(x) x$volatility)
  regime_paths <- sapply(simulation_results, function(x) x$regimes)

  # Расчет квантилей
  quantiles <- apply(simulation_paths, 1, function(x) {
    quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
  })

  # Анализ режимов в симуляциях
  regime_frequencies <- apply(regime_paths, 1, function(x) {
    table(factor(x, levels = 1:7)) / length(x)
  })

  results <- list(
    simulations = simulation_paths,
    volatility_simulations = volatility_paths,
    regime_simulations = regime_paths,
    quantiles = t(quantiles),
    regime_frequencies = t(regime_frequencies),
    parameters = list(
      stochastic = stochastic_params,
      taylor_base = base_model$coefficients,
      taylor_extended = coef(extended_model),
      hmm = hmm_model
    ),
    forecast_dates = forecast_data$date,
    model_formulas = list(
      base_taylor = "i_t = ρ*i_{t-1} + (1-ρ)*[r* + φ_π*(π_t-π*) + φ_y*y_t] + ε_t",
      extended_taylor = "i_t = φ_ar*i_{t-1} + (1-φ_ar)*[r*_t + β'X_t] + ε_t",
      stochastic_vol = "dlog(V_t) = κ(θ - log(V_t))dt + σ_v dW_t",
      stable_errors = "ε_t ~ S_α(β, γ√V_t, δ)",
      har_integration = "vol_t = 0.6 * HAR_forecast_t + 0.4 * stochastic_evolution_t",
      hmm_integration = "regime_t ~ Markov_transition (смешивание уровней через вероятности)",
      combined_approach = "Детерминистические прогнозы (HAR+HMM+SVD-VAR) + стохастические шоки (Pareto+α-stable)"
    )
  )

  return(results)
}

# Подготовка моделей для интегрированной симуляции
taylor_models <- list(
  base = taylor_cb_model,
  extended = extended_taylor_model
)

# Запуск интегрированной модели с исправленной симуляцией
timer_integrated_model <- Sys.time()
final_integrated_model <- IntegratedStochasticCBModel(
  train_data = train_data_with_regimes,
  forecast_data = forecast_data,
  stochastic_params = stochastic_params,
  taylor_models = taylor_models,
  hmm_model = hmm_model,
  horizon = HORIZON,
  N_SIMULATIONS
)
log_performance("Integrated stochastic simulation", timer_integrated_model, "O(H * N_sim)")

#### 10. Анализ результатов и визуализация ####

# Диагностика модели
diagnostics <- data.table(
  horizon = 1:nrow(final_integrated_model$quantiles),
  date = forecast_data$date,
  q05 = final_integrated_model$quantiles[, 1],
  q10 = final_integrated_model$quantiles[, 2],
  q25 = final_integrated_model$quantiles[, 3],
  median = final_integrated_model$quantiles[, 4],
  q75 = final_integrated_model$quantiles[, 5],
  q90 = final_integrated_model$quantiles[, 6],
  q95 = final_integrated_model$quantiles[, 7]
)

# Анализ волатильности и режимов
volatility_analysis <- list(
  mean_vol_path = apply(final_integrated_model$volatility_simulations, 1, mean),
  vol_quantiles = apply(final_integrated_model$volatility_simulations, 1, function(x) {
    quantile(x, probs = c(0.25, 0.5, 0.75))
  }),
  regime_stability = apply(final_integrated_model$regime_frequencies, 1, function(x) {
    max(x) - min(x)  # Мера концентрации режимов
  })
)

quantile_ok <- all(apply(final_integrated_model$quantiles, 1, function(row) {
  all(diff(row) >= -1e-8)
}))
vol_ok <- all(volatility_analysis$mean_vol_path >= -1e-8, na.rm = TRUE)
regime_diff <- diff(volatility_analysis$regime_stability)
regime_trend <- if (length(regime_diff)) all(regime_diff <= 1e-8 | is.na(regime_diff)) else TRUE

validation_report <- data.table(
  Test = c(
    "Порядок квантилей прогнозов",
    "Неотрицательность волатильности",
    "Нисходящая концентрация режимов"
  ),
  Status = c(
    if (quantile_ok) "PASS" else "WARN",
    if (vol_ok) "PASS" else "WARN",
    if (regime_trend) "PASS" else "WARN"
  ),
  Runtime = c(NA_real_, NA_real_, NA_real_),
  Notes = c(
    if (quantile_ok) "Квантильные ленты упорядочены" else "Обнаружены нарушения порядка квантилей",
    if (vol_ok) "Средняя волатильность не опускается ниже нуля" else "Зафиксированы отрицательные значения волатильности",
    if (regime_trend) "Регимная концентрация снижается либо стабильна" else "Концентрация режимов увеличилась на отдельных горизонтах"
  )
)
validation_report[, `:=`(
  Test = enc2utf8(Test),
  Notes = enc2utf8(Notes)
)]
Encoding(validation_report$Test) <- "UTF-8"
Encoding(validation_report$Notes) <- "UTF-8"

cat("Медианный прогноз ставки через 12 месяцев:", round(diagnostics$median[12], 2), "%\n")

#### 10.1. Хвостовые метрики риска по горизонтам (VaR / CTE, верхний хвост) ####
# Квантили и средние хвоста считаются по всем путям для каждого горизонта
TAIL_CONFIDENCE <- 0.95
tail_risk <- ComputeTailRiskByHorizon(final_integrated_model$simulations, confidence = TAIL_CONFIDENCE)
tail_risk[, date := forecast_data$date]
print(tail_risk)

#### 10.2. Оценка прогнозного распределения (CRPS / PIT) ####
# Фактические месячные средние RUONIA за прогнозные месяцы берутся из того же файла,
# что и в prep_data(); последний (возможно, неполный) месяц файла не используется.
ruonia_daily <- read_xlsx("3.Data/RC_F11_01_2010_T04_07_2025.xlsx") %>%
  as.data.table() %>%
  .[, .(date = as.Date(DT), ruo)]
last_ruonia_month <- floor_date(max(ruonia_daily$date, na.rm = TRUE), "month")
realised_ruonia <- ruonia_daily[floor_date(date, "month") < last_ruonia_month,
                                .(ruonia_rate = mean(ruo, na.rm = TRUE)),
                                by = .(date = floor_date(date, "month"))]
realised_for_forecast <- realised_ruonia$ruonia_rate[match(forecast_data$date, realised_ruonia$date)]
forecast_evaluation <- EvaluatePredictiveSample(final_integrated_model$simulations, realised_for_forecast)
if (nrow(forecast_evaluation) == 0) {
  message("CRPS/PIT пропущены: в 3.Data нет полных месяцев RUONIA в прогнозном периоде (",
          format(min(forecast_data$date)), " – ", format(max(forecast_data$date)),
          "; последний полный месяц в файле — ", format(max(realised_ruonia$date)), ").")
} else {
  print(forecast_evaluation)
}

#### 11. Визуализация результатов ####
# Создание продвинутого графика
# Создание красивой визуализации Monte Carlo
beautiful_plot <- CreateBeautifulMonteCarloPlot(final_integrated_model, train_data_with_regimes, forecast_dates)

beautiful_plot

# Создание графика макроэкономических трендов
macro_trends_plot <- CreateMacroTrendsPlot(train_data_with_regimes, forecast_data)

macro_trends_plot

#### 11.25. Сравнительные визуализации моделей прогнозирования ####

# Функция для создания базового прогноза (детерминистические модели)
CreateBaselineForecasts <- function(train_data, forecast_horizon = 12) {

  # Подготовка данных для прогнозирования
  clean_data <- train_data[!is.na(ruonia_rate)]
  n_obs <- nrow(clean_data)

  # Метод 1: Простое среднее последних значений (консервативный подход)
  recent_window <- max(1, n_obs - 20):n_obs  # Последние 12 месяцев
  mean_forecast <- mean(clean_data$ruonia_rate[recent_window], na.rm = TRUE)

  # Метод 2: Линейная экстраполяция тренда
  if (n_obs >= 24) {
    # Используем последние 24 месяца для оценки тренда
    trend_window <- max(1, n_obs - 23):n_obs
    trend_data <- clean_data[trend_window]
    trend_model <- lm(ruonia_rate ~ seq_along(ruonia_rate), data = trend_data)
    trend_coef <- coef(trend_model)[2]  # Коэффициент тренда

    # Линейная экстраполяция
    linear_forecasts <- numeric(forecast_horizon)
    for (h in 1:forecast_horizon) {
      linear_forecasts[h] <- last(clean_data$ruonia_rate) + trend_coef * h
    }
  } else {
    # Если недостаточно данных, используем среднее
    linear_forecasts <- rep(mean_forecast, forecast_horizon)
  }

  # Метод 3: Улучшенное правило Тейлора с затуханием
  last_rate <- last(clean_data$ruonia_rate)
  last_inflation <- last(clean_data$inflation_deviation)
  last_gdp_gap <- ifelse("gdp_gap" %in% names(clean_data),
                         mean(clean_data$gdp_gap), -0.5)  # Дефолтное значение

  # Параметры затухания (более реалистичные)
  inflation_decay <- 0.15  # Скорость затухания инфляционного давления
  gdp_decay <- 0.08       # Скорость затухания разрыва выпуска

  taylor_forecasts <- numeric(forecast_horizon)
  current_inflation <- last_inflation
  current_gdp <- last_gdp_gap

  for (h in 1:forecast_horizon) {
    # Экспоненциальное затухание отклонений
    current_inflation <- current_inflation * exp(-inflation_decay * h)
    current_gdp <- current_gdp * exp(-gdp_decay * h)

    # Правило Тейлора с затуханием
    taylor_forecasts[h] <- TAYLOR_NEUTRAL_RATE +
                          TAYLOR_INFLATION_WEIGHT * current_inflation +
                          TAYLOR_OUTPUT_WEIGHT * current_gdp
  }

  # Комбинированный прогноз (взвешенное среднее методов)
  baseline_forecasts <- 0.3 * rep(mean_forecast, forecast_horizon) +
                       0.4 * linear_forecasts +
                       0.3 * taylor_forecasts

  # Ограничения на разумные значения
  baseline_forecasts <- pmax(pmin(baseline_forecasts, 25), 1)  # От 1% до 25%

  return(baseline_forecasts)
}

# Функция для создания красивого сравнительного графика моделей
CreateModelComparisonPlot <- function(baseline_forecasts, stochastic_results, train_data, forecast_dates) {

  # Подготовка исторических данных (последние 36 месяцев для лучшего контекста)
  historical_data <- train_data[!is.na(ruonia_rate)]
  last_36_months <- max(1, nrow(historical_data) - 35):nrow(historical_data)
  hist_subset <- historical_data[last_36_months]

  # Подготовка данных для сравнения
  comparison_data <- data.table(
    date = forecast_dates,
    baseline = baseline_forecasts,
    stochastic_median = stochastic_results$quantiles[, 4],  # Медиана
    stochastic_q25 = stochastic_results$quantiles[, 3],     # 25% квантиль
    stochastic_q75 = stochastic_results$quantiles[, 5],     # 75% квантиль
    stochastic_q05 = stochastic_results$quantiles[, 1],     # 5% квантиль
    stochastic_q95 = stochastic_results$quantiles[, 7]      # 95% квантиль
  )

  # Расчет статистик для аннотаций
  final_baseline <- last(baseline_forecasts)
  final_stochastic <- last(stochastic_results$quantiles[, 4])
  uncertainty_range <- last(stochastic_results$quantiles[, 5]) - last(stochastic_results$quantiles[, 3])

  # Создание красивого графика с белым фоном
  p <- plot_ly()

  # Стохастическая модель - внешний доверительный интервал (90%)
  p <- p %>% add_ribbons(
    data = comparison_data,
    x = ~date, ymin = ~stochastic_q05, ymax = ~stochastic_q95,
    fillcolor = 'rgba(255, 140, 0, 0.15)',
    line = list(color = 'rgba(255, 140, 0, 0.3)', width = 1),
    name = '90% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>90% ДИ: %{ymin:.1f}% - %{ymax:.1f}%<extra></extra>'
  )

  # Стохастическая модель - внутренний доверительный интервал (50%)
  p <- p %>% add_ribbons(
    data = comparison_data,
    x = ~date, ymin = ~stochastic_q25, ymax = ~stochastic_q75,
    fillcolor = 'rgba(255, 180, 0, 0.25)',
    line = list(color = 'rgba(255, 180, 0, 0.4)', width = 1),
    name = '50% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>50% ДИ: %{ymin:.1f}% - %{ymax:.1f}%<extra></extra>'
  )

    # Исторические данные с градиентной линией
  p <- p %>% add_lines(
    data = hist_subset,
    x = ~date, y = ~ruonia_rate,
    line = list(color = '#2E86AB', width = 3.5, shape = 'spline'),
    name = 'Исторические данные',
    hovertemplate = '<b>%{x}</b><br>Фактическая ставка: %{y:.2f}%<extra></extra>'
  )

  # Медиана стохастической модели
  p <- p %>% add_lines(
    data = comparison_data,
    x = ~date, y = ~stochastic_median,
    line = list(color = '#F7931E', width = 4, shape = 'spline'),
    name = 'Стохастическая модель (медиана)',
    hovertemplate = '<b>%{x}</b><br>Стохастический прогноз: %{y:.2f}%<extra></extra>'
  )

  # Базовая модель
  p <- p %>% add_lines(
    data = comparison_data,
    x = ~date, y = ~baseline,
    line = list(color = '#E74C3C', width = 3.5, dash = 'dash'),
    name = 'Базовая модель (детерминистическая)',
    hovertemplate = '<b>%{x}</b><br>Базовый прогноз: %{y:.2f}%<extra></extra>'
  )

    # Точка начала прогноза с красивым маркером
  p <- p %>% add_markers(
    x = last(hist_subset$date), y = last(hist_subset$ruonia_rate),
    marker = list(

      size = 20,
      symbol = 'diamond',
      line = list(color = '#E74C3C', width = 3)
    ),
    name = 'Начало прогноза',
    hovertemplate = '<b>Начало прогноза</b><br>%{x}<br>Ставка: %{y:.2f}%<extra></extra>'
  )

  # Аннотации с ключевыми метриками
  p <- p %>% add_annotations(
    x = 0.02, y = 0.98,
    xref = 'paper', yref = 'paper',
    text = paste0(
      '<b>Ключевые метрики:</b><br>',
      '• Базовая модель: ', round(final_baseline, 1), '%<br>',
      '• Стохастическая модель: ', round(final_stochastic, 1), '%<br>',
      '• Неопределенность: ±', round(uncertainty_range/2, 1), '%'
    ),
    showarrow = FALSE,
    font = list(color = '#2C3E50', size = 12),
    bgcolor = 'rgba(255,255,255,0.9)',
    bordercolor = 'rgba(44,62,80,0.3)',
    borderwidth = 1,
    align = 'left'
  )

  # Оформление с современным дизайном
  p <- p %>% layout(
    title = list(
      text = '<b>Сравнение моделей прогнозирования ставок ЦБ РФ</b><br><sub>Интегрированный подход: HAR + HMM + α-stable распределения + Taleb R²</sub>',
      font = list(color = '#2C3E50', size = 22, family = 'Arial Black'),
      x = 0.5,
      y = 0.98
    ),
    xaxis = list(
      title = list(
        text = '<b>Временной горизонт</b>',
        font = list(color = '#2C3E50', size = 16, family = 'Arial')
      ),
      color = '#2C3E50',
      gridcolor = 'rgba(200,200,200,0.5)',
      gridwidth = 1,
      tickfont = list(color = '#2C3E50', size = 12),
      zeroline = FALSE,
      showline = TRUE,
      linecolor = 'rgba(44,62,80,0.3)',
      linewidth = 1
    ),
    yaxis = list(
      title = list(
        text = '<b>Ключевая ставка ЦБ РФ, %</b>',
        font = list(color = '#2C3E50', size = 16, family = 'Arial')
      ),
      color = '#2C3E50',
      gridcolor = 'rgba(200,200,200,0.5)',
      gridwidth = 1,
      tickfont = list(color = '#2C3E50', size = 12),
      zeroline = FALSE,
      showline = TRUE,
      linecolor = 'rgba(44,62,80,0.3)',
      linewidth = 1,
      range = c(0, max(25, max(comparison_data$stochastic_q95) + 1))
    ),
    font = list(color = '#2C3E50', family = 'Arial'),
    margin = list(l = 100, r = 100, t = 150, b = 100),
    showlegend = TRUE
  )

  return(p)
}

# Создание базового прогноза для сравнения
baseline_forecasts <- CreateBaselineForecasts(train_data_with_regimes, HORIZON)

# Создание сравнительного графика
model_comparison_plot <- CreateModelComparisonPlot(
  baseline_forecasts,
  final_integrated_model,
  train_data_with_regimes,
  forecast_dates
)
model_comparison_plot

# Создание интерактивной сети марковских цепей
# markov_network <- CreateMarkovNetworkVisualization(hmm_model)
# markov_network

# Количественное сравнение моделей
forecast_comparison <- data.table(
  Horizon = seq_len(HORIZON),
  date = forecast_dates[seq_len(HORIZON)],
  Baseline = baseline_forecasts,
  Stochastic_Median = final_integrated_model$quantiles[, 4],
  Stochastic_Q05 = final_integrated_model$quantiles[, 1],
  Stochastic_Q25 = final_integrated_model$quantiles[, 3],
  Stochastic_Q75 = final_integrated_model$quantiles[, 5],
  Stochastic_Q95 = final_integrated_model$quantiles[, 7],
  Difference = final_integrated_model$quantiles[, 4] - baseline_forecasts,
  Uncertainty_Range = final_integrated_model$quantiles[, 5] - final_integrated_model$quantiles[, 3]
)

# Количественное сравнение моделей выполнено

print(forecast_comparison[, .(Horizon, Baseline, Stochastic_Median, Difference, Uncertainty_Range)])

#### 11.5. Визуализации марковских цепей и режимов ####

# Функция для создания интерактивной сети марковских цепей
CreateMarkovNetworkVisualization <- function(hmm_model, regime_names = NULL) {

  transition_matrix <- hmm_model$transition_matrix
  if (is.null(regime_names)) regime_names <- paste0("regime", 1:nrow(transition_matrix))
  N_REGIMES <- length(regime_names)

  # Цветовая палитра для режимов
  palette <- c("#3498DB", "#E74C3C", "#F39C12", "#2ECC71", "#9B59B6")

  # Узлы сети
  nodes <- data.frame(
    id = 1:N_REGIMES,
    label = regime_names,
    title = paste0(
      "<div style='font-family:Helvetica; font-size:20px;'>",
      "<b>", regime_names, "</b><br>",
      "Устойчивость: ", round(diag(transition_matrix), 3),
      "</div>"
    ),
    value = 1 + diag(transition_matrix) * 5,
    color = list(
      background = palette[1:N_REGIMES],
      border = "#2C3E50",
      highlight = list(background = "#F39C12", border = "#D35400")
    ),
    font = list(
      size = 18,
      face = "arial",
      color = "#2C3E50"
    ),
    shadow = TRUE
  )

  # Рёбра сети
  edges <- data.table()
  for (i in seq_len(N_REGIMES)) {
    for (j in seq_len(N_REGIMES)) {
      prob <- transition_matrix[i, j]
      if (prob > 0.01) {  # порог видимости
        edges <- rbind(edges, data.frame(
          from = i,
          to = j,
          label = sprintf("%.2f", prob),
          value = log(prob) * 2,
          title = paste0(
            "<div style='font-family:Helvetica; font-size:16px;'>",
            regime_names[i], " → ", regime_names[j], "<br>",
            "P = ", round(prob, 3),
            "</div>"
          ),
          color = if (i == j) "#95A5A6" else "#BDC3C7",
          arrows = "to"
        ))
      }
    }
  }

  # Построение интерактивной сети
  net <- visNetwork(nodes, edges) %>%
    visLayout(randomSeed = 2025) %>%
    visNodes(
      shape = "dot",
      scaling = list(min = 20, max = 60),
      borderWidthSelected = 5
    ) %>%
    visEdges(
      smooth = list(enabled = TRUE, type = "continuous"),
      width = htmlwidgets::JS("function(edge) { return Math.max(1, edge.value); }"),
      font = list(size = 16, align = "middle"),
      shadow = TRUE
    ) %>%
    visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)
    ) %>%
    visLegend(
      useGroups = FALSE,
      addNodes = data.frame(
        label = regime_names,
        shape = "dot",
        color = palette[1:N_REGIMES]
      ),
      width = 0.2,
      position = "right"
    ) %>%
    visPhysics(
      stabilization = list(enabled = TRUE, iterations = 250),
      barnesHut = list(
        gravitationalConstant = -30000,
        springLength = 120,
        springConstant = 0.002
      )
    ) %>%
    visInteraction(
      dragNodes = TRUE,
      dragView = TRUE,
      zoomView = TRUE,
      tooltipDelay = 150
    )

  return(net)
}

# Функция для создания тепловой карты матрицы переходов
CreateTransitionHeatmap <- function(transition_matrix, regime_names = NULL) {

  if (is.null(regime_names)) regime_names <- paste0("regime", 1:nrow(transition_matrix))

  # Подготовка данных для тепловой карты
  transition_data <- expand.grid(
    From = factor(regime_names, levels = regime_names),
    To = factor(regime_names, levels = rev(regime_names))
  )
  transition_data$Probability <- as.vector(t(transition_matrix[length(regime_names):1, ]))

  # Создание тепловой карты
  p <- plot_ly(
    data = transition_data,
    x = ~From, y = ~To, z = ~Probability,
    type = "heatmap",
    colorscale = list(
      c(0, "rgba(236,240,241,0.9)"),
      c(0.3, "rgba(52,152,219,0.9)"),
      c(0.6, "rgba(243,156,18,0.9)"),
      c(1, "rgba(231,76,60,0.95)")
    ),
    text = ~paste0(round(Probability * 100, 1), "%"),
    texttemplate = "%{text}",
    textfont = list(color = "#2C3E50", size = 16, family = "Arial Black"),
    hovertemplate = '<b>%{x} → %{y}</b><br>Вероятность: %{z:.3f}<br>%{text}<extra></extra>',
    showscale = TRUE,
    colorbar = list(
      title = list(text = "<b>Вероятность<br>перехода</b>", font = list(color = "#2C3E50", size = 14)),
      tickfont = list(color = "#2C3E50", size = 12),
      tickcolor = "#2C3E50",
      ticklen = 5
    )
  ) %>%
  layout(
    title = list(
      text = '<b>Матрица переходов между режимами монетарной политики ЦБ РФ</b><br><sub>HMM анализ скрытых состояний экономики</sub>',
      font = list(color = '#2C3E50', size = 18),
      x = 0.5
    ),
    xaxis = list(
      title = list(text = '<b>Исходный режим</b>', font = list(color = '#2C3E50', size = 16)),
      tickfont = list(color = '#2C3E50', size = 14),
      gridcolor = 'rgba(200,200,200,0.3)'
    ),
    yaxis = list(
      title = list(text = '<b>Целевой режим</b>', font = list(color = '#2C3E50', size = 16)),
      tickfont = list(color = '#2C3E50', size = 14),
      gridcolor = 'rgba(200,200,200,0.3)'
    ),
    plot_bgcolor = 'rgba(255,255,255,1)',
    paper_bgcolor = 'rgba(255,255,255,1)',
    font = list(color = '#2C3E50', family = 'Arial'),
    margin = list(l = 90, r = 140, t = 120, b = 90)
  )

  return(p)
}

# Функция для визуализации эволюции вероятностей режимов
CreateRegimeEvolutionPlot <- function(hmm_model, train_data, forecast_data) {

  # Исторические режимы
  historical_states <- ViterbiDecoding(train_data$regime_indicators, hmm_model)
  historical_probs <- data.table(
    date = train_data$date,
    regime = historical_states,
    type = "Historical"
  )

  # Прогнозные вероятности режимов
  forecast_probs <- data.table(
    date = forecast_data$date,
    normal_prob = forecast_data$regime_1,
    crisis_prob = forecast_data$regime_2,
    recovery_prob = forecast_data$regime_3,
    type = "Forecast"
  )

  # Создание графика исторических режимов
  p_historical <- plot_ly(
    data = historical_probs,
    x = ~date, y = ~regime,
    type = "scatter", mode = "markers+lines",
    marker = list(
      color = ~regime,
      colorscale = list(
        c(0, "#2E8B57"),    # Normal - SeaGreen
        c(0.5, "#FF6347"),  # Crisis - Tomato
        c(1, "#4169E1")     # Recovery - RoyalBlue
      ),
      size = 8,
      line = list(color = "white", width = 1)
    ),
    line = list(color = "rgba(255,255,255,0.3)", width = 2),
    name = "Исторические режимы",
    hovertemplate = '<b>%{x}</b><br>Режим: %{y}<br>Тип: Исторический<extra></extra>'
  )

  # Добавление прогнозных вероятностей как областей
  p_forecast <- plot_ly() %>%
    add_ribbons(
      data = forecast_probs,
      x = ~date, ymin = 0, ymax = ~normal_prob,
      fillcolor = 'rgba(46, 139, 87, 0.6)',
      line = list(color = 'transparent'),
      name = 'Вероятность нормального режима',
      hovertemplate = '<b>%{x}</b><br>Normal: %{ymax:.3f}<extra></extra>'
    ) %>%
    add_ribbons(
      data = forecast_probs,
      x = ~date, ymin = ~normal_prob, ymax = ~(normal_prob + crisis_prob),
      fillcolor = 'rgba(255, 99, 71, 0.6)',
      line = list(color = 'transparent'),
      name = 'Вероятность кризисного режима',
      hovertemplate = '<b>%{x}</b><br>Crisis: %{ymax:.3f}<extra></extra>'
    ) %>%
    add_ribbons(
      data = forecast_probs,
      x = ~date, ymin = ~(normal_prob + crisis_prob), ymax = 1,
      fillcolor = 'rgba(65, 105, 225, 0.6)',
      line = list(color = 'transparent'),
      name = 'Вероятность восстановительного режима',
      hovertemplate = '<b>%{x}</b><br>Recovery: %{ymax:.3f}<extra></extra>'
    )

  # Объединение графиков
  combined_plot <- subplot(
    p_historical, p_forecast,
    nrows = 2, heights = c(0.4, 0.6),
    subplot_titles = c(
      "<b>Исторические режимы (Витерби декодирование)</b>",
      "<b>Прогноз вероятностей режимов</b>"
    ),
    margin = 0.08,
    titleX = TRUE, titleY = TRUE
  ) %>%
  layout(
    title = list(
      text = '<b>Эволюция экономических режимов в модели ЦБ РФ</b><br><sub>Марковские цепи: переходы между нормальным, кризисным и восстановительным состояниями</sub>',
      font = list( size = 18),
      x = 0.5
    ),
    plot_bgcolor = 'rgba(15,15,25,0.95)',
    paper_bgcolor = 'rgba(15,15,25,0.95)',
    font = list( family = 'Arial'),
    showlegend = TRUE,
    legend = list(
      font = list( size = 12),
      bgcolor = 'rgba(0,0,0,0.8)',
      bordercolor = 'rgba(255,255,255,0.3)',
      borderwidth = 2
    ),
    margin = list(l = 90, r = 90, t = 140, b = 90)
  )

  return(combined_plot)
}

# 1. Тепловая карта матрицы переходов
transition_heatmap <- CreateTransitionHeatmap(hmm_model$transition_matrix, regime_names = paste0("regime", 1:nrow(hmm_model$transition_matrix)))
transition_heatmap

# 2. Эволюция режимов во времени
tryCatch({
  regime_evolution <- CreateRegimeEvolutionPlot(hmm_model, train_data_with_regimes, forecast_data)
  regime_evolution
}, error = function(e) {
  # Используется упрощенная визуализация режимов

  # Упрощенная альтернативная визуализация
  historical_states <- ViterbiDecoding(train_data_with_regimes$regime_indicators, hmm_model)

  # Проверяем совместимость размеров данных
  if (length(historical_states) == nrow(train_data_with_regimes)) {
    simple_regime_plot <- plot_ly(
      x = train_data_with_regimes$date,
      y = historical_states,
      type = "scatter", mode = "markers+lines",
      name = "Исторические режимы",
      marker = list(size = 8, color = historical_states, colorscale = "Viridis")
    ) %>%
    layout(
      title = "Эволюция экономических режимов ЦБ РФ",
      xaxis = list(title = "Дата"),
      yaxis = list(title = "Режим"),
      plot_bgcolor = 'rgba(15,15,25,0.95)',
      paper_bgcolor = 'rgba(15,15,25,0.95)',
      font = list( family = 'Arial')
    )

    simple_regime_plot
  } else {
    # Пропуск визуализации из-за несовместимых размеров данных
  }
})

# 3. Статистики устойчивости режимов
regime_stats <- data.table(
  Self_Transition = diag(hmm_model$transition_matrix),
  Mean_Level = hmm_model$emission_params$mean,
  Volatility = hmm_model$emission_params$sd,
  Expected_Duration = 1 / (1 - diag(hmm_model$transition_matrix))
)
print(regime_stats)

#### 11.8. Дополнительные визуализации ####

# Функция для создания графиков прогнозов макропоказателей
CreateMacroForecastPlots <- function(train_data, forecast_data, stochastic_results) {

  # Подготовка исторических данных (последние 24 месяца)
  historical_data <- train_data[!is.na(ruonia_rate)]
  last_24_months <- max(1, nrow(historical_data) - 23):nrow(historical_data)
  hist_subset <- historical_data[last_24_months]

  # 1. График инфляции
  inflation_plot <- plot_ly() %>%
    add_lines(
      data = hist_subset,
      x = ~date, y = ~cpi_inflation,
      line = list(color = '#FF6B6B', width = 3),
      name = 'Историческая инфляция',
      hovertemplate = '<b>%{x}</b><br>Инфляция: %{y:.2f}%<extra></extra>'
    ) %>%
    add_lines(
      data = forecast_data,
      x = ~date, y = ~(cpi_inflation),
      line = list(color = '#FFD700', width = 3, dash = 'dash'),
      name = 'Прогноз инфляции (SVD-VAR)',
      hovertemplate = '<b>%{x}</b><br>Прогноз: %{y:.2f}%<extra></extra>'
    ) %>%
    layout(
      title = list(text = '<b>Прогноз инфляции ЦБ РФ</b><br><sub>SVD-VAR экстраполяция + цель ЦБ</sub>',
                   font = list(size = 18)),
      xaxis = list(title = 'Дата'),
      yaxis = list(title = 'Инфляция (%)'),
      legend = list(orientation = "l", x = 0.1, y = 0.1)
    )

  # 2. График разрыва ВВП
  gdp_plot <- plot_ly() %>%
    add_lines(
      data = hist_subset,
      x = ~date, y = ~gdp_gap,
      line = list(color = '#4ECDC4', width = 3),
      name = 'Исторический разрыв ВВП',
      hovertemplate = '<b>%{x}</b><br>Разрыв ВВП: %{y:.2f}%<extra></extra>'
    ) %>%
    add_lines(
      data = forecast_data,
      x = ~date, y = ~gdp_gap,
      line = list(color = '#FFD700', width = 3, dash = 'dash'),
      name = 'Прогноз разрыва ВВП (SVD-VAR)',
      hovertemplate = '<b>%{x}</b><br>Прогноз: %{y:.2f}%<extra></extra>'
    ) %>%
    layout(
      title = list(text = '<b>Прогноз разрыва выпуска</b><br><sub>Отклонение от потенциального ВВП</sub>',
                   font = list( size = 18)),
      xaxis = list(title = 'Дата'),
      yaxis = list(title = 'Разрыв ВВП (%)'),
      legend = list(orientation = "l", x = 0.5, y = 0.9)
    )

  # 3. График волатильности (HAR vs фактическая)
  volatility_plot <- plot_ly() %>%
    add_lines(
      data = hist_subset,
      x = ~date, y = ~rate_volatility,
      line = list(color = '#96CEB4', width = 3),
      name = 'Историческая волатильность',
      hovertemplate = '<b>%{x}</b><br>Волатильность: %{y:.4f}<extra></extra>'
    ) %>%
    add_lines(
      data = forecast_data,
      x = ~date, y = ~har_volatility,
      line = list(color = '#FFD700', width = 3, dash = 'dash'),
      name = 'HAR прогноз волатильности',
      hovertemplate = '<b>%{x}</b><br>HAR прогноз: %{y:.4f}<extra></extra>'
    ) %>%
    # Добавление стохастических путей волатильности
    add_ribbons(
      x = forecast_data$date,
      ymin = apply(stochastic_results$volatility_simulations, 1, quantile, 0.1),
      ymax = apply(stochastic_results$volatility_simulations, 1, quantile, 0.9),
      line = list(color = 'transparent'),
      name = '80% ДИ стохастической волатильности',
      hovertemplate = '<b>%{x}</b><br>80% ДИ: %{ymin:.4f} - %{ymax:.4f}<extra></extra>'
    ) %>%
    layout(
      title = list(text = '<b>Прогноз волатильности процентных ставок</b><br><sub>HAR модель + стохастическая эволюция</sub>',
                   font = list( size = 18)),
      xaxis = list(title = 'Дата'),
      yaxis = list(title = 'Волатильность'),
      legend = list(orientation = "l", x = 0.5, y = 0.9)
    )

  # Комбинированный график
  combined_macro <- subplot(
    inflation_plot, gdp_plot, volatility_plot,
    nrows = 3, margin = 0.08,
    titleX = TRUE, titleY = TRUE
  ) %>%
  layout(
    title = list(
      text = '<b>Прогнозы макроэкономических индикаторов ЦБ РФ</b><br><sub>SVD-VAR экстраполяция + HAR волатильность + стохастические компоненты</sub>',
      font = list( size = 20),
      x = 0.5
    ),
    margin = list(l = 90, r = 90, t = 140, b = 90)
  )

  return(list(
    inflation = inflation_plot,
    gdp = gdp_plot,
    volatility = volatility_plot,
    combined = combined_macro
  ))
}

# Функция для создания продвинутой визуализации марковских цепей
CreateAdvancedMarkovChainPlot <- function(hmm_model, train_data, forecast_data) {

  # 1. Анимированная тепловая карта матрицы переходов (в стиле 5.HMM)
  transition_matrix <- hmm_model$transition_matrix
  regime_names <- paste0("regime", 1:nrow(transition_matrix))

  # Подготовка данных для тепловой карты
  transition_data <- expand.grid(
    From = factor(regime_names, levels = regime_names),
    To = factor(regime_names, levels = rev(regime_names))
  )
  transition_data$Probability <- as.vector(t(transition_matrix[length(regime_names):1, ]))

  # Создание продвинутой тепловой карты
  advanced_heatmap <- plot_ly(
    data = transition_data,
    x = ~From, y = ~To, z = ~Probability,
    type = "heatmap",
    colorscale = list(
      c(0, "rgba(200,200,200,0.95)"),
      c(0.3, "rgba(255,244,210,0.95)"),
      c(0.6, "rgba(255,234,130,0.95)"),
      c(1, "rgba(255,220,50,0.95)")
    ),
    text = ~paste0(round(Probability * 100, 1), "%"),
    texttemplate = "%{text}",
    hovertemplate = '<b>%{x} → %{y}</b><br>Вероятность: %{z:.3f}<br>%{text}<extra></extra>',
    showscale = TRUE,
    colorbar = list(
      title = list(text = "<b>Вероятность<br>перехода</b>"),
      ticklen = 8,
      thickness = 20
    )
  ) %>%
  layout(
    title = list(
      text = '<b>Матрица переходов марковских режимов монетарной политики ЦБ РФ</b><br><sub>Вероятности смены экономических режимов (HMM анализ)</sub>',
      x = 0.5
    ),
    xaxis = list(
      title = list(text = '<b>Исходный режим</b>', font = list( size = 18)),
      tickfont = list( size = 16)
    ),
    yaxis = list(
      title = list(text = '<b>Целевой режим</b>', font = list( size = 18)),
      tickfont = list( size = 16)

    ),
    margin = list(l = 120, r = 160, t = 140, b = 90)
  )

  # 2. График эволюции режимных вероятностей (стек по всем режимам)
  # Собираем матрицу вероятностей режимов regime_1..regime_N
  regime_cols <- paste0("regime_", 1:length(regime_names))
  available_cols <- regime_cols[regime_cols %in% names(forecast_data)]
  N <- length(available_cols)

  regime_evolution <- plot_ly()
  if (N > 0) {
    prob_mat <- as.matrix(forecast_data[, ..available_cols])
    cum_mat <- t(apply(prob_mat, 1, cumsum))

    base_palette <- c('#2E8B57', '#FF6347', '#4169E1', '#F1C40F', '#8E44AD', '#1ABC9C', '#E67E22')
    for (k in seq_len(N)) {
      y0 <- if (k == 1) rep(0, nrow(forecast_data)) else cum_mat[, k - 1]
      y1 <- cum_mat[, k]
      regime_evolution <- regime_evolution %>% add_ribbons(
        data = forecast_data,
        x = ~date, ymin = y0, ymax = y1,
        line = list(color = 'transparent'),
        fillcolor = base_palette[(k - 1) %% length(base_palette) + 1],
        name = available_cols[k],
        hovertemplate = '<b>%{x}</b><br>Накопл. вероятность: %{ymax:.3f}<extra></extra>'
      )
    }
    regime_evolution <- regime_evolution %>% layout(
      title = list(
        text = '<b>Прогнозная эволюция режимных вероятностей</b><br><sub>Марковские цепи: динамика экономических состояний</sub>',
        font = list( size = 18),
        x = 0.5
      ),
      xaxis = list(title = 'Дата'),
      yaxis = list(title = 'Накопленная вероятность'),
      legend = list(x = 0.02, y = 0.98),
      margin = list(l = 90, r = 90, t = 120, b = 90)
    )
  }

  # 3. Комбинированная визуализация
  combined_markov <- subplot(
    advanced_heatmap, regime_evolution,
    nrows = 2, heights = c(0.6, 0.4),
    margin = 0.08,
    titleX = TRUE, titleY = TRUE
  ) %>%
  layout(
    title = list(
      text = '<b>Анализ марковских режимов денежно-кредитной политики ЦБ РФ</b><br><sub>Переходы между состояниями экономики и прогнозная динамика</sub>',
      font = list( size = 22),
      x = 0.5
    ),
    margin = list(l = 100, r = 100, t = 160, b = 100)
  )

  return(list(
    heatmap = advanced_heatmap,
    evolution = regime_evolution,
    combined = combined_markov
  ))
}

# Создание дополнительных визуализаций

# Макропоказатели
macro_plots <- NULL
tryCatch({
  macro_plots <- CreateMacroForecastPlots(train_data_with_regimes, forecast_data, final_integrated_model)
  macro_plots$combined
}, error = function(e) {
  cat("Визуализация макропоказателей пропущена из-за ошибки:", e$message, "\n")
})

# Продвинутая марковская визуализация
advanced_markov_plots <- CreateAdvancedMarkovChainPlot(hmm_model, train_data_with_regimes, forecast_data)
advanced_markov_plots$combined

# Подготовка экспортируемых режимных данных
regime_sequence_export <- train_data_with_regimes[, .(date, regime = as.integer(as.character(regime)))]
regime_prob_cols <- grep("^regime_", names(forecast_data), value = TRUE)
forecast_regime_probs_export <- if (length(regime_prob_cols)) {
  forecast_data[, c("date", regime_prob_cols), with = FALSE]
} else {
  NULL
}

#### 11.9. Обновленная документация модели ####

# Финальная документация и архитектура модели доступна в final_results
hmm_model$smoothed_states <- hidden_states
# Сохранение всех результатов включая новые визуализации
final_results <- list(
  train_data = train_data_with_regimes,
  forecast_data = forecast_data,
  stochastic_parameters = stochastic_params,
  har_forecast = forecast_vals,
  hmm_model = hmm_model,
  taylor_models = taylor_models,
  taylor_comparison = taylor_comparison,
  baseline_forecasts = baseline_forecasts,
  forecast_comparison = forecast_comparison,
  regime_stats = regime_stats,
  shapley_analysis = shapley_r2,
  integrated_model = final_integrated_model,
  diagnostics = diagnostics,
  tail_risk = tail_risk,
  forecast_evaluation = forecast_evaluation,
  volatility_analysis = volatility_analysis,
  model_formulas = final_integrated_model$model_formulas,
  calibrated_shocks_info = list(
    approach = "Калибровка через оптимизацию МНП для каждой переменной",
    distributions = stochastic_params$calibration_summary,
    description = "Автоматический выбор оптимального жирохвостого распределения (Pareto/t-Student/α-stable/Normal) для каждой макропеременной"
  ),

  # Новые компоненты
  macro_visualizations = macro_plots,
  markov_visualizations = advanced_markov_plots,
  improved_monte_carlo = list(
    enhanced_paths = TRUE,
    extreme_scenarios = TRUE,
    colorful_visualization = TRUE,
    fixed_date_gap = TRUE
  ),
  regime_sequence = regime_sequence_export,
  forecast_regime_probabilities = forecast_regime_probs_export,
  performance_profile = performance_profile,
  regularisation_log = regularisation_log,
  validation_report = validation_report,

  timestamp = Sys.time()
)

SaveProjectResults(final_results, "3.Data/final_integrated_cb_model.qs")
