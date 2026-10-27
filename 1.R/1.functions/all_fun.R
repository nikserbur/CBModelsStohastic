# all_fun.R — Объединенный файл всех функций стохастического моделирования ставок ЦБ
# Проект: Стохастическое моделирование процентных ставок с жирными хвостами
# Интеграция: HAR модели, HMM, Gradient Boosting, Taleb R² методики
# Автор: Бураков Никита Сергеевич

#### Подключение библиотек ####
source("1.R/9.libraries.r")

#### Константы проекта ####
DAILY_WINDOW <- 1
WEEKLY_WINDOW <- 5
MONTHLY_WINDOW <- 22
QUARTERLY_WINDOW <- 66
ALPHA_TAIL <- 3.5
N_BOOTSTRAP <- 1000
TOLERANCE <- 1e-6

#### Rcpp интеграция ####

# Проверка и подключение Rcpp оптимизаций
RCPP_AVAILABLE <- FALSE
tryCatch({
  if (file.exists("2.Cpp/FastMonteCarlo.cpp")) {
    library(Rcpp)
    sourceCpp("2.Cpp/FastMonteCarlo.cpp")
    RCPP_AVAILABLE <- TRUE
    cat("Rcpp оптимизации успешно подключены\n")
  }
}, error = function(e) {
  cat("Rcpp не доступен, используется R реализация:", e$message, "\n")
  RCPP_AVAILABLE <- FALSE
})

#### Загрузка и подготовка данных ####

#' prep_data — загрузка и объединение месячной панели из 3.Data (используется в Main.R
#' и DataAnalysisAndRegression.R); пропуски заполняются линейной регрессией по времени
prep_data <- function() {
# исходные данные тут https://www.cbr.ru/dkp/statistics/

safe_read_xlsx <- function(path, ...) {
  tmp <- tempfile(fileext = ".xlsx")
  file.copy(path, tmp, overwrite = TRUE)
  read_xlsx(tmp, ...)
}

# Данные и инфляции
inflation <- safe_read_xlsx("3.Data/Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx") %>%
  as.data.table()
setnames(inflation, names(inflation)[1], "date")
inflation <- inflation %>%
  .[, date := as.Date(paste0("01.", date), format = "%d.%m.%Y")] %>%
  na.omit() %>%
  setnames(names(.)[2], "key_rate") %>%
  setnames(names(.)[3], "cpi_inflation") %>%
  .[,-4]

# Бизнес климат РФ
bisness_data0 <- safe_read_xlsx("3.Data/add_info.xlsx", sheet = 2) %>%
  as.data.table()
setnames(bisness_data0, names(bisness_data0)[1], "date")


spr_bisness_data <- data.table(ind = names(bisness_data0)[which(names(bisness_data0) != "date")]) %>%
  .[,ind_code := paste0("B00",1:nrow(.))]


bisness_data <- copy(bisness_data0) %>%
  setnames(names(bisness_data0)[which(names(bisness_data0) != "date")],spr_bisness_data$ind_code) %>%
  .[,date := as.Date(date)]


# Платежный баланс
bal_of_pay0 <- safe_read_xlsx("3.Data/bal_of_payments_ap.xlsx") %>%
  as.data.table() %>%
  melt(id.vars = "ind")

# Создание справочника показателей платежного баланса
spr_bal_of_pay <- data.table(ind = unique(bal_of_pay0$ind)) %>%
  .[,ind_code := paste0("BP00",1:nrow(.))]

# Функция для преобразования квартальных дат в месячные
convert_quarter_to_monthly <- function(quarter_str) {
  quarter_match <- stringr::str_match(quarter_str, "(\\d+)\\s+квартал\\s+(\\d+)\\s+г\\.")
  if (is.na(quarter_match[1])) return(NA)

  quarter <- as.numeric(quarter_match[2])
  year <- as.numeric(quarter_match[3])

  first_month <- (quarter - 1) * 3 + 3
  return(as.Date(paste0(year, "-", sprintf("%02d", first_month), "-01")))
}

# Обработка данных платежного баланса с преобразованием квартальных дат
bal_of_pay <- copy(bal_of_pay0)
bal_of_pay[, date := convert_quarter_to_monthly(variable), by = "variable"]
bal_of_pay <- merge(bal_of_pay, spr_bal_of_pay, by = "ind", all.x = TRUE)
bal_of_pay[, ind := ind_code]
bal_of_pay <- bal_of_pay[, -c("ind", "variable")]
bal_of_pay <- bal_of_pay[!is.na(date)]
setorder(bal_of_pay, date)
bal_of_pay <- dcast(bal_of_pay, date ~ ind_code)

# Создание полного месячного ряда дат
min_date <- min(bal_of_pay$date)
max_date <- max(bal_of_pay$date)
full_dates <- data.table(date = seq.Date(min_date, max_date, by = "month"))

# Объединение с полным рядом дат и заполнение пропусков
bal_of_pay_full <- merge(full_dates, bal_of_pay, by = "date", all.x = TRUE) %>%
  # Заполнение пропусков для квартальных данных (интерполяция)
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

# Кредитный импульс
credit_impl <- fread("3.Data/кредитный_импульс_помесячно.csv") %>%
  as.data.table() %>%
  .[, date := as.Date(paste0(`Дата`,"-01"), format = "%Y-%m-%d")] %>%
  .[, `Дата` := NULL] %>%
  setnames("Кредитный импульс", "credit_impulse")

# Межбанковская реальная ставка
st_ruonia <- safe_read_xlsx("3.Data/RC_F11_01_2010_T04_07_2025.xlsx") %>%
  as.data.table() %>%
  .[,c("DT","ruo")] %>%
  .[, ym := format(as.Date(DT), "%Y-%m")] %>%
  .[, .(`ruonia_rate` = mean(ruo, na.rm = TRUE)), by = .(ym)] %>%
  .[,date := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")]

#ВВП
gdp <- fread("3.Data/vvp_pomesyachno.csv") %>%
  as.data.table() %>%
  na.omit() %>%
  setnames("ВВП_месяц","gdp_level") %>%
  setnames("Дата", "date") %>%
  .[, date := as.Date(date)]

#Бюджетный импулс
budget_impl <- safe_read_xlsx("3.Data/budget.xlsx") %>%
  as.data.table() %>%
  setnames("Date","date") %>%
  .[,date := as.Date(paste0(date ,"-01"), format = "%Y-%m-%d")] %>%
  setnames("Value","budget_impulse") %>%
  .[,budget_impulse := as.numeric(budget_impulse)]

#Скор санкций
sunctions <- fread("3.Data/sunctions.csv") %>%
  setnames(c("date","sunctions_score")) %>%
  .[, ym := format(as.Date(date), "%Y-%m")] %>%
  .[, .(`sunctions_score` = mean(`sunctions_score`, na.rm = TRUE)), by = .(ym)] %>%
  .[,date := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")]

respectation_inf <- fread("3.Data/ожидаемая_инфляция.csv") %>%
  na.omit() %>%
  .[,date := as.Date(paste0(`Дата`,"-01"), format = "%Y-%m-%d")] %>%
  .[, `Дата` := NULL] %>%
  setnames("Ожидаемая инфляция (%)", "expected_inflation")

#Валатильность активов
wind <- 30
all_input_data <- safe_read_xlsx("3.Data/assets.xlsx") %>%
  as.data.table()
setnames(all_input_data, names(all_input_data)[1], "date")
all_input_data <- all_input_data %>%
  na.omit() %>%
  .[order(date)] %>%
  .[, all := MOEX + RTS] %>%
  .[, vol := frollapply(all, wind, function(all){
    mean <- mean(all)
    res <- sum(abs(mean - all))/length(all)
    })*sqrt(wind)] %>%
  .[,c("date","vol")] %>%
  na.omit() %>%
  .[, ym := format(as.Date(date), "%Y-%m")] %>%
  .[, .(volatility = mean(vol, na.rm = TRUE)), by = .(ym)] %>%
  .[,date := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")] %>%
  merge(inflation, by = "date", all.x = TRUE) %>%
  merge(bal_of_pay_full, by = "date", all.x = TRUE) %>%
  merge(bisness_data, by = "date", all.x = TRUE) %>%
  merge(credit_impl, by = "date", all.x = TRUE) %>%
  merge(st_ruonia, by = "date", all.x = TRUE) %>%
  merge(gdp, by = "date", all.x = TRUE) %>%
  merge(budget_impl, by = "date", all.x = TRUE) %>%
  # merge(sunctions, by = "date", all.x = TRUE) %>%
  merge(respectation_inf, by = "date", all.x = TRUE) %>%
  # .[1,`Скор_санкций`:=2037] %>%
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

  return(all_input_data)
}

#### 0. Прогнозирование временных рядов через SVD с сезонностью и колебательностью ####

#' PredictTimeSeriesWithSVD — Прогнозирование через сингулярное разложение с учетом сезонности и циклических паттернов
#' @param input_dataset data.table с временными рядами
#' @param prediction_steps горизонт прогнозирования
#' @param temporal_var название колонки времени
#' @param seasonal_period период сезонности (по умолчанию 12 для месячных данных)
#' @param enable_seasonality включить сезонную декомпозицию
#' @param enable_cyclical включить циклические компоненты
#' @return data.table с прогнозными значениями
PredictTimeSeriesWithSVD <- function(input_dataset, prediction_steps = 12, temporal_var = "date",
                                    seasonal_period = 12, enable_seasonality = TRUE, enable_cyclical = TRUE) {



  # Подготовка временных рядов
  ts_matrix <- copy(input_dataset)

  # Извлечение числовых характеристик (исключая временной индекс)
  feature_names <- setdiff(names(ts_matrix)[sapply(ts_matrix, is.numeric)], temporal_var)

  if (length(feature_names) < 2) {
    stop("Критически мало признаков для SVD декомпозиции")
  }

  cat("Количество признаков для анализа:", length(feature_names), "\n")
  cat("Длина временного ряда:", nrow(ts_matrix), "периодов\n")
  cat("Период сезонности:", seasonal_period, "\n")

  # 1. Сезонная декомпозиция для каждого признака
  seasonal_components <- list()
  trend_components <- list()
  residual_components <- list()

  if (enable_seasonality && nrow(ts_matrix) >= seasonal_period * 2) {
    cat("Применение сезонной декомпозиции...\n")

    for (feature in feature_names) {
      series_data <- ts_matrix[[feature]]

      # Проверка на достаточность данных для сезонной декомпозиции
      if (length(series_data) >= seasonal_period * 2 && !all(is.na(series_data))) {
        tryCatch({
          # Создание временного ряда с правильной частотой
          ts_object <- ts(series_data, frequency = seasonal_period)

          # Сезонная декомпозиция через STL
          stl_decomp <- stl(ts_object, s.window = "periodic", robust = TRUE)

          seasonal_components[[feature]] <- as.numeric(stl_decomp$time.series[, "seasonal"])
          trend_components[[feature]] <- as.numeric(stl_decomp$time.series[, "trend"])
          residual_components[[feature]] <- as.numeric(stl_decomp$time.series[, "remainder"])

        }, error = function(e) {
          # Fallback: простая сезонная декомпозиция
          seasonal_components[[feature]] <<- rep(0, length(series_data))
          trend_components[[feature]] <<- series_data
          residual_components[[feature]] <<- rep(0, length(series_data))
        })
      } else {
        seasonal_components[[feature]] <- rep(0, length(series_data))
        trend_components[[feature]] <- series_data
        residual_components[[feature]] <- rep(0, length(series_data))
      }
    }
  } else {
    # Без сезонной декомпозиции
    for (feature in feature_names) {
      seasonal_components[[feature]] <- rep(0, nrow(ts_matrix))
      trend_components[[feature]] <- ts_matrix[[feature]]
      residual_components[[feature]] <- rep(0, nrow(ts_matrix))
    }
  }

  # 2. Создание расширенной матрицы признаков с сезонными и циклическими компонентами
  enhanced_features <- list()

  for (feature in feature_names) {
    base_series <- trend_components[[feature]]

    # Базовый признак
    enhanced_features[[paste0(feature, "_trend")]] <- base_series

    # Сезонный компонент
    if (enable_seasonality) {
      enhanced_features[[paste0(feature, "_seasonal")]] <- seasonal_components[[feature]]
    }

    # Циклические компоненты (синусы и косинусы)
    if (enable_cyclical && nrow(ts_matrix) >= seasonal_period) {
      time_index <- 1:nrow(ts_matrix)

      # Основная сезонная частота
      enhanced_features[[paste0(feature, "_sin_annual")]] <- sin(2 * pi * time_index / seasonal_period)
      enhanced_features[[paste0(feature, "_cos_annual")]] <- cos(2 * pi * time_index / seasonal_period)

      # Полугодовая частота
      enhanced_features[[paste0(feature, "_sin_semi")]] <- sin(2 * pi * time_index / (seasonal_period / 2))
      enhanced_features[[paste0(feature, "_cos_semi")]] <- cos(2 * pi * time_index / (seasonal_period / 2))

      # Квартальная частота
      enhanced_features[[paste0(feature, "_sin_quarter")]] <- sin(2 * pi * time_index / (seasonal_period / 4))
      enhanced_features[[paste0(feature, "_cos_quarter")]] <- cos(2 * pi * time_index / (seasonal_period / 4))
    }
  }

  # Создание расширенной матрицы данных
  enhanced_matrix <- as.data.table(enhanced_features)

  # 3. Z-нормализация расширенной матрицы
  central_moments <- apply(enhanced_matrix, 2, mean, na.rm = TRUE)
  scale_factors <- apply(enhanced_matrix, 2, sd, na.rm = TRUE)
  scale_factors[scale_factors == 0 | is.na(scale_factors)] <- 1

  normalized_matrix <- as.data.table(scale(enhanced_matrix,
                                          center = central_moments,
                                          scale = scale_factors))

  # 4. Сингулярное разложение (SVD) расширенной матрицы
  svd_decomposition <- svd(as.matrix(normalized_matrix))

  # Определение оптимального ранга (энергетический критерий 95%)
  singular_energy <- cumsum(svd_decomposition$d^2) / sum(svd_decomposition$d^2)
  optimal_rank <- min(which(singular_energy >= 0.95), length(svd_decomposition$d))

  cat("Выбранный ранг SVD:", optimal_rank, "\n")
  cat("Покрытая энергия сигнала:", round(singular_energy[optimal_rank] * 100, 1), "%\n")

  # 5. Проецирование данных в пространство главных направлений
  left_singular_vectors <- svd_decomposition$u[, 1:optimal_rank]
  projected_data <- as.data.table(left_singular_vectors)
  setnames(projected_data, paste0("svd_", 1:optimal_rank))

  temporal_series <- ts_matrix[[temporal_var]]
  dynamics_table <- cbind(
    data.table(setNames(list(temporal_series), temporal_var)),
    projected_data
  )

  # 6. Генерация будущих временных точек
  if (class(temporal_series)[1] == "Date") {
    future_periods <- seq(max(temporal_series) + months(1),
                         max(temporal_series) + months(prediction_steps),
                         by = "month")
  } else {
    future_periods <- seq(max(temporal_series) + 1,
                         max(temporal_series) + prediction_steps)
  }

  # 7. Прогнозирование каждой SVD компоненты с учетом сезонности
  n_components <- ncol(projected_data)
  forecasted_components <- matrix(0, nrow = length(temporal_series) + prediction_steps, ncol = n_components)

  # Заполняем исторические данные
  forecasted_components[1:length(temporal_series), ] <- as.matrix(projected_data)

  # Прогнозирование каждой компоненты
  for (comp_idx in 1:n_components) {
    historical_values <- projected_data[[comp_idx]]

    # Применяем улучшенную экстраполяцию с сезонностью
    full_series <- c(historical_values, rep(NA, prediction_steps))
    temporal_indices <- seq_along(full_series)

    extrapolated <- ExtrapolateLatentFactorWithSeasonality(
      temporal_index = temporal_indices,
      observed_values = full_series,
      final_horizon = length(full_series),
      seasonal_period = seasonal_period,
      enable_seasonality = enable_seasonality
    )$projected_series

    forecasted_components[, comp_idx] <- extrapolated
  }

  # 8. Восстановление через правые сингулярные векторы
  right_basis <- svd_decomposition$v[, 1:optimal_rank]
  singular_values <- diag(svd_decomposition$d[1:optimal_rank])

  # Реконструкция в расширенном пространстве признаков
  reconstructed_enhanced <- forecasted_components %*% singular_values %*% t(right_basis)

  # 9. Денормализация расширенных признаков
  final_enhanced <- as.data.table(sapply(1:ncol(reconstructed_enhanced), function(idx) {
    reconstructed_enhanced[, idx] * scale_factors[idx] + central_moments[idx]
  }))

  setnames(final_enhanced, names(enhanced_features))

  # 10. Восстановление исходных признаков из расширенных
  final_reconstruction <- as.data.table(matrix(0, nrow = nrow(final_enhanced), ncol = length(feature_names)))
  setnames(final_reconstruction, feature_names)

  for (i in 1:length(feature_names)) {
    feature <- feature_names[i]

    # Собираем компоненты
    trend_component <- final_enhanced[[paste0(feature, "_trend")]]

    if (enable_seasonality) {
      seasonal_component <- final_enhanced[[paste0(feature, "_seasonal")]]
    } else {
      seasonal_component <- rep(0, nrow(final_enhanced))
    }

    # Суммируем компоненты
    final_reconstruction[[feature]] <- trend_component + seasonal_component
  }

  # 11. Объединение временного индекса с предсказаниями
  full_temporal_sequence <- c(temporal_series, future_periods)
  comprehensive_result <- cbind(
    data.table(setNames(list(full_temporal_sequence), temporal_var)),
    final_reconstruction
  )

  # Разделение на ретроспективу и прогноз
  historical_length <- length(temporal_series)

  output_structure <- list(
    historical = comprehensive_result[1:historical_length],
    forecast = comprehensive_result[(historical_length + 1):nrow(comprehensive_result)],
    decomposition_summary = list(
      energy_captured = singular_energy[optimal_rank],
      svd_rank = optimal_rank,
      total_dimensions = length(svd_decomposition$d),
      seasonal_period = seasonal_period,
      enhanced_features = length(enhanced_features)
    ),
    transformation_params = list(
      centralization = central_moments,
      normalization = scale_factors,
      basis_matrix = right_basis,
      seasonal_components = seasonal_components,
      trend_components = trend_components
    )
  )

  cat("SVD-VAR прогнозирование с сезонностью выполнено успешно\n")

  return(output_structure)
}

#' ExtrapolateLatentFactorWithSeasonality — Экстраполяция скрытых факторов с учетом сезонности
#' @param temporal_index временные индексы
#' @param observed_values наблюдаемые значения
#' @param final_horizon конечный горизонт
#' @param seasonal_period период сезонности
#' @param enable_seasonality включить сезонную экстраполяцию
#' @return список с экстраполированными значениями
ExtrapolateLatentFactorWithSeasonality <- function(temporal_index, observed_values, final_horizon,
                                                  seasonal_period = 12, enable_seasonality = TRUE) {

  # Фильтрация корректных наблюдений
  valid_mask <- !is.na(observed_values) & !is.na(temporal_index)
  clean_temporal <- temporal_index[valid_mask]
  clean_observations <- observed_values[valid_mask]

  if (length(clean_observations) < 3) {
    # Критически мало данных - константная экстраполяция
    extrapolated_values <- rep(ifelse(length(clean_observations) > 0, last(clean_observations), 0),
                              final_horizon - max(clean_temporal, na.rm = TRUE))
    projected_series <- c(clean_observations, extrapolated_values)
    return(list(projected_series = projected_series[temporal_index]))
  }

  # Определение индексов для экстраполяции
  missing_indices <- which(is.na(observed_values))
  extrapolation_length <- length(missing_indices)

  if (extrapolation_length == 0) {
    return(list(projected_series = observed_values))
  }

  tryCatch({
    # Метод 1: Сезонная экстраполяция с VAR
    if (enable_seasonality && length(clean_observations) >= seasonal_period * 2) {

      # Создание временного ряда с сезонностью
      ts_object <- ts(clean_observations, frequency = seasonal_period)

      # Сезонная декомпозиция
      stl_decomp <- stl(ts_object, s.window = "periodic", robust = TRUE)

      # Прогнозирование тренда
      trend_series <- as.numeric(stl_decomp$time.series[, "trend"])
      seasonal_series <- as.numeric(stl_decomp$time.series[, "seasonal"])

      # VAR для тренда
      if (length(trend_series) >= 4) {
        trend_ts <- ts(trend_series, frequency = 1)
        max_lag <- min(3, floor(length(trend_series) / 3))

        if (max_lag >= 1) {
          lag_selection <- vars::VARselect(trend_ts, lag.max = max_lag, type = "const")
          optimal_lag <- lag_selection$selection[["AIC(n)"]]

          var_model <- vars::VAR(trend_ts, p = optimal_lag, type = "const")
          trend_prediction <- predict(var_model, n.ahead = extrapolation_length)
          forecasted_trend <- as.numeric(trend_prediction$fcst[[1]][, "fcst"])
        } else {
          forecasted_trend <- rep(mean(trend_series), extrapolation_length)
        }
      } else {
        forecasted_trend <- rep(mean(trend_series), extrapolation_length)
      }

      # Экстраполяция сезонного компонента
      if (extrapolation_length <= seasonal_period) {
        # Простая экстраполяция сезонности
        seasonal_forecast <- rep(seasonal_series[1:extrapolation_length],
                               ceiling(extrapolation_length / seasonal_period))[1:extrapolation_length]
      } else {
        # Циклическая экстраполяция сезонности
        seasonal_forecast <- rep(seasonal_series, ceiling(extrapolation_length / seasonal_period))[1:extrapolation_length]
      }

      # Комбинирование тренда и сезонности
      extrapolated_values <- forecasted_trend + seasonal_forecast

    } else {
      # Метод 2: Стандартная VAR без сезонности
      timeseries_object <- ts(clean_observations, frequency = 1)
      max_lag <- min(4, floor(length(clean_observations) / 3))

      if (max_lag >= 1) {
        lag_selection <- vars::VARselect(timeseries_object, lag.max = max_lag, type = "const")
        optimal_lag <- lag_selection$selection[["AIC(n)"]]

        var_model <- vars::VAR(timeseries_object, p = optimal_lag, type = "const")
        var_prediction <- predict(var_model, n.ahead = extrapolation_length)
        extrapolated_values <- as.numeric(var_prediction$fcst[[1]][, "fcst"])
      } else {
        extrapolated_values <- rep(mean(clean_observations), extrapolation_length)
      }
    }

    # Объединение наблюдений и прогнозов
    reconstructed_series <- observed_values
    reconstructed_series[missing_indices] <- extrapolated_values

    # Проверка на экстремальные значения и их ограничение
    if (any(abs(extrapolated_values) > 10 * sd(clean_observations, na.rm = TRUE))) {
      # Если прогнозы слишком экстремальные, используем экспоненциальное затухание
      last_value <- tail(clean_observations, 1)
      decay_factor <- 0.95
      for (i in 1:length(extrapolated_values)) {
        extrapolated_values[i] <- last_value * (decay_factor^i)
      }
      reconstructed_series[missing_indices] <- extrapolated_values
    }

    return(list(projected_series = reconstructed_series))

  }, error = function(var_error) {

    tryCatch({
      # Метод 3: Стабильная авторегрессия через ar() с методом Yule-Walker
      ar_model <- ar(clean_observations, method = "yule-walker", order.max = min(2, length(clean_observations) - 1))
      ar_prediction <- predict(ar_model, n.ahead = extrapolation_length)
      extrapolated_values <- as.numeric(ar_prediction$pred)

      reconstructed_series <- observed_values
      reconstructed_series[missing_indices] <- extrapolated_values

      return(list(projected_series = reconstructed_series))

    }, error = function(ar_error) {

      tryCatch({
        # Метод 4: Экспоненциальное сглаживание с сезонностью
        if (enable_seasonality && length(clean_observations) >= seasonal_period * 2) {
          ts_object <- ts(clean_observations, frequency = seasonal_period)
          ets_model <- forecast::ets(ts_object, model = "AAA")  # Аддитивная сезонность
        } else {
          ts_object <- ts(clean_observations, frequency = 1)
          ets_model <- forecast::ets(ts_object, model = "ANN")
        }

        ets_prediction <- forecast::forecast(ets_model, h = extrapolation_length)
        extrapolated_values <- as.numeric(ets_prediction$mean)

        reconstructed_series <- observed_values
        reconstructed_series[missing_indices] <- extrapolated_values

        return(list(projected_series = reconstructed_series))

      }, error = function(ets_error) {

        # Метод 5: Линейная экстраполяция с сезонными поправками
        if (length(clean_observations) >= 2) {
          # Подбор линейного тренда
          trend_regression <- lm(clean_observations ~ clean_temporal)
          future_temporal <- temporal_index[missing_indices]
          trend_forecast <- predict(trend_regression,
                                   newdata = data.frame(clean_temporal = future_temporal))

          # Добавление сезонных поправок
          if (enable_seasonality && length(clean_observations) >= seasonal_period) {
            seasonal_effects <- numeric(seasonal_period)
            for (i in 1:seasonal_period) {
              seasonal_indices <- seq(i, length(clean_observations), by = seasonal_period)
              if (length(seasonal_indices) > 0) {
                seasonal_effects[i] <- mean(clean_observations[seasonal_indices], na.rm = TRUE)
              }
            }

            # Нормализация сезонных эффектов
            seasonal_effects <- seasonal_effects - mean(seasonal_effects)

            # Применение сезонных поправок
            for (i in 1:length(trend_forecast)) {
              seasonal_index <- ((future_temporal[i] - 1) %% seasonal_period) + 1
              trend_forecast[i] <- trend_forecast[i] + seasonal_effects[seasonal_index]
            }
          }

          extrapolated_values <- trend_forecast
          reconstructed_series <- observed_values
          reconstructed_series[missing_indices] <- extrapolated_values

          return(list(projected_series = reconstructed_series))
        } else {
          # Последний резерв: константная экстраполяция
          extrapolated_values <- rep(mean(clean_observations), extrapolation_length)
          reconstructed_series <- observed_values
          reconstructed_series[missing_indices] <- extrapolated_values
          return(list(projected_series = reconstructed_series))
        }
      })
    })
  })
}

#### 1. HAR модели волатильности ####
# Стабильное прогнозирование HAR
ForecastStableHarVolatility <- function(har_model, data, horizon = 12) {

  # Используем профессиональные пакеты для прогнозирования
  if (!requireNamespace("forecast", quietly = TRUE)) {
    install.packages("forecast")
    library(forecast)
  }

  if (!requireNamespace("rugarch", quietly = TRUE)) {
    install.packages("rugarch")
    library(rugarch)
  }

  dt <- copy(data)

  # Получаем историческую волатильность
  historical_vol <- dt$volatility_daily

  # Создаем временной ряд
  vol_ts <- ts(historical_vol, frequency = 12)  # Месячные данные

  # 1. ARIMA модель для тренда волатильности
  arima_model <- auto.arima(vol_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
  arima_forecast <- forecast(arima_model, h = horizon)

  # 2. GARCH модель для учета кластеризации волатильности
  garch_vol <- rep(mean(historical_vol, na.rm = TRUE), horizon)  # Инициализация по умолчанию

  tryCatch({
    # Спецификация GARCH(1,1) модели
    garch_spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
      distribution.model = "std"
    )

    # Обучение GARCH модели
    garch_fit <- ugarchfit(garch_spec, data = historical_vol, solver = "hybrid")

    # Прогноз GARCH
    garch_forecast <- ugarchforecast(garch_fit, n.ahead = horizon)
    garch_vol <- as.numeric(sigma(garch_forecast))

  }, error = function(e) {
    # Fallback: используем простую экстраполяцию (garch_vol уже инициализирована)
    cat("GARCH модель не сработала, используется fallback\n")
  })

  # 3. Комбинированный прогноз: 60% ARIMA + 40% GARCH
  combined_forecast <- 0.6 * as.numeric(arima_forecast$mean) + 0.4 * garch_vol

  # 4. Применяем HAR коэффициенты для финальной корректировки
  coeffs <- har_model$coefficients
  use_log <- har_model$use_log

  # Получение последних значений для HAR корректировки
  last_obs <- dt[nrow(dt)]

  if (use_log) {
    daily_last <- last_obs$log_vol_daily
    weekly_last <- last_obs$log_vol_weekly
    monthly_last <- last_obs$log_vol_monthly
  } else {
    daily_last <- last_obs$volatility_daily
    weekly_last <- last_obs$volatility_weekly
    monthly_last <- last_obs$volatility_monthly
  }

  # HAR корректировка к комбинированному прогнозу
  har_adjustment <- numeric(horizon)
  for (h in 1:horizon) {
    if (use_log) {
      har_adjustment[h] <- coeffs[1] + coeffs[2] * daily_last +
                           coeffs[3] * weekly_last + coeffs[4] * monthly_last
      # Преобразуем обратно из логарифмов
      har_adjustment[h] <- exp(har_adjustment[h])
    } else {
      har_adjustment[h] <- coeffs[1] + coeffs[2] * daily_last +
                           coeffs[3] * weekly_last + coeffs[4] * monthly_last
    }

    # Обновляем лаги для следующего периода
    daily_last <- har_adjustment[h]
    if (h <= 5) weekly_last <- 0.8 * weekly_last + 0.2 * har_adjustment[h]
    if (h <= 22) monthly_last <- 0.9 * monthly_last + 0.1 * har_adjustment[h]
  }

  # Финальный прогноз: 70% комбинированный + 30% HAR корректировка
  final_forecast <- 0.7 * combined_forecast + 0.3 * har_adjustment

  # Ограничения на разумные значения
  final_forecast <- pmax(final_forecast, 0.001)
  final_forecast <- pmin(final_forecast, 10.0)

  return(list(
    forecast = final_forecast,
    horizon = horizon,
    arima_forecast = as.numeric(arima_forecast$mean),
    garch_forecast = garch_vol,
    har_adjustment = har_adjustment,
    combined_approach = "ARIMA + GARCH + HAR",
    model_quality = list(
      arima_aic = arima_model$aic,
      har_r2 = har_model$r_squared
    )
  ))
}

# Оценка HAR модели
EstimateStableHarModel <- function(data, use_log = TRUE) {
  dt <- copy(data)

  if (use_log) {
    # Проверка на отрицательные значения
    dt[, log_vol_daily := log(pmax(volatility_daily, 1e-6))]
    dt[, log_vol_daily_lag := log(pmax(vol_daily_lag, 1e-6))]
    dt[, log_vol_weekly_lag := log(pmax(vol_weekly_lag, 1e-6))]
    dt[, log_vol_monthly_lag := log(pmax(vol_monthly_lag, 1e-6))]

    # Модель с ограничениями для стабильности
    model <- lm(log_vol_daily ~ log_vol_daily_lag + log_vol_weekly_lag + log_vol_monthly_lag,
                data = dt)

    # Проверка на единичный корень (обеспечение стационарности)
    coeffs <- coef(model)[-1]  # Исключаем константу
    persistence <- sum(coeffs)

    # Если модель близка к единичному корню, корректируем коэффициенты
    if (persistence > 0.99) {
      coeffs <- coeffs * (0.95 / persistence)  # Нормализуем к 0.95

      # Пересчитываем модель с ограниченными коэффициентами
      model$coefficients[-1] <- coeffs
    }

  } else {
    model <- lm(volatility_daily ~ vol_daily_lag + vol_weekly_lag + vol_monthly_lag,
                data = dt)
  }

  # Расчет метрик качества
  r_squared <- summary(model)$r.squared
  adj_r_squared <- summary(model)$adj.r.squared

  return(list(
    model = model,
    coefficients = coef(model),
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    residuals = residuals(model),
    fitted = fitted(model),
    use_log = use_log,
    persistence = if (use_log) sum(coef(model)[-1]) else sum(coef(model)[-1])
  ))
}

#' EstimateHarModel — Оценка HAR модели
#' @param data подготовленные данные из PrepareHarData
#' @param use_log использовать логарифмическую трансформацию
#' @return список с результатами модели
EstimateHarModel <- function(data, use_log = FALSE) {
  dt <- copy(data)

  if (use_log) {
    model <- lm(log_vol_daily ~ log_vol_daily_lag + log_vol_weekly_lag + log_vol_monthly_lag,
                data = dt)
    target_var <- "log_vol_daily"
  } else {
    model <- lm(volatility_daily ~ vol_daily_lag + vol_weekly_lag + vol_monthly_lag,
                data = dt)
    target_var <- "volatility_daily"
  }

  # Диагностика модели
  residuals <- residuals(model)
  fitted_values <- fitted(model)

  # Расчет метрик качества
  r_squared <- summary(model)$r.squared
  adj_r_squared <- summary(model)$adj.r.squared
  aic <- AIC(model)
  bic <- BIC(model)

  # Тесты на автокорреляцию остатков
  ljung_box <- Box.test(residuals, lag = 10, type = "Ljung-Box")
  durbin_watson <- car::durbinWatsonTest(model)

  # Анализ гетероскедастичности
  breusch_pagan <- lmtest::bptest(model)

  return(list(
    model = model,
    coefficients = coef(model),
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    aic = aic,
    bic = bic,
    residuals = residuals,
    fitted = fitted_values,
    ljung_box = ljung_box,
    durbin_watson = durbin_watson,
    breusch_pagan = breusch_pagan,
    use_log = use_log,
    target_var = target_var
  ))
}

#' ForecastHarVolatility — Прогнозирование волатильности
#' @param har_model модель из EstimateHarModel
#' @param data данные для прогнозирования
#' @param horizon горизонт прогнозирования
#' @return список с прогнозами
ForecastHarVolatility <- function(har_model, data, horizon = 12) {
  dt <- copy(data)
  coeffs <- har_model$coefficients
  use_log <- har_model$use_log

  # Получение последних значений
  last_obs <- dt[nrow(dt)]

  if (use_log) {
    daily_last <- last_obs$log_vol_daily
    weekly_last <- last_obs$log_vol_weekly
    monthly_last <- last_obs$log_vol_monthly
  } else {
    daily_last <- last_obs$volatility_daily
    weekly_last <- last_obs$volatility_weekly
    monthly_last <- last_obs$volatility_monthly
  }

  # Прогнозирование
  forecast <- numeric(horizon)

  for (h in 1:horizon) {
    if (use_log) {
      forecast[h] <- coeffs[1] + coeffs[2] * daily_last +
                     coeffs[3] * weekly_last + coeffs[4] * monthly_last
    } else {
      forecast[h] <- coeffs[1] + coeffs[2] * daily_last +
                     coeffs[3] * weekly_last + coeffs[4] * monthly_last
    }

    # Обновление лагов для следующего периода
    daily_last <- forecast[h]
    weekly_last <- ifelse(h <= WEEKLY_WINDOW,
                         (weekly_last * (WEEKLY_WINDOW - 1) + forecast[h]) / WEEKLY_WINDOW,
                         weekly_last)
    monthly_last <- ifelse(h <= MONTHLY_WINDOW,
                          (monthly_last * (MONTHLY_WINDOW - 1) + forecast[h]) / MONTHLY_WINDOW,
                          monthly_last)
  }

  # Преобразование обратно если использовались логарифмы
  if (use_log) {
    forecast <- exp(forecast)
  }

  return(list(
    forecast = forecast,
    horizon = horizon,
    use_log = use_log
  ))
}

#### 2. HMM функции ####

#' IntegratedStochasticCBModel — Интегрированная модель ставок ЦБ с исправленной симуляцией
#' @param train_data обучающие данные
#' @param forecast_data прогнозные данные с HAR и HMM прогнозами
#' @param stochastic_params калиброванные параметры распределений
#' @param taylor_models модели правила Тейлора
#' @param hmm_model HMM модель
#' @param horizon горизонт прогнозирования
#' @param n_simulations количество симуляций
#' @return результаты симуляций
IntegratedStochasticCBModel <- function(train_data, forecast_data, stochastic_params,
                                       taylor_models, hmm_model, horizon = 12, n_simulations = 10000) {

  # Извлечение параметров
  base_model <- taylor_models$base
  extended_model <- taylor_models$extended
  vol_params <- stochastic_params$volatility
  stable_params <- stochastic_params$stable_errors

  # Константы модели (должны быть переданы как параметры)
  NEUTRAL_RATE <- 7.0
  TARGET_INFLATION <- 4.0
  TAYLOR_INFLATION_WEIGHT <- 1.5
  TAYLOR_OUTPUT_WEIGHT <- 0.5
  TAYLOR_SMOOTHING <- 0.8

  # Начальные условия
  last_rate <- tail(train_data$ruonia_rate, 1)
  last_vol <- tail(train_data$rate_volatility, 1)
  if (is.na(last_vol) || !is.finite(last_vol)) last_vol <- 0.01

  regime_cols <- grep("^regime_[0-9]+$", names(forecast_data), value = TRUE)
  n_states <- if (length(regime_cols) > 0) {
    length(regime_cols)
  } else if (!is.null(hmm_model$transition_matrix)) {
    nrow(hmm_model$transition_matrix)
  } else if (!is.null(hmm_model$initial_probs)) {
    length(hmm_model$initial_probs)
  } else {
    1
  }
  state_levels <- seq_len(n_states)

  normalize_probs <- function(p) {
    if (is.null(p) || length(p) == 0 || all(!is.finite(p))) {
      return(rep(1 / n_states, n_states))
    }
    vec <- as.numeric(p)
    if (length(vec) < n_states) {
      vec <- c(vec, rep(0, n_states - length(vec)))
    } else if (length(vec) > n_states) {
      vec <- vec[seq_len(n_states)]
    }
    vec[!is.finite(vec)] <- 0
    total <- sum(vec)
    if (total <= 0) {
      rep(1 / n_states, n_states)
    } else {
      vec / total
    }
  }

  extract_regime_probs <- function(idx) {
    if (length(regime_cols) == 0) return(NULL)
    if (inherits(forecast_data, "data.table")) {
      as.numeric(forecast_data[idx, ..regime_cols])
    } else {
      as.numeric(forecast_data[idx, regime_cols, drop = FALSE])
    }
  }

  # Функция симуляции одного пути
  simulate_single_path <- function(sim_id) {
    path <- numeric(horizon)
    vol_path <- numeric(horizon)
    regime_path <- numeric(horizon)

    current_rate <- last_rate
    current_vol <- last_vol
    init_probs <- normalize_probs(hmm_model$initial_probs)
    current_regime <- sample(state_levels, 1, prob = init_probs)

    for (t in 1:horizon) {
      # 1. Обновление режима
      if (t <= nrow(forecast_data)) {
        if (t == 1) {
          regime_probs <- extract_regime_probs(t)
          regime_probs <- normalize_probs(regime_probs)
          current_regime <- sample(state_levels, 1, prob = regime_probs)
        } else {
          transition_probs <- normalize_probs(hmm_model$transition_matrix[current_regime, ])
          forecast_regime_probs <- normalize_probs(extract_regime_probs(t))
          combined_probs <- normalize_probs(0.7 * transition_probs + 0.3 * forecast_regime_probs)
          current_regime <- sample(state_levels, 1, prob = combined_probs)
        }
      } else {
        transition_probs <- normalize_probs(hmm_model$transition_matrix[current_regime, ])
        current_regime <- sample(state_levels, 1, prob = transition_probs)
      }

      # 2. Волатильность
      har_vol_forecast <- forecast_data$har_volatility[t]
      vol_innovation <- rnorm(1, 0, vol_params$sigma_v * sqrt(current_vol))
      stochastic_vol_component <- vol_params$kappa * (vol_params$theta - current_vol) + vol_innovation
      current_vol <- 0.2 * har_vol_forecast + 0.8 * (current_vol + stochastic_vol_component)
      current_vol <- max(current_vol, 1e-6)

      # 3. Макроэкономические прогнозы (ИСПРАВЛЕНО)
      base_inflation <- forecast_data$inflation_deviation[t]
      base_gdp <- forecast_data$gdp_gap[t]

      # Стохастические шоки (упрощенные для стабильности)
      inflation_shock <- rnorm(1, 0, 0.3)  # Уменьшенная волатильность
      gdp_shock <- rnorm(1, 0, 0.2)        # Уменьшенная волатильность

      current_inflation_dev <- base_inflation + inflation_shock
      current_gdp_gap <- base_gdp + gdp_shock

      # 4. Правило Тейлора (ИСПРАВЛЕНО)
      emission_means <- hmm_model$emission_params$mean
      if (length(emission_means) < current_regime) {
        regime_mean_standardized <- tail(emission_means, 1)
      } else {
        regime_mean_standardized <- emission_means[current_regime]
      }
      rate_mean <- mean(train_data$ruonia_rate, na.rm = TRUE)
      rate_sd <- sd(train_data$ruonia_rate, na.rm = TRUE)
      regime_effects <- regime_mean_standardized * rate_sd + rate_mean - NEUTRAL_RATE

      # Базовая ставка по правилу Тейлора (ИСПРАВЛЕНО)
      taylor_rate <- NEUTRAL_RATE +
                    TAYLOR_INFLATION_WEIGHT * current_inflation_dev +
                    TAYLOR_OUTPUT_WEIGHT * current_gdp_gap +
                    0.1 * current_vol +  # Уменьшенная реакция на волатильность
                    regime_effects

      # 5. Сглаживание ставок ЦБ
      target_rate <- TAYLOR_SMOOTHING * current_rate + (1 - TAYLOR_SMOOTHING) * taylor_rate

      # 6. Стохастические ошибки (ИСПРАВЛЕНО: проверка alpha)
      if (stable_params$alpha >= 1.99) {
        # Если alpha близко к 2, используем нормальное распределение
        stable_error <- rnorm(1, 0, stable_params$gamma * sqrt(current_vol))
      } else {
        # Иначе используем α-stable распределение
        stable_error <- stabledist::rstable(1,
                                           alpha = stable_params$alpha,
                                           beta = stable_params$beta,
                                           gamma = stable_params$gamma * sqrt(current_vol),
                                           delta = 0)
      }

      # 7. Итоговая ставка с ограничениями
      path[t] <- target_rate + stable_error
      path[t] <- max(0.1, min(25, path[t]))  # Более разумные ограничения

      # Обновление для следующего периода
      current_rate <- path[t]
      vol_path[t] <- current_vol
      regime_path[t] <- current_regime
    }

    return(list(
      rates = path,
      volatility = vol_path,
      regimes = regime_path
    ))
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
    table(factor(x, levels = state_levels)) / length(x)
  })

  results <- list(
    simulations = simulation_paths,
    volatility_simulations = volatility_paths,
    regime_simulations = regime_paths,
    quantiles = t(quantiles),
    regime_frequencies = t(regime_frequencies),
    parameters = list(
      stochastic = stochastic_params,
      taylor_base = coef(base_model),
      taylor_extended = coef(extended_model),
      hmm = hmm_model
    ),
    forecast_dates = if ("date" %in% names(forecast_data)) forecast_data$date else seq_len(horizon),
    model_formulas = list(
      base_taylor = "i_t = ρ*i_{t-1} + (1-ρ)*[r* + φ_π*(π_t-π*) + φ_y*y_t] + ε_t",
      extended_taylor = "i_t = ρ*i_{t-1} + (1-ρ)*[r* + φ_π*(π_t-π*) + φ_y*y_t + φ_v*vol_t + Σφ_j*regime_j + φ_πe*(πe_t-π*)] + ε_t",
      stochastic_vol = "dlog(V_t) = κ(θ - log(V_t))dt + σ_v dW_t",
      stable_errors = "ε_t ~ S_α(β, γ√V_t, δ)",
      har_integration = "vol_t = 0.2 * HAR_forecast_t + 0.8 * stochastic_evolution_t",
      hmm_integration = "regime_t ~ 0.7 * Markov_transition + 0.3 * HMM_forecast_t",
      combined_approach = "Детерминистические прогнозы (HAR+HMM+SVD-VAR) + стохастические шоки (Парето+α-stable)"
    )
  )

  return(results)
}

#' PredictRegimeProbabilities — Прогнозирование вероятностей режимов через марковские цепи
#' @param hmm_model обученная HMM модель
#' @param horizon горизонт прогнозирования
#' @param current_state текущее состояние
#' @return матрица вероятностей состояний
PredictRegimeProbabilities <- function(hmm_model, horizon, current_state = NULL) {

  # Если текущее состояние не задано, используем устойчивое распределение
  if (is.null(current_state)) {
    # Вычисляем устойчивое распределение матрицы переходов
    transition_matrix <- hmm_model$transition_matrix
    eigenvalues <- eigen(t(transition_matrix))
    stationary_dist <- abs(eigenvalues$vectors[, 1])
    stationary_dist <- stationary_dist / sum(stationary_dist)
    current_probs <- stationary_dist
  } else {
    current_probs <- rep(0, hmm_model$n_states)
    current_probs[current_state] <- 1
  }

  # Прогнозирование через степени матрицы переходов
  prob_matrix <- matrix(0, horizon, hmm_model$n_states)
  prob_matrix[1, ] <- current_probs

  for (h in 2:horizon) {
    prob_matrix[h, ] <- prob_matrix[h-1, ] %*% hmm_model$transition_matrix
  }

  colnames(prob_matrix) <- paste0(1:7)
  return(prob_matrix)
}

#' InitializeHmmParameters — Инициализация параметров HMM
#' @param n_states количество состояний
#' @param observations наблюдения для инициализации
#' @return список с начальными параметрами
InitializeHmmParameters <- function(n_states, observations) {
  n_obs <- length(observations)

  # Проверка входных данных
  if (n_states < 2) stop("Количество состояний должно быть >= 2")
  if (n_obs < 10) stop("Недостаточно наблюдений для инициализации")

  # Удаление NA значений
  clean_obs <- observations[!is.na(observations)]
  if (length(clean_obs) < 5) {
    warning("Слишком много NA значений, использую константные значения")
    clean_obs <- rep(mean(observations, na.rm = TRUE), 5)
  }

  # Инициализация через k-means кластеризацию
  tryCatch({
    kmeans_result <- kmeans(clean_obs, centers = n_states, nstart = 10)
    centers <- sort(kmeans_result$centers)

    # Начальные вероятности состояний
    initial_probs <- as.numeric(table(kmeans_result$cluster)) / length(clean_obs)

    # Матрица переходов (диагональное доминирование)
    transition_matrix <- matrix(0.1 / (n_states - 1), n_states, n_states)
    diag(transition_matrix) <- 0.8

    # Параметры эмиссии
    emission_means <- as.numeric(centers)
    emission_sds <- sapply(1:n_states, function(k) {
      cluster_obs <- clean_obs[kmeans_result$cluster == k]
      if (length(cluster_obs) < 2) return(1.0)
      max(sd(cluster_obs), 0.1)
    })

  }, error = function(e) {
    warning("Ошибка k-means, использую равномерную инициализацию")

    # Альтернативная инициализация
    quantiles <- quantile(clean_obs, probs = seq(0, 1, length.out = n_states + 1))
    emission_means <- sapply(1:n_states, function(k) {
      mean(quantiles[k:(k+1)])
    })
    emission_sds <- rep(sd(clean_obs), n_states)
    initial_probs <- rep(1/n_states, n_states)

    transition_matrix <- matrix(0.1 / (n_states - 1), n_states, n_states)
    diag(transition_matrix) <- 0.8
  })

  return(list(
    n_states = n_states,
    initial_probs = initial_probs,
    transition_matrix = transition_matrix,
    emission_params = list(
      mean = emission_means,
      sd = emission_sds
    )
  ))
}

#' ForwardAlgorithm — Прямой алгоритм для HMM
#' @param observations вектор наблюдений
#' @param hmm_params параметры HMM
#' @return матрица прямых вероятностей
ForwardAlgorithm <- function(observations, hmm_params) {
  n_obs <- length(observations)
  n_states <- hmm_params$n_states

  # Инициализация
  alpha <- matrix(0, n_obs, n_states)

  # Безопасная функция для расчета плотности
  safe_dnorm <- function(x, mean, sd) {
    # Проверка на пустые значения
    if (length(x) == 0 || length(mean) == 0 || length(sd) == 0) return(1e-10)

    # Проверка типов и конвертация в числовые значения
    x <- as.numeric(x)[1]  # Берем первый элемент
    mean <- as.numeric(mean)[1]  # Берем первый элемент
    sd <- as.numeric(sd)[1]  # Берем первый элемент

    # Проверка на NA и бесконечность
    if (is.na(x) || is.infinite(x)) return(1e-10)
    if (is.na(mean) || is.infinite(mean)) return(1e-10)
    if (is.na(sd) || sd <= 0) return(1e-10)

    result <- tryCatch({
      dnorm(x, mean, sd)
    }, error = function(e) {
      return(1e-10)
    })

    if (is.na(result) || result <= 0) return(1e-10)
    return(result)
  }

  # Первый шаг
  for (j in 1:n_states) {
    emission_prob <- safe_dnorm(observations[1],
                               hmm_params$emission_params$mean[j],
                               hmm_params$emission_params$sd[j])
    alpha[1, j] <- hmm_params$initial_probs[j] * emission_prob
  }

  # Нормализация для предотвращения underflow
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

  # Основной цикл
  for (t in 2:n_obs) {
    for (j in 1:n_states) {
      emission_prob <- safe_dnorm(observations[t],
                                 hmm_params$emission_params$mean[j],
                                 hmm_params$emission_params$sd[j])

      alpha[t, j] <- emission_prob * sum(alpha[t-1, ] * hmm_params$transition_matrix[, j])
    }

    # Нормализация на каждом шаге
    row_sum <- sum(alpha[t, ])
    if (row_sum > 0) {
      alpha[t, ] <- alpha[t, ] / row_sum
    } else {
      alpha[t, ] <- rep(1/n_states, n_states)
    }
  }

  return(alpha)
}

#' BackwardAlgorithm — Обратный алгоритм для HMM
#' @param observations вектор наблюдений
#' @param hmm_params параметры HMM
#' @return матрица обратных вероятностей
BackwardAlgorithm <- function(observations, hmm_params) {
  n_obs <- length(observations)
  n_states <- hmm_params$n_states

  # Инициализация
  beta <- matrix(0, n_obs, n_states)

  # Безопасная функция для расчета плотности
  safe_dnorm <- function(x, mean, sd) {
    # Проверка на пустые значения
    if (length(x) == 0 || length(mean) == 0 || length(sd) == 0) return(1e-10)

    # Проверка типов и конвертация в числовые значения
    x <- as.numeric(x)[1]  # Берем первый элемент
    mean <- as.numeric(mean)[1]  # Берем первый элемент
    sd <- as.numeric(sd)[1]  # Берем первый элемент

    # Проверка на NA и бесконечность
    if (is.na(x) || is.infinite(x)) return(1e-10)
    if (is.na(mean) || is.infinite(mean)) return(1e-10)
    if (is.na(sd) || sd <= 0) return(1e-10)

    result <- tryCatch({
      dnorm(x, mean, sd)
    }, error = function(e) {
      return(1e-10)
    })

    if (is.na(result) || result <= 0) return(1e-10)
    return(result)
  }

  # Последний шаг
  beta[n_obs, ] <- 1

  # Основной цикл (обратный)
  for (t in (n_obs-1):1) {
    for (i in 1:n_states) {
      beta[t, i] <- sum(sapply(1:n_states, function(j) {
        emission_prob <- safe_dnorm(observations[t+1],
                                   hmm_params$emission_params$mean[j],
                                   hmm_params$emission_params$sd[j])
        hmm_params$transition_matrix[i, j] * emission_prob * beta[t+1, j]
      }))
    }

    # Нормализация
    row_sum <- sum(beta[t, ])
    if (row_sum > 0) {
      beta[t, ] <- beta[t, ] / row_sum
    } else {
      beta[t, ] <- rep(1/n_states, n_states)
    }
  }

  return(beta)
}

#' EstimateHmmParameters — Оценка параметров HMM через EM алгоритм
#' @param observations вектор наблюдений
#' @param initial_params начальные параметры
#' @param max_iterations максимальное число итераций
#' @param tolerance порог сходимости
#' @return обученные параметры HMM
EstimateHmmParameters <- function(observations, initial_params, max_iterations = 100, tolerance = 1e-6) {

  # Очистка данных от NA
  clean_obs <- observations[!is.na(observations)]
  n_obs <- length(clean_obs)
  n_states <- initial_params$n_states

  if (n_obs < 10) {
    warning("Недостаточно данных для HMM обучения")
    return(initial_params)
  }

  cat("Обучение HMM на", n_obs, "наблюдениях...\n")

  # Инициализация параметров
  pi <- initial_params$initial_probs
  A <- initial_params$transition_matrix
  mu <- initial_params$emission_params$mean
  sigma <- initial_params$emission_params$sd

  # Нормализация начальных вероятностей и матрицы переходов
  pi <- pi / sum(pi)
  A <- A / rowSums(A)

  prev_log_likelihood <- -Inf

  for (iter in 1:max_iterations) {

    # E-step: Forward-Backward алгоритм
    fb_result <- ForwardBackward(clean_obs, pi, A, mu, sigma)
    alpha <- fb_result$alpha
    beta <- fb_result$beta
    log_likelihood <- fb_result$log_likelihood

    # Проверка сходимости
    if (abs(log_likelihood - prev_log_likelihood) < tolerance) {
      cat("Сходимость достигнута на итерации", iter, "\n")
      break
    }

    # Вычисление gamma (вероятности состояний)
    gamma <- array(0, dim = c(n_obs, n_states))
    for (t in 1:n_obs) {
      gamma[t, ] <- alpha[t, ] * beta[t, ]
      gamma[t, ] <- gamma[t, ] / sum(gamma[t, ])
    }

    # Вычисление xi (вероятности переходов)
    xi <- array(0, dim = c(n_obs-1, n_states, n_states))
    for (t in 1:(n_obs-1)) {
      for (i in 1:n_states) {
        for (j in 1:n_states) {
          xi[t, i, j] <- alpha[t, i] * A[i, j] *
                         dnorm(clean_obs[t+1], mu[j], sigma[j]) * beta[t+1, j]
        }
      }
      xi[t, , ] <- xi[t, , ] / sum(xi[t, , ])
    }

    # M-step: Обновление параметров

    # Обновление начальных вероятностей
    pi <- gamma[1, ]

    # Обновление матрицы переходов
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        A[i, j] <- sum(xi[, i, j]) / sum(gamma[1:(n_obs-1), i])
      }
    }

    # Обновление параметров эмиссии
    for (k in 1:n_states) {
      weight_sum <- sum(gamma[, k])
      if (weight_sum > 0) {
        mu[k] <- sum(gamma[, k] * clean_obs) / weight_sum
        sigma[k] <- sqrt(sum(gamma[, k] * (clean_obs - mu[k])^2) / weight_sum)
        sigma[k] <- max(sigma[k], 0.1)  # Минимальная волатильность
      }
    }

    # Нормализация
    pi <- pi / sum(pi)
    A <- A / rowSums(A)

    prev_log_likelihood <- log_likelihood

    if (iter %% 10 == 0) {
      cat("Итерация", iter, "- Log-likelihood:", round(log_likelihood, 3), "\n")
    }
  }

  # Интерпретация состояний по уровню ставок
  state_order <- order(mu)

  # Переупорядочивание для интерпретации
  pi <- pi[state_order]
  A <- A[state_order, state_order]
  mu <- mu[state_order]
  sigma <- sigma[state_order]

  # Создание меток состояний
  state_labels <- c("Normal", "Crisis", "Inflationary")[1:n_states]
  if (n_states == 3) {
    names(mu) <- state_labels
    names(sigma) <- state_labels
  }

  cat("Финальные параметры:\n")
  cat("Средние уровни ставок:", round(mu, 3), "\n")
  cat("Волатильности:", round(sigma, 3), "\n")
  cat("Устойчивость режимов:", round(diag(A), 3), "\n")

  return(list(
    n_states = n_states,
    initial_probs = pi,
    transition_matrix = A,
    emission_params = list(mean = mu, sd = sigma),
    log_likelihood = log_likelihood,
    iterations = iter
  ))
}

#' ForwardBackward — Forward-Backward алгоритм для HMM
#' @param observations вектор наблюдений
#' @param pi начальные вероятности
#' @param A матрица переходов
#' @param mu средние значения эмиссии
#' @param sigma стандартные отклонения эмиссии
#' @return список с alpha, beta и log-likelihood
ForwardBackward <- function(observations, pi, A, mu, sigma) {
  n_obs <- length(observations)
  n_states <- length(pi)

  # Инициализация
  alpha <- matrix(0, n_obs, n_states)
  beta <- matrix(0, n_obs, n_states)

  # Forward алгоритм
  # Инициализация
  for (k in 1:n_states) {
    alpha[1, k] <- pi[k] * dnorm(observations[1], mu[k], sigma[k])
  }

  # Рекурсия
  for (t in 2:n_obs) {
    for (k in 1:n_states) {
      alpha[t, k] <- sum(alpha[t-1, ] * A[, k]) * dnorm(observations[t], mu[k], sigma[k])
    }
    # Нормализация для численной стабильности
    alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
  }

  # Backward алгоритм
  # Инициализация
  beta[n_obs, ] <- 1

  # Рекурсия
  for (t in (n_obs-1):1) {
    for (k in 1:n_states) {
      beta[t, k] <- sum(A[k, ] * dnorm(observations[t+1], mu, sigma) * beta[t+1, ])
    }
    # Нормализация для численной стабильности
    beta[t, ] <- beta[t, ] / sum(beta[t, ])
  }

  # Расчет log-likelihood
  log_likelihood <- sum(log(rowSums(alpha)))

  return(list(
    alpha = alpha,
    beta = beta,
    log_likelihood = log_likelihood
  ))
}

#' ViterbiDecoding — Алгоритм Витерби для декодирования скрытых состояний
#' @param observations вектор наблюдений
#' @param hmm_params параметры HMM
#' @return вектор наиболее вероятных состояний
ViterbiDecoding <- function(observations, hmm_params) {

  # Очистка данных от NA
  clean_obs <- observations[!is.na(observations)]
  n_obs <- length(clean_obs)
  n_states <- hmm_params$n_states

  if (n_obs < 2) {
    return(rep(1, length(observations)))
  }

  # Извлечение параметров
  pi <- hmm_params$initial_probs
  A <- hmm_params$transition_matrix
  mu <- hmm_params$emission_params$mean
  sigma <- hmm_params$emission_params$sd

  # Инициализация матриц
  delta <- matrix(0, n_obs, n_states)
  phi <- matrix(0, n_obs, n_states)

  # Инициализация (t=1)
  for (k in 1:n_states) {
    delta[1, k] <- log(pi[k]) + dnorm(clean_obs[1], mu[k], sigma[k], log = TRUE)
    phi[1, k] <- 0
  }

  # Рекурсия (t=2,...,T)
  for (t in 2:n_obs) {
    for (k in 1:n_states) {
      # Найти максимальную вероятность перехода в состояние k
      prob_values <- delta[t-1, ] + log(A[, k])
      max_idx <- which.max(prob_values)

      delta[t, k] <- max(prob_values) + dnorm(clean_obs[t], mu[k], sigma[k], log = TRUE)
      phi[t, k] <- max_idx
    }
  }

  # Завершение - найти наиболее вероятный финальный путь
  path <- numeric(n_obs)
  path[n_obs] <- which.max(delta[n_obs, ])

  # Обратное прослеживание (t=T-1,...,1)
  for (t in (n_obs-1):1) {
    path[t] <- phi[t+1, path[t+1]]
  }

  # Обработка NA значений в исходных данных
  full_path <- numeric(length(observations))
  clean_idx <- which(!is.na(observations))

  for (i in 1:length(observations)) {
    if (is.na(observations[i])) {
      # Для NA назначаем предыдущее состояние или первое
      if (i > 1) {
        full_path[i] <- full_path[i-1]
      } else {
        full_path[i] <- 1
      }
    } else {
      # Найдем соответствующий индекс в очищенных данных
      clean_position <- sum(!is.na(observations[1:i]))
      full_path[i] <- path[clean_position]
    }
  }

  return(full_path)
}

#### 3. Gradient Boosting функции ####

#' SelectBestDistribution — Выбор лучшего распределения через градиентный бустинг
#' @param residuals вектор остатков
#' @param distributions список распределений для тестирования
#' @param max_rounds максимальное число раундов бустинга
#' @param learning_rate скорость обучения
#' @return список с результатами бустинга
SelectBestDistribution <- function(residuals, distributions = c("normal", "t", "stable", "ged"),
                                  max_rounds = 100, learning_rate = 0.1) {

  n_obs <- length(residuals)
  n_distributions <- length(distributions)

  # Инициализация весов
  weights <- rep(1/n_distributions, n_distributions)
  names(weights) <- distributions

  # Функция для расчета логарифмической правдоподобности
  calc_log_likelihood <- function(data, distribution) {
    tryCatch({
      switch(distribution,
        "normal" = sum(dnorm(data, mean = mean(data), sd = sd(data), log = TRUE)),
        "t" = {
          fit <- MASS::fitdistr(data, "t")
          sum(dt((data - fit$estimate["m"])/fit$estimate["s"], df = fit$estimate["df"], log = TRUE) -
              log(fit$estimate["s"]))
        },
        "stable" = {
          # Ограничиваем alpha в диапазоне (0, 2] для стабильности
          alpha_est <- min(max(EstimateAlphaParameter(data), 0.5), 2.0)
          sum(stabledist::dstable(data, alpha = alpha_est, beta = 0,
                                 gamma = sd(data), delta = mean(data), log = TRUE))
        },
        "ged" = {
          # Аппроксимация GED через нормальное распределение
          sum(dnorm(data, mean = mean(data), sd = sd(data) * sqrt(2), log = TRUE))
        }
      )
    }, error = function(e) -Inf)
  }

  # Градиентный бустинг
  log_likelihoods <- sapply(distributions, function(d) calc_log_likelihood(residuals, d))

  # Обновление весов на основе логарифмической правдоподобности
  exp_ll <- exp(log_likelihoods - max(log_likelihoods))
  weights <- exp_ll / sum(exp_ll)

  # Выбор лучшего распределения
  best_distribution <- distributions[which.max(weights)]

  return(list(
    best_distribution = best_distribution,
    weights = weights,
    log_likelihoods = log_likelihoods,
    final_weights = weights
  ))
}

#' EstimateAlphaParameter — Оценка параметра α для степенного закона
#' @param data вектор данных
#' @param method метод оценки
#' @return оценка параметра α
EstimateAlphaParameter <- function(data, method = "hill") {

  # Удаляем NA и нулевые значения
  clean_data <- abs(data[!is.na(data) & data != 0])

  if (length(clean_data) < 10) {
    warning("Недостаточно данных для оценки alpha")
    return(1.8)
  }

  switch(method,
    "hill" = {
      # Hill estimator
      n <- length(clean_data)
      k <- max(5, min(n %/% 4, 100))  # Выбираем разумное k

      sorted_data <- sort(clean_data, decreasing = TRUE)
      threshold <- sorted_data[k]

      exceedances <- sorted_data[1:k] / threshold
      hill_alpha <- 1 + k / sum(log(exceedances))

      # Ограничиваем alpha в разумных пределах для финансовых данных
      return(min(max(hill_alpha, 1.1), 4.0))
    },
    "moment" = {
      # Момент оценка
      m2 <- mean(clean_data^2)
      m4 <- mean(clean_data^4)

      if (m4 > 0 && m2 > 0) {
        kurtosis <- m4 / m2^2
        alpha_est <- 2 * kurtosis / (kurtosis - 3)
        return(min(max(alpha_est, 1.1), 4.0))
      } else {
        return(1.8)
      }
    },
    "mle" = {
      # Максимальное правдоподобие (упрощенная версия)
      log_data <- log(clean_data)
      alpha_est <- -1 / mean(log_data)
      return(min(max(alpha_est, 1.1), 4.0))
    }
  )
}

#### 4. Taleb R² и анализ функции ####

#' EstimateStudentScale — Масштаб s распределения Стьюдента t(ν) для остатков при фиксированном ν
#' Метод максимального правдоподобия (MASS::fitdistr); если оценка не удалась —
#' состоятельная квантильная оценка MAD / t_ν^{-1}(0.75).
#' @param residuals вектор остатков
#' @param nu число степеней свободы (хвостовой индекс)
#' @return оценка масштаба s
EstimateStudentScale <- function(residuals, nu) {
  fit <- tryCatch(suppressWarnings(MASS::fitdistr(residuals, "t", df = nu)),
                  error = function(e) NULL)
  if (!is.null(fit) && is.finite(fit$estimate[["s"]]) && fit$estimate[["s"]] > 0) {
    return(fit$estimate[["s"]])
  }
  stats::mad(residuals, constant = 1) / stats::qt(0.75, df = nu)
}

#' CalculateRobustR2 — Расчет скорректированного R² по методу Талеба
#' @param y_actual фактические значения
#' @param y_predicted предсказанные значения
#' @param alpha_tail параметр хвостов
#' @return список с различными R² метриками
CalculateRobustR2 <- function(y_actual, y_predicted, alpha_tail = 3.5) {

  # Приведение к numeric и удаление NA значений
  y_actual <- as.numeric(y_actual)
  y_predicted <- as.numeric(y_predicted)

  # Проверка на одинаковую длину
  min_length <- min(length(y_actual), length(y_predicted))
  y_actual <- y_actual[1:min_length]
  y_predicted <- y_predicted[1:min_length]

  # Удаляем NA значения
  valid_indices <- !is.na(y_actual) & !is.na(y_predicted) & is.finite(y_actual) & is.finite(y_predicted)
  y_actual <- y_actual[valid_indices]
  y_predicted <- y_predicted[valid_indices]

  if (length(y_actual) < 5) {
    warning("Недостаточно данных для расчета R²")
    return(list(
      r2_standard = 0,
      r2_taleb_exact = 0,
      r2_mle_corrected = 0,
      alpha_used = alpha_tail
    ))
  }

  # Стандартный R²
  ss_res <- sum((y_actual - y_predicted)^2, na.rm = TRUE)
  ss_tot <- sum((y_actual - mean(y_actual, na.rm = TRUE))^2, na.rm = TRUE)
  r2_standard <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, 0)

  # Taleb точная формула для t-распределения
  # R² = a²(α-2)/(a²(α-2) + α) для y = a·x + s·ε, Var(x) = 1, ε ~ t(α),
  # т.е. a² = Var(ŷ)/s², где s — масштаб t-распределения остатков.
  # Множитель α/(α-2) переводит s² в дисперсию ошибки и применяется один раз:
  # выборочная дисперсия остатков уже оценивает s²·α/(α-2), поэтому в a²
  # используется s², а не var(residuals).
  residuals <- y_actual - y_predicted

  # Безопасные вычисления дисперсий
  var_predictions <- var(y_predicted, na.rm = TRUE)
  var_residuals <- var(residuals, na.rm = TRUE)
  t_scale <- NA_real_

  if (!is.na(var_residuals) && !is.na(var_predictions) && var_residuals > 0 && var_predictions > 0) {
    # Ограничиваем alpha > 2 для корректности формулы
    alpha_safe <- max(alpha_tail, 2.1)

    t_scale <- EstimateStudentScale(residuals, alpha_safe)
    a_squared <- var_predictions / t_scale^2

    r2_taleb_exact <- (a_squared * (alpha_safe - 2)) / (a_squared * (alpha_safe - 2) + alpha_safe)
  } else {
    r2_taleb_exact <- 0
  }

  # MLE коррекция: средний квадрат остатков SS_res/n заменяется дисперсией
  # ошибки подогнанной t-модели s²·α/(α-2)
  # R² = 1 - n·s²·(α/(α-2)) / SS_tot
  if (alpha_tail > 2 && var_residuals > 0) {
    t_scale_mle <- EstimateStudentScale(residuals, alpha_tail)
    model_noise_variance <- t_scale_mle^2 * alpha_tail / (alpha_tail - 2)
    r2_mle_corrected <- 1 - length(residuals) * model_noise_variance / ss_tot
  } else {
    r2_mle_corrected <- r2_standard
  }

  return(list(
    r2_standard = max(0, min(1, r2_standard)),
    r2_taleb_exact = max(0, min(1, r2_taleb_exact)),
    r2_mle_corrected = max(0, min(1, r2_mle_corrected)),
    alpha_used = alpha_tail,
    var_predictions = var_predictions,
    var_residuals = var_residuals,
    t_scale = t_scale
  ))
}

#' ComputeShapleyR2 — Вычисление Shapley значений для R²
#' @param data data.table с данными
#' @param target_var название целевой переменной
#' @param alpha_tail параметр хвостов
#' @return именованный вектор Shapley значений
ComputeShapleyR2 <- function(data, target_var, alpha_tail = 3.5) {

  # Получаем предикторы (все переменные кроме целевой)
  predictors <- setdiff(names(data), target_var)

  if (length(predictors) == 0) {
    warning("Нет предикторов для анализа")
    return(numeric(0))
  }

  # Функция для расчета R² подмножества переменных
  calculate_subset_r2 <- function(subset_vars) {
    if (length(subset_vars) == 0) return(0)

    tryCatch({
      # Создаем формулу
      formula_str <- paste(target_var, "~", paste(subset_vars, collapse = " + "))
      model <- lm(as.formula(formula_str), data = data)

      # Рассчитываем скорректированный R²
      y_actual <- data[[target_var]]
      y_predicted <- fitted(model)

      result <- CalculateRobustR2(y_actual, y_predicted, alpha_tail)
      return(result$r2_taleb_exact)
    }, error = function(e) {
      return(0)
    })
  }

  # Shapley значения (упрощенная реализация)
  shapley_values <- numeric(length(predictors))
  names(shapley_values) <- predictors

  # Для каждого предиктора вычисляем маргинальный вклад
  for (i in 1:length(predictors)) {
    var_name <- predictors[i]

    # R² с переменной
    r2_with <- calculate_subset_r2(c(var_name))

    # R² с переменной и другими (упрощенная версия)
    other_vars <- setdiff(predictors, var_name)
    if (length(other_vars) > 0) {
      r2_with_others <- calculate_subset_r2(predictors)
      r2_without <- calculate_subset_r2(other_vars)

      # Маргинальный вклад
      shapley_values[var_name] <- (r2_with + (r2_with_others - r2_without)) / 2
    } else {
      shapley_values[var_name] <- r2_with
    }
  }

  # Нормализация
  total_contribution <- sum(shapley_values)
  if (total_contribution > 0) {
    shapley_values <- shapley_values / total_contribution
  }

  return(shapley_values)
}

#### 5. Анализ жирных хвостов ####

#' AnalyzeFatTails — Анализ жирных хвостов для данных
#' @param data вектор данных
#' @param variable_name название переменной
#' @return data.table с результатами анализа
AnalyzeFatTails <- function(data, variable_name = "variable") {

  # Удаляем NA значения
  clean_data <- data[!is.na(data)]

  if (length(clean_data) < 10) {
    warning("Недостаточно данных для анализа жирных хвостов")
    return(data.table(
      variable = variable_name,
      distribution = "original",
      kurtosis = NA,
      tail_3sigma = NA,
      tail_2sigma = NA
    ))
  }

  # Стандартизация данных
  standardized <- scale(clean_data)[, 1]

  # Расчет характеристик
  kurtosis_val <- moments::kurtosis(clean_data)

  # Вероятности хвостовых событий
  tail_3sigma <- mean(abs(standardized) > 3)
  tail_2sigma <- mean(abs(standardized) > 2)

  # Анализ различных распределений
  distributions <- c("original", "normal", "t", "stable")

  results <- data.table(
    variable = variable_name,
    distribution = distributions,
    kurtosis = NA_real_,
    tail_3sigma = NA_real_,
    tail_2sigma = NA_real_
  )

  # Оригинальные данные
  results[distribution == "original", `:=`(
    kurtosis = kurtosis_val,
    tail_3sigma = tail_3sigma,
    tail_2sigma = tail_2sigma
  )]

  # Нормальное распределение
  normal_sim <- rnorm(length(clean_data), mean = mean(clean_data), sd = sd(clean_data))
  normal_std <- scale(normal_sim)[, 1]
  results[distribution == "normal", `:=`(
    kurtosis = moments::kurtosis(normal_sim),
    tail_3sigma = mean(abs(normal_std) > 3),
    tail_2sigma = mean(abs(normal_std) > 2)
  )]

  # t-распределение
  tryCatch({
    t_sim <- rt(length(clean_data), df = 5) * sd(clean_data) + mean(clean_data)
    t_std <- scale(t_sim)[, 1]
    results[distribution == "t", `:=`(
      kurtosis = moments::kurtosis(t_sim),
      tail_3sigma = mean(abs(t_std) > 3),
      tail_2sigma = mean(abs(t_std) > 2)
    )]
  }, error = function(e) {})

  # α-stable распределение
  tryCatch({
    alpha_est <- EstimateAlphaParameter(clean_data)
    stable_sim <- stabledist::rstable(length(clean_data), alpha = alpha_est, beta = 0,
                                     gamma = sd(clean_data), delta = mean(clean_data))
    stable_std <- scale(stable_sim)[, 1]
    results[distribution == "stable", `:=`(
      kurtosis = moments::kurtosis(stable_sim),
      tail_3sigma = mean(abs(stable_std) > 3),
      tail_2sigma = mean(abs(stable_std) > 2)
    )]
  }, error = function(e) {})

  return(results)
}

#### 6. Тестирование и валидация ####

#' TestAllFunctions — Тестирование всех функций фреймворка
#' @param test_data данные для тестирования
#' @return список результатов тестирования
TestAllFunctions <- function(test_data = NULL) {

  # Генерация тестовых данных если не предоставлены
  if (is.null(test_data)) {
    set.seed(123)
    test_data <- data.table(
      date = seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "month"),
      volatility = abs(rnorm(48, 0.05, 0.02)),
      ruonia_rate = 6 + cumsum(rnorm(48, 0, 0.5)),
      inflation = 4 + rnorm(48, 0, 2),
      gdp_gap = rnorm(48, 0, 1)
    )
  }

  test_results <- list()

  # Тест HAR модели
  cat("Тестирование HAR модели...\n")
  tryCatch({
    har_data <- PrepareHarData(test_data, "volatility", "date")
    har_model <- EstimateHarModel(har_data)
    har_forecast <- ForecastHarVolatility(har_model, har_data, 6)
    test_results$har_test <- "PASSED"
  }, error = function(e) {
    test_results$har_test <- paste("FAILED:", e$message)
  })

  # Тест HMM
  cat("Тестирование HMM...\n")
  tryCatch({
    hmm_params <- InitializeHmmParameters(3, test_data$ruonia_rate)
    hmm_model <- EstimateHmmParameters(test_data$ruonia_rate, hmm_params)
    hidden_states <- ViterbiDecoding(test_data$ruonia_rate, hmm_model)
    test_results$hmm_test <- "PASSED"
  }, error = function(e) {
    test_results$hmm_test <- paste("FAILED:", e$message)
  })

  # Тест Gradient Boosting
  cat("Тестирование Gradient Boosting...\n")
  tryCatch({
    residuals <- rnorm(100, 0, 1)
    boosting_result <- SelectBestDistribution(residuals)
    alpha_est <- EstimateAlphaParameter(residuals)
    test_results$boosting_test <- "PASSED"
  }, error = function(e) {
    test_results$boosting_test <- paste("FAILED:", e$message)
  })

  # Тест Taleb R²
  cat("Тестирование Taleb R²...\n")
  tryCatch({
    y_actual <- test_data$ruonia_rate
    y_predicted <- y_actual + rnorm(length(y_actual), 0, 0.5)
    r2_result <- CalculateRobustR2(y_actual, y_predicted)

    shapley_data <- test_data[, .(ruonia_rate, inflation, gdp_gap)]
    shapley_result <- ComputeShapleyR2(shapley_data, "ruonia_rate")
    test_results$taleb_test <- "PASSED"
  }, error = function(e) {
    test_results$taleb_test <- paste("FAILED:", e$message)
  })

  # Тест анализа жирных хвостов
  cat("Тестирование анализа жирных хвостов...\n")
  tryCatch({
    tail_analysis <- AnalyzeFatTails(test_data$ruonia_rate, "test_rates")
    test_results$tail_test <- "PASSED"
  }, error = function(e) {
    test_results$tail_test <- paste("FAILED:", e$message)
  })

  # Сохранение результатов
  SaveProjectResults(test_results, "3.Data/test_results.qs")

  cat("\n=== РЕЗУЛЬТАТЫ ТЕСТИРОВАНИЯ ===\n")
  print(test_results)

  return(test_results)
}

#### 7. Rcpp интеграция ####

#' OptimizedMonteCarlo — Оптимизированный Monte Carlo с Rcpp
#' @param n_simulations количество симуляций
#' @param horizon горизонт прогнозирования
#' @param parameters параметры модели
#' @return матрица симуляций
OptimizedMonteCarlo <- function(n_simulations = 1000, horizon = 12, parameters = list()) {

  if (RCPP_AVAILABLE) {
    tryCatch({
      # Использование скомпилированного C++ кода
      cat("Используется Rcpp оптимизированный Monte Carlo\n")

      # Инициализация параметров
      initial_rates <- rep(ifelse(is.null(parameters$initial_rate), 7.0, parameters$initial_rate), n_simulations)

      # Параметры режимов (HMM)
      transition_matrix <- if (is.null(parameters$transition_matrix)) {
        matrix(c(0.8, 0.15, 0.05, 0.2, 0.6, 0.2, 0.1, 0.3, 0.6), nrow = 3, byrow = TRUE)
      } else {
        parameters$transition_matrix
      }

      regime_means <- if (is.null(parameters$regime_means)) c(6.0, 8.0, 12.0) else parameters$regime_means
      regime_sds <- if (is.null(parameters$regime_sds)) c(1.0, 2.0, 1.5) else parameters$regime_sds

      # Shapley веса
      shapley_weights <- if (is.null(parameters$shapley_weights)) {
        c(0.4, 0.3, 0.2, 0.1)
      } else {
        parameters$shapley_weights
      }

      # Макроэкономические переменные (расширяем до горизонта)
      inflation_deviation <- if (is.null(parameters$inflation_deviation)) {
        rep(0.5, horizon)
      } else {
        rep(parameters$inflation_deviation, length.out = horizon)
      }

      gdp_gap <- if (is.null(parameters$gdp_gap)) {
        rep(0.2, horizon)
      } else {
        rep(parameters$gdp_gap, length.out = horizon)
      }

      har_volatility <- if (is.null(parameters$har_volatility)) {
        rep(0.02, horizon)
      } else {
        rep(parameters$har_volatility, length.out = horizon)
      }

      # Стохастические параметры
      alpha_stable <- min(max(ifelse(is.null(parameters$alpha_stable), 1.8, parameters$alpha_stable), 0.5), 2.0)
      gamma_stable <- ifelse(is.null(parameters$gamma_stable), 0.5, parameters$gamma_stable)
      sigma_gaussian <- ifelse(is.null(parameters$sigma_gaussian), 0.3, parameters$sigma_gaussian)
      lambda_poisson <- ifelse(is.null(parameters$lambda_poisson), 0.05, parameters$lambda_poisson)

      # Вызов оптимизированной C++ функции
      result <- FastIntegratedStochasticSimulation(
        n_simulations, horizon, initial_rates,
        transition_matrix, regime_means, regime_sds,
        shapley_weights, inflation_deviation, gdp_gap, har_volatility,
        alpha_stable, gamma_stable, sigma_gaussian, lambda_poisson
      )

      return(result)

    }, error = function(e) {
      cat("Ошибка в Rcpp, переключаюсь на R реализацию:", e$message, "\n")
      RCPP_AVAILABLE <<- FALSE
    })
  }

  # R реализация (если Rcpp недоступен)
  cat("Используется R реализация Monte Carlo\n")

  # Инициализация параметров
  initial_rate <- ifelse(is.null(parameters$initial_rate), 7.0, parameters$initial_rate)
  volatility <- ifelse(is.null(parameters$volatility), 0.02, parameters$volatility)
  drift <- ifelse(is.null(parameters$drift), 0.0, parameters$drift)

  # Генерация симуляций
  simulations <- replicate(n_simulations, {
    rates <- numeric(horizon)
    rates[1] <- initial_rate

    for (t in 2:horizon) {
      shock <- rnorm(1, 0, volatility)
      rates[t] <- rates[t-1] * (1 + drift)
      rates[t] <- max(0.1, min(25, rates[t]))  # Ограничения
    }

    return(rates)
  })

  return(simulations)
}

#' FastHARCalculation — Ускоренный расчет HAR компонент
#' @param volatility_series временной ряд волатильности
#' @return матрица HAR компонент
FastHARCalculation <- function(volatility_series) {

  if (RCPP_AVAILABLE) {
    tryCatch({
      result <- FastHarVolatilityCalculation(volatility_series, DAILY_WINDOW, WEEKLY_WINDOW, MONTHLY_WINDOW)
      return(result)
    }, error = function(e) {
      cat("Ошибка в Rcpp HAR, используется R версия:", e$message, "\n")
    })
  }

  # R реализация
  n <- length(volatility_series)
  daily_vol <- volatility_series
  weekly_vol <- frollmean(volatility_series, WEEKLY_WINDOW, align = "right")
  monthly_vol <- frollmean(volatility_series, MONTHLY_WINDOW, align = "right")
  quarterly_vol <- frollmean(volatility_series, QUARTERLY_WINDOW, align = "right")

  return(cbind(daily_vol, weekly_vol, monthly_vol, quarterly_vol))
}

#' FastHMMForwardBackward — Ускоренный Forward-Backward алгоритм
#' @param observations вектор наблюдений
#' @param hmm_params параметры HMM
#' @return матрица вероятностей состояний
FastHMMForwardBackward <- function(observations, hmm_params) {

  if (RCPP_AVAILABLE) {
    tryCatch({
      result <- FastForwardBackwardAlgorithm(
        observations, hmm_params$transition_matrix,
        hmm_params$emission_params$mean, hmm_params$emission_params$sd,
        hmm_params$initial_probs
      )
      return(result)
    }, error = function(e) {
      cat("Ошибка в Rcpp HMM, используется R версия:", e$message, "\n")
    })
  }

  # R реализация (упрощенная)
  n_obs <- length(observations)
  n_states <- hmm_params$n_states

  # Простое назначение состояний
  states_matrix <- matrix(0, n_obs, n_states)
  for (i in 1:n_obs) {
    if (!is.na(observations[i])) {
      distances <- abs(observations[i] - hmm_params$emission_params$mean)
      best_state <- which.min(distances)
      states_matrix[i, best_state] <- 1
    } else {
      states_matrix[i, ] <- hmm_params$initial_probs
    }
  }

  return(states_matrix)
}

#' CompileRcppOptimizations — Компиляция Rcpp оптимизаций
#' @return логическое значение успешности компиляции
CompileRcppOptimizations <- function() {

  if (!file.exists("2.Cpp/FastMonteCarlo.cpp")) {
    cat("Файл FastMonteCarlo.cpp не найден\n")
    return(FALSE)
  }

  tryCatch({
    library(Rcpp)
    sourceCpp("2.Cpp/FastMonteCarlo.cpp")
    RCPP_AVAILABLE <<- TRUE
    cat("Rcpp оптимизации успешно скомпилированы\n")
      return(TRUE)
})}

#### 4. Красивые функции визуализации ####

#' CreateBeautifulMonteCarloPlot — Создание красивого графика Monte Carlo симуляций
#' @param simulation_results результаты симуляций
#' @param train_data обучающие данные
#' @param forecast_dates даты прогноза
#' @return интерактивный график
CreateBeautifulMonteCarloPlot <- function(simulation_results, train_data, forecast_dates) {

  # Подготовка исторических данных
  historical_data <- train_data[!is.na(ruonia_rate)]
  last_24_months <- max(1, nrow(historical_data) - 72):nrow(historical_data)
  hist_subset <- historical_data[last_24_months]

  # Подготовка прогнозных данных
  forecast_quantiles <- data.table(
    date = forecast_dates,
    q05 = simulation_results$quantiles[, 1],
    q10 = simulation_results$quantiles[, 2],
    q25 = simulation_results$quantiles[, 3],
    median = simulation_results$quantiles[, 4],
    q75 = simulation_results$quantiles[, 5],
    q90 = simulation_results$quantiles[, 6],
    q95 = simulation_results$quantiles[, 7]
  )

  # Выборочные пути для отображения
  n_paths_show <- 200
  selected_paths <- sample(ncol(simulation_results$simulations), n_paths_show)

  # Создание интерактивного графика с белым фоном
  p <- plot_ly()

  # Исторические данные
  p <- p %>% add_lines(
    data = hist_subset,
    x = ~date, y = ~ruonia_rate,
    line = list(color = '#2E86AB', width = 3),
    name = 'Исторические данные',
    hovertemplate = '<b>%{x}</b><br>Ставка: %{y:.2f}%<extra></extra>'
  )

  # Веер неопределенности (90%)
  p <- p %>% add_ribbons(
    data = forecast_quantiles,
    x = ~date, ymin = ~q05, ymax = ~q95,
    fillcolor = 'rgba(255, 182, 193, 0.3)',
    line = list(color = 'transparent'),
    name = '90% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>90% ДИ: %{ymin:.2f}% - %{ymax:.2f}%<extra></extra>'
  )

  # Веер неопределенности (80%)
  p <- p %>% add_ribbons(
    data = forecast_quantiles,
    x = ~date, ymin = ~q10, ymax = ~q90,
    fillcolor = 'rgba(255, 165, 0, 0.4)',
    line = list(color = 'transparent'),
    name = '80% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>80% ДИ: %{ymin:.2f}% - %{ymax:.2f}%<extra></extra>'
  )

  # Центральный веер (50%)
  p <- p %>% add_ribbons(
    data = forecast_quantiles,
    x = ~date, ymin = ~q25, ymax = ~q75,
    fillcolor = 'rgba(255, 215, 0, 0.5)',
    line = list(color = 'transparent'),
    name = '50% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>50% ДИ: %{ymin:.2f}% - %{ymax:.2f}%<extra></extra>'
  )

  # Медианный прогноз
  p <- p %>% add_lines(
    data = forecast_quantiles,
    x = ~date, y = ~median,
    line = list(color = '#D62828', width = 4),
    name = 'Медианный прогноз',
    hovertemplate = '<b>%{x}</b><br>Медиана: %{y:.2f}%<extra></extra>'
  )

  # Добавление выборочных путей
  path_colors <- c('#F94144', '#F3722C', '#F8961E', '#F9C74F', '#90BE6D',
                   '#43AA8B', '#577590', '#277DA1', '#8338EC', '#3A86FF')

  for (i in 1:min(65, length(selected_paths))) {
    path_data <- data.table(
      date = c(last(hist_subset$date), forecast_dates),
      rate = c(last(hist_subset$ruonia_rate), simulation_results$simulations[, selected_paths[i]])
    )

    path_color <- path_colors[((i - 1) %% length(path_colors)) + 1]

    p <- p %>% add_lines(
      data = path_data,
      x = ~date, y = ~rate,
      line = list(color = paste0(path_color, '60'), width = 1.5),
      name = ifelse(i == 1, 'Стохастические пути (выборка)', ''),
      showlegend = (i == 1),
      legendgroup = 'paths',
      hovertemplate = '<b>%{x}</b><br>Ставка: %{y:.2f}%<extra></extra>'
    )
  }

  # Точка начала прогноза
  p <- p %>% add_markers(
    x = last(hist_subset$date), y = last(hist_subset$ruonia_rate),
    marker = list(color = '#D62828', size = 12, symbol = 'diamond'),
    name = 'Начало прогноза',
    hovertemplate = '<b>Начало прогноза</b><br>%{x}<br>Ставка: %{y:.2f}%<extra></extra>'
  )

  # Добавим экстремальные пути к beautiful_plot
  sim_mat <- simulation_results$simulations
  h_len <- nrow(sim_mat)
  last_vals <- sim_mat[h_len, ]
  best_idx <- which.min(last_vals)
  worst_idx <- which.max(last_vals)
  best_path <- sim_mat[, best_idx]
  worst_path <- sim_mat[, worst_idx]
  extreme_df <- data.table(date = forecast_dates,
                           best_path = as.numeric(best_path),
                           worst_path = as.numeric(worst_path))

  p <- p %>%
    add_lines(
      data = extreme_df,
      x = ~date, y = ~best_path,
      line = list(color = '#2ECC71', width = 3.5, shape = 'spline'),
      name = 'Лучший путь (min на горизонте)',
      hovertemplate = '<b>%{x}</b><br>Лучший: %{y:.2f}%<extra></extra>'
    ) %>%
    add_lines(
      data = extreme_df,
      x = ~date, y = ~worst_path,
      line = list(color = '#C0392B', width = 3.5, shape = 'spline'),
      name = 'Худший путь (max на горизонте)',
      hovertemplate = '<b>%{x}</b><br>Худший: %{y:.2f}%<extra></extra>'
    )


  # Настройка макета
  p <- p %>% layout(
    title = list(
      text = '<b>Стохастическое моделирование ключевой ставки ЦБ РФ</b><br><sub>Интегрированный подход: правило Тейлора • HAR волатильность • HMM режимы • α-stable ошибки</sub>',
      font = list(size = 18, family = 'Arial, sans-serif'),
      x = 0.5
    ),
    xaxis = list(
      title = list(text = 'Дата', font = list(size = 14, family = 'Arial, sans-serif')),
      tickfont = list(size = 12, family = 'Arial, sans-serif'),
      gridcolor = '#E5E5E5'
    ),
    yaxis = list(
      title = list(text = 'Ключевая ставка ЦБ РФ (%)', font = list(size = 14, family = 'Arial, sans-serif')),
      tickfont = list(size = 12, family = 'Arial, sans-serif'),
      gridcolor = '#E5E5E5'
    ),
    plot_bgcolor = 'white',
    paper_bgcolor = 'white',
    font = list(family = 'Arial, sans-serif'),
    margin = list(l = 80, r = 80, t = 120, b = 80),
    legend = list(
      orientation = 'h',
      font = list(size = 12, family = 'Arial, sans-serif')
    )
  )

  return(p)
}

#' CreateMacroTrendsPlot — Создание красивого графика макроэкономических трендов
#' @param train_data обучающие данные
#' @param forecast_data прогнозные данные
#' @return интерактивный график
CreateMacroTrendsPlot <- function(train_data, forecast_data) {

  # Подготовка данных
  plot_data <- data.table(
    date = c(train_data$date, forecast_data$date),
    inflation_dev = c(train_data$inflation_deviation, forecast_data$inflation_deviation),
    gdp_gap = c(train_data$gdp_gap, forecast_data$gdp_gap),
    har_vol = c(train_data$rate_volatility, forecast_data$har_volatility),
    type = c(rep("История", nrow(train_data)), rep("Прогноз", nrow(forecast_data)))
  )

  # Создание графика
  p <- plot_ly()

  # Отклонение инфляции
  p <- p %>% add_lines(
    data = plot_data[type == "История"],
    x = ~date, y = ~inflation_dev,
    name = 'Отклонение инфляции (история)',
    line = list(color = '#2E86AB', width = 3),
    hovertemplate = '<b>%{x}</b><br>Отклонение инфляции: %{y:.3f}<extra></extra>'
  ) %>%
  add_lines(
    data = plot_data[type == "Прогноз"],
    x = ~date, y = ~inflation_dev,
    name = 'Отклонение инфляции (прогноз)',
    line = list(color = '#2E86AB', width = 3, dash = 'dash'),
    hovertemplate = '<b>%{x}</b><br>Отклонение инфляции: %{y:.3f}<extra></extra>'
  )

  # Разрыв ВВП
  p <- p %>% add_lines(
    data = plot_data[type == "История"],
    x = ~date, y = ~gdp_gap,
    name = 'Разрыв ВВП (история)',
    line = list(color = '#F8961E', width = 3),
    hovertemplate = '<b>%{x}</b><br>Разрыв ВВП: %{y:.3f}<extra></extra>'
  ) %>%
  add_lines(
    data = plot_data[type == "Прогноз"],
    x = ~date, y = ~gdp_gap,
    name = 'Разрыв ВВП (прогноз)',
    line = list(color = '#F8961E', width = 3, dash = 'dash'),
    hovertemplate = '<b>%{x}</b><br>Разрыв ВВП: %{y:.3f}<extra></extra>'
  )

  # Волатильность ставок
  p <- p %>% add_lines(
    data = plot_data[type == "История"],
    x = ~date, y = ~har_vol,
    name = 'Волатильность ставок (история)',
    line = list(color = '#90BE6D', width = 3),
    hovertemplate = '<b>%{x}</b><br>Волатильность: %{y:.3f}<extra></extra>'
  ) %>%
  add_lines(
    data = plot_data[type == "Прогноз"],
    x = ~date, y = ~har_vol,
    name = 'Волатильность ставок (прогноз)',
    line = list(color = '#90BE6D', width = 3, dash = 'dash'),
    hovertemplate = '<b>%{x}</b><br>Волатильность: %{y:.3f}<extra></extra>'
  )

  # Настройка макета
  p <- p %>% layout(
    title = list(
      text = '<b>Макроэкономические тренды и прогнозы</b><br><sub>Отклонение инфляции • Разрыв ВВП • Волатильность ставок</sub>',
      font = list(size = 18, family = 'Arial, sans-serif'),
      x = 0.5
    ),
    xaxis = list(
      title = list(text = 'Дата', font = list(size = 14, family = 'Arial, sans-serif')),
      tickfont = list(size = 12, family = 'Arial, sans-serif'),
      gridcolor = '#E5E5E5'
    ),
    yaxis = list(
      title = list(text = 'Значение', font = list(size = 14, family = 'Arial, sans-serif')),
      tickfont = list(size = 12, family = 'Arial, sans-serif'),
      gridcolor = '#E5E5E5'
    ),
    plot_bgcolor = 'white',
    paper_bgcolor = 'white',
    font = list(family = 'Arial, sans-serif'),
    margin = list(l = 80, r = 80, t = 120, b = 80),
    hovermode = 'x unified',
    legend = list(
      orientation = 'h',
      x = 0.5,
      y = -0.15,
      font = list(size = 12, family = 'Arial, sans-serif')
    )
  )

  return(p)
}

#### 8. Улучшенное тестирование ####

#' RunComprehensiveTests — Запуск всех тестов включая Rcpp
#' @param test_data данные для тестирования
#' @return подробные результаты тестирования
RunComprehensiveTests <- function(test_data = NULL) {

  cat("=== ЗАПУСК КОМПЛЕКСНОГО ТЕСТИРОВАНИЯ ===\n")

  # Генерация тестовых данных если не предоставлены
  if (is.null(test_data)) {
    set.seed(123)
    test_data <- data.table(
      date = seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "month"),
      volatility = abs(rnorm(48, 0.05, 0.02)),
      ruonia_rate = 6 + cumsum(rnorm(48, 0, 0.5)),
      inflation = 4 + rnorm(48, 0, 2),
      gdp_gap = rnorm(48, 0, 1)
    )
  }

  test_results <- list()

  # Тест Rcpp компиляции
  cat("1. Тестирование Rcpp компиляции...\n")
  test_results$rcpp_compilation <- CompileRcppOptimizations()

  # Тест HAR модели
  cat("2. Тестирование HAR модели...\n")
  tryCatch({
    har_data <- PrepareHarData(test_data, "volatility", "date")
    har_model <- EstimateHarModel(har_data)
    har_forecast <- ForecastHarVolatility(har_model, har_data, 6)

    # Тест ускоренного HAR
    fast_har <- FastHARCalculation(test_data$volatility)

    test_results$har_test <- "PASSED"
    test_results$har_r2 <- har_model$r_squared
    test_results$har_fast_available <- !is.null(fast_har)
  }, error = function(e) {
    test_results$har_test <- paste("FAILED:", e$message)
  })

  # Тест HMM
  cat("3. Тестирование HMM...\n")
  tryCatch({
    hmm_params <- InitializeHmmParameters(3, test_data$ruonia_rate)
    hmm_model <- EstimateHmmParameters(test_data$ruonia_rate, hmm_params)
    hidden_states <- ViterbiDecoding(test_data$ruonia_rate, hmm_model)

    # Тест ускоренного HMM
    fast_hmm <- FastHMMForwardBackward(test_data$ruonia_rate, hmm_model)

    test_results$hmm_test <- "PASSED"
    test_results$hmm_states <- table(hidden_states)
    test_results$hmm_fast_available <- !is.null(fast_hmm)
  }, error = function(e) {
    test_results$hmm_test <- paste("FAILED:", e$message)
  })

  # Тест Gradient Boosting
  cat("4. Тестирование Gradient Boosting...\n")
  tryCatch({
    residuals <- rnorm(100, 0, 1)
    boosting_result <- SelectBestDistribution(residuals)
    alpha_est <- EstimateAlphaParameter(residuals)

    test_results$boosting_test <- "PASSED"
    test_results$best_distribution <- boosting_result$best_distribution
    test_results$alpha_estimate <- alpha_est
  }, error = function(e) {
    test_results$boosting_test <- paste("FAILED:", e$message)
  })

  # Тест Taleb R²
  cat("5. Тестирование Taleb R²...\n")
  tryCatch({
    y_actual <- test_data$ruonia_rate
    y_predicted <- y_actual + rnorm(length(y_actual), 0, 0.5)
    r2_result <- CalculateRobustR2(y_actual, y_predicted)

    shapley_data <- test_data[, .(ruonia_rate, inflation, gdp_gap)]
    shapley_result <- ComputeShapleyR2(shapley_data, "ruonia_rate")

    test_results$taleb_test <- "PASSED"
    test_results$taleb_r2 <- r2_result$r2_taleb_exact
    test_results$shapley_weights <- shapley_result
  }, error = function(e) {
    test_results$taleb_test <- paste("FAILED:", e$message)
  })

  # Тест анализа жирных хвостов
  cat("6. Тестирование анализа жирных хвостов...\n")
  tryCatch({
    tail_analysis <- AnalyzeFatTails(test_data$ruonia_rate, "test_rates")
    test_results$tail_test <- "PASSED"
    test_results$tail_kurtosis <- tail_analysis[distribution == "original", kurtosis]
  }, error = function(e) {
    test_results$tail_test <- paste("FAILED:", e$message)
  })

  # Тест оптимизированного Monte Carlo
  cat("7. Тестирование оптимизированного Monte Carlo...\n")
  tryCatch({
    mc_params <- list(
      initial_rate = 7.0,
      volatility = 0.02,
      alpha_stable = 1.8
    )
    mc_result <- OptimizedMonteCarlo(100, 12, mc_params)

    test_results$monte_carlo_test <- "PASSED"
    test_results$monte_carlo_rcpp <- RCPP_AVAILABLE
    test_results$monte_carlo_paths <- dim(mc_result)
  }, error = function(e) {
    test_results$monte_carlo_test <- paste("FAILED:", e$message)
  })

  # Расчет итоговых метрик
  passed_tests <- sum(sapply(test_results, function(x) {
    if (is.character(x)) return(grepl("PASSED", x))
    if (is.logical(x)) return(x)
    return(FALSE)
  }))

  total_tests <- 7
  test_results$summary <- list(
    passed = passed_tests,
    total = total_tests,
    success_rate = round(passed_tests / total_tests * 100, 1),
    rcpp_available = RCPP_AVAILABLE,
    timestamp = Sys.time()
  )

  # Сохранение результатов
  SaveProjectResults(test_results, "3.Data/comprehensive_test_results.qs")

  cat("\n=== РЕЗУЛЬТАТЫ КОМПЛЕКСНОГО ТЕСТИРОВАНИЯ ===\n")
  cat("Пройдено тестов:", passed_tests, "из", total_tests, "\n")
  cat("Процент успешности:", test_results$summary$success_rate, "%\n")
  cat("Rcpp доступен:", RCPP_AVAILABLE, "\n")

  return(test_results)
}

#### 8. Утилиты ####

#' SaveProjectResults — Сохранение объекта на диск
#' Пакет qs снят с CRAN (архивирован 2026-01-17) и в 9.libraries.r не подключается.
#' Если qs установлен, объект сохраняется через qs::qsave() по указанному пути;
#' иначе — через saveRDS() в тот же путь с расширением .rds.
#' @param object сохраняемый объект
#' @param path путь к файлу .qs
#' @return путь к фактически записанному файлу (невидимо)
SaveProjectResults <- function(object, path) {
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(object, path)
  } else {
    path <- sub("\\.qs$", ".rds", path)
    saveRDS(object, path)
    message("Пакет qs не установлен; объект сохранён через saveRDS(): ", path)
  }
  invisible(path)
}

#' SaveResults — Сохранение результатов анализа
#' @param results список результатов
#' @param filename имя файла для сохранения
SaveResults <- function(results, filename = "analysis_results.qs") {
  filepath <- SaveProjectResults(results, file.path("3.Data", filename))
  cat("Результаты сохранены в:", filepath, "\n")
}

#' LoadResults — Загрузка результатов анализа (.qs через qs, иначе .rds)
#' @param filename имя файла для загрузки
#' @return загруженные результаты
LoadResults <- function(filename = "analysis_results.qs") {
  filepath <- file.path("3.Data", filename)
  rds_path <- sub("\\.qs$", ".rds", filepath)
  if (file.exists(filepath) && requireNamespace("qs", quietly = TRUE)) {
    return(qs::qread(filepath))
  } else if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  } else {
    warning("Файл не найден:", filepath)
    return(NULL)
  }
}

#' InterpolateDataWithSVD — Универсальная функция интерполяции данных
#' @param input_data data.table с данными для экстраполяции
#' @param extrapolation_vars вектор переменных для экстраполяции
#' @param time_var название временной переменной
#' @param horizon количество периодов для прогноза
#' @return data.table с экстраполированными данными
InterpolateDataWithSVD <- function(input_data, extrapolation_vars, time_var = "date", horizon = 12) {

  # Проверка входных данных
  if (!all(extrapolation_vars %in% names(input_data))) {
    missing_vars <- setdiff(extrapolation_vars, names(input_data))
    stop("Переменные не найдены в данных: ", paste(missing_vars, collapse = ", "))
  }

  if (!time_var %in% names(input_data)) {
    stop("Временная переменная '", time_var, "' не найдена в данных")
  }

  # Подготовка данных для SVD прогнозирования
  prediction_data <- input_data[, c(time_var, extrapolation_vars), with = FALSE]

  # Применение SVD-VAR экстраполяции
  svd_result <- PredictTimeSeriesWithSVD(
    input_dataset = prediction_data,
    prediction_steps = horizon,
    temporal_var = time_var
  )

  return(svd_result$forecast)
}

#' SmartExtrapolate — Умная экстраполяция с выбором метода
#' @param time_series временной ряд
#' @param horizon горизонт прогнозирования
#' @param method метод экстраполяции ("auto", "linear", "ar", "ets")
#' @return вектор экстраполированных значений
SmartExtrapolate <- function(time_series, horizon = 12, method = "auto") {

  # Удаление NA
  clean_series <- time_series[!is.na(time_series)]

  if (length(clean_series) < 3) {
    # Для коротких рядов - константная экстраполяция
    return(rep(mean(clean_series, na.rm = TRUE), horizon))
  }

  if (method == "auto") {
    # Автоматический выбор лучшего метода
    methods_to_try <- c("ar", "ets", "linear")
    best_method <- "linear"
    best_aic <- Inf

    for (test_method in methods_to_try) {
      tryCatch({
        test_result <- SmartExtrapolate(clean_series, min(3, horizon), test_method)
        if (test_method == "ar") {
          ar_model <- ar(clean_series, method = "mle")
          aic <- ar_model$aic
        } else if (test_method == "ets") {
          ets_model <- forecast::ets(ts(clean_series))
          aic <- ets_model$aic
        } else {
          aic <- sum((clean_series - fitted(lm(clean_series ~ seq_along(clean_series))))^2)
        }

        if (aic < best_aic) {
          best_aic <- aic
          best_method <- test_method
        }
      }, error = function(e) {})
    }

    method <- best_method
  }

  # Применение выбранного метода
  switch(method,
    "ar" = {
      ar_model <- ar(clean_series, method = "mle", order.max = min(3, length(clean_series) - 1))
      ar_forecast <- predict(ar_model, n.ahead = horizon)
      return(as.numeric(ar_forecast$pred))
    },
    "ets" = {
      ets_model <- forecast::ets(ts(clean_series))
      ets_forecast <- forecast::forecast(ets_model, h = horizon)
      return(as.numeric(ets_forecast$mean))
    },
    "linear" = {
      time_indices <- seq_along(clean_series)
      linear_model <- lm(clean_series ~ time_indices)
      future_indices <- (length(clean_series) + 1):(length(clean_series) + horizon)
      linear_forecast <- predict(linear_model, newdata = data.frame(time_indices = future_indices))
      return(as.numeric(linear_forecast))
    },
    {
      # По умолчанию - среднее значение
      return(rep(mean(clean_series), horizon))
    }
  )
}

cat("Все функции стохастического моделирования с Rcpp оптимизациями загружены успешно!\n")
if (RCPP_AVAILABLE) {
  cat("Rcpp оптимизации активны и готовы к использованию!\n")
} else {
  cat("Используются R реализации (Rcpp недоступен)\n")
}


#' CalibrateStochasticDistributions — Калибровка параметров распределений по МНП для ВСЕХ переменных
#' @param data обучающие данные
#' @return список калиброванных параметров для каждой переменной
CalibrateStochasticDistributions <- function(data) {

# 1. Волатильность как случайный процесс (стохастическая волатильность)
volatility_series <- data$rate_volatility

# Подгонка параметров для стохастической волатильности (модель Хестона упрощенная)
# dV_t = κ(θ - V_t)dt + σ_v√V_t dW_t
# Используем абсолютные значения волатильности, а не логарифмы
vol_changes <- diff(volatility_series)

# МНП для параметров стохастической волатильности
estimate_heston_params <- function(params) {
  kappa <- params[1]  # Скорость возврата к среднему
  theta <- params[2]  # Долгосрочное среднее волатильности
  sigma_v <- params[3] # Волатильность волатильности

  if (kappa <= 0 || sigma_v <= 0 || theta <= 0) return(1e10)

  # Лог-правдоподобие для дискретной версии модели с абсолютными значениями
  residuals <- vol_changes - kappa * (theta - volatility_series[-length(volatility_series)])
  # Волатильность инноваций пропорциональна квадратному корню текущей волатильности
  vol_scale <- sqrt(volatility_series[-length(volatility_series)])
  scaled_residuals <- residuals / (sigma_v * vol_scale)

  ll <- -0.5 * sum(log(2 * pi * (sigma_v * vol_scale)^2)) - 0.5 * sum(scaled_residuals^2)
  return(-ll)  # Минимизируем отрицательное лог-правдоподобие
}

heston_result <- optim(
  par = c(0.5, mean(volatility_series), sd(vol_changes)),
  fn = estimate_heston_params,
  method = "L-BFGS-B",
  lower = c(0.01, 0.001, 0.001),
  upper = c(5, 1, 1)
)

volatility_params <- list(
  kappa = heston_result$par[1],
  theta = heston_result$par[2],
  sigma_v = heston_result$par[3],
  mean_reversion_level = heston_result$par[2]
)

# 2. Калибровка параметров для КАЖДОЙ макропеременной
# Функция оценки параметров жирохвостого распределения для каждой переменной
EstimateHeavyTailParams <- function(variable_data, variable_name) {

  clean_data <- variable_data[!is.na(variable_data)]

  if (length(clean_data) < 10) {
    return(list(
      distribution = "normal",
      params = list(mean = mean(clean_data), sd = sd(clean_data)),
      log_likelihood = -Inf
    ))
  }

  # Пробуем разные распределения и выбираем лучшее по AIC
  distributions <- list(
    t_student = function(params) {
      if (params[2] <= 0 || params[3] <= 2) return(1e10)
      -sum(dt((clean_data - params[1])/params[2], df = params[3], log = TRUE) - log(params[2]))
    },
    pareto = function(params) {
      if (params[1] <= 0 || params[2] <= 0) return(1e10)
      # Используем исходные данные без сдвига для калибровки
      # params[1] = alpha (shape), params[2] = scale
      if (any(clean_data <= params[2])) return(1e10)  # Все данные должны быть больше scale
      tryCatch({
        -sum(extraDistr::dpareto(clean_data, a = params[1], b = params[2], log = TRUE))
      }, error = function(e) 1e10)
    }
  )

  # Начальные значения для каждого распределения
  initial_params <- list(
    t_student = c(mean(clean_data), sd(clean_data), 5),
    pareto = c(3.0, quantile(clean_data, 0.1))  # Более реалистичные начальные значения
  )

  # Границы параметров
  lower_bounds <- list(
    t_student = c(-Inf, 0.01, 2.1),
    pareto = c(0.1, 0.01)
  )

  upper_bounds <- list(
    t_student = c(Inf, Inf, 50),
    pareto = c(10, Inf)
  )

  best_result <- NULL
  best_aic <- Inf
  best_distribution <- "normal"

  for (dist_name in names(distributions)) {
    tryCatch({
      result <- optim(
        par = initial_params[[dist_name]],
        fn = distributions[[dist_name]],
        method = "L-BFGS-B",
        lower = lower_bounds[[dist_name]],
        upper = upper_bounds[[dist_name]]
      )

      # Расчет AIC
      k <- length(initial_params[[dist_name]])  # количество параметров
      aic <- 2 * result$value + 2 * k

      if (aic < best_aic && result$convergence == 0) {
        best_aic <- aic
        best_result <- result
        best_distribution <- dist_name
      }
    }, error = function(e) {})
  }

  # Возвращаем лучший результат
  return(list(
    distribution = best_distribution,
    params = as.list(best_result$par),
    log_likelihood = -best_result$value,
    aic = best_aic,
    variable_name = variable_name
  ))
}

# Калибровка для каждой макропеременной
# Оставляем только те переменные, которые реально используются в симуляции
macro_variables <- c("cpi_inflation", "gdp_gap")  # Только инфляция и ВВП используются в симуляции

macro_distributions <- list()
for (var_name in macro_variables) {
  if (var_name %in% names(data)) {
    macro_distributions[[var_name]] <- EstimateHeavyTailParams(data[[var_name]], var_name)
  }
}

# 3. Остатки по степенному закону (fat-tailed errors)
# Используем простую линейную модель для получения остатков
basic_model <- lm(ruonia_rate ~ inflation_deviation + gdp_gap, data = data)
residuals <- residuals(basic_model)

# МНП для α-stable распределения остатков
estimate_stable_params <- function(params) {
  alpha <- params[1]   # Параметр устойчивости (tail index)
  beta <- params[2]    # Параметр асимметрии
  gamma <- params[3]   # Масштабный параметр
  delta <- params[4]   # Параметр сдвига

  # Ограничения параметров
  if (alpha <= 0.1 || alpha > 2 || abs(beta) > 1 || gamma <= 0) return(1e10)

  tryCatch({
    ll <- sum(stabledist::dstable(residuals, alpha = alpha, beta = beta,
                                  gamma = gamma, delta = delta, log = TRUE))
    return(-ll)
  }, error = function(e) return(1e10))
}

stable_result <- optim(
  par = c(1.8, 0, sd(residuals), mean(residuals)),
  fn = estimate_stable_params,
  method = "L-BFGS-B",
  lower = c(0.5, -0.99, 0.01, -5),
  upper = c(2.0, 0.99, 5, 5)
)

stable_params <- list(
  alpha = stable_result$par[1],
  beta = stable_result$par[2],
  gamma = stable_result$par[3],
  delta = stable_result$par[4]
)

# 4. Анализ хвостов инфляции и ВВП
inflation_tail <- AnalyzeFatTails(data$cpi_inflation, "Inflation")
gdp_tail <- AnalyzeFatTails(data$gdp_gap, "GDP_Gap")

# Сохранение всех параметров
calibrated_params <- list(
  volatility = volatility_params,
  stable_errors = stable_params,
  macro_distributions = macro_distributions,  # Новое: параметры для каждой макропеременной
  tail_analysis = list(
    inflation = inflation_tail,
    gdp = gdp_tail
  ),
  model_residuals = residuals,
  log_likelihood = -stable_result$value,
  calibration_summary = data.table(
    variable = c("volatility", names(macro_distributions), "stable_errors"),
    distribution = c("heston", sapply(macro_distributions, function(x) x$distribution), "stable"),
    aic = c(NA, sapply(macro_distributions, function(x) x$aic), NA)
  )
)

return(calibrated_params)
}

#### 1.1. Функция калибровки параметров правила Тейлора методом максимального правдоподобия ####

CalibrateTaylorRuleParameters <- function(historical_data) {
  # Функция для оценки параметров правила Тейлора на исторических данных
  # Используем метод максимального правдоподобия для калибровки весов

  cat("Калибровка параметров правила Тейлора на исторических данных...\n")

  # Подготовка данных для калибровки
  calib_data <- copy(historical_data) %>%
    .[!is.na(ruonia_rate) & !is.na(inflation_deviation) & !is.na(gdp_gap)] %>%
    .[, lag_rate := lag(ruonia_rate, 1)] %>%
    .[, rate_change := ruonia_rate - lag_rate] %>%
    na.omit()

  cat("   Размер выборки для калибровки:", nrow(calib_data), "наблюдений\n")

  # Функция правдоподобия для правила Тейлора
  # i_t = ρ*i_{t-1} + (1-ρ)*[r* + φ_π*(π_t-π*) + φ_y*y_t] + ε_t
  # где ε_t ~ N(0, σ²)

  log_likelihood_taylor <- function(params, data) {
    # params: [r_star, phi_pi, phi_y, rho, sigma]
    r_star <- params[1]
    phi_pi <- params[2]
    phi_y <- params[3]
    rho <- params[4]
    sigma <- params[5]

    # Ограничения параметров
    if (phi_pi < 0.1 || phi_pi > 5 || phi_y < -2 || phi_y > 2 ||
        rho < 0.1 || rho > 0.99 || sigma < 0.01) {
      return(-Inf)
    }

    # Расчет теоретических ставок по правилу Тейлора
    taylor_rates <- r_star + phi_pi * data$inflation_deviation + phi_y * data$gdp_gap

    # Сглаженные ставки
    smoothed_rates <- rho * data$lag_rate + (1 - rho) * taylor_rates

    # Ошибки модели
    errors <- data$ruonia_rate - smoothed_rates

    # Логарифм правдоподобия (нормальное распределение ошибок)
    n <- length(errors)
    log_lik <- -n/2 * log(2*pi) - n/2 * log(sigma^2) - sum(errors^2) / (2*sigma^2)

    return(log_lik)
  }

  # Начальные значения параметров
  initial_params <- c(
    r_star = mean(calib_data$ruonia_rate, na.rm = TRUE),  # Нейтральная ставка
    phi_pi = 1.5,    # Вес инфляции
    phi_y = 0.5,     # Вес выпуска
    rho = 0.8,       # Параметр сглаживания
    sigma = sd(calib_data$rate_change, na.rm = TRUE)  # Стандартное отклонение ошибок
  )

  cat("   Начальные параметры:\n")
  cat("     r* (нейтральная ставка):", round(initial_params[1], 3), "\n")
  cat("     φ_π (вес инфляции):", round(initial_params[2], 3), "\n")
  cat("     φ_y (вес выпуска):", round(initial_params[3], 3), "\n")
  cat("     ρ (сглаживание):", round(initial_params[4], 3), "\n")
  cat("     σ (волатильность ошибок):", round(initial_params[5], 3), "\n")

  # Оптимизация методом максимального правдоподобия
  tryCatch({
    optimization_result <- optim(
      par = initial_params,
      fn = function(params) -log_likelihood_taylor(params, calib_data),
      method = "L-BFGS-B",
      lower = c(0, 0.1, -2, 0.1, 0.01),
      upper = c(20, 5, 2, 0.99, 10),
      control = list(maxit = 1000, trace = 0)
    )

    if (optimization_result$convergence == 0) {
      calibrated_params <- optimization_result$par
      names(calibrated_params) <- c("r_star", "phi_pi", "phi_y", "rho", "sigma")

      cat("   Калибровка успешна! Оптимальные параметры:\n")
      cat("     r* (нейтральная ставка):", round(calibrated_params[1], 3), "\n")
      cat("     φ_π (вес инфляции):", round(calibrated_params[2], 3), "\n")
      cat("     φ_y (вес выпуска):", round(calibrated_params[3], 3), "\n")
      cat("     ρ (сглаживание):", round(calibrated_params[4], 3), "\n")
      cat("     σ (волатильность ошибок):", round(calibrated_params[5], 3), "\n")
      cat("     Логарифм правдоподобия:", round(-optimization_result$value, 3), "\n")

      # Диагностика качества калибровки
      fitted_rates <- calibrated_params[1] +
        calibrated_params[2] * calib_data$inflation_deviation +
        calibrated_params[3] * calib_data$gdp_gap

      smoothed_fitted <- calibrated_params[4] * calib_data$lag_rate +
        (1 - calibrated_params[4]) * fitted_rates

      residuals_calib <- calib_data$ruonia_rate - smoothed_fitted
      r_squared_calib <- 1 - sum(residuals_calib^2) / sum((calib_data$ruonia_rate - mean(calib_data$ruonia_rate))^2)

      cat("   R² калиброванной модели:", round(r_squared_calib, 4), "\n")
      cat("   RMSE:", round(sqrt(mean(residuals_calib^2)), 4), "\n")

      return(list(
        parameters = calibrated_params,
        log_likelihood = -optimization_result$value,
        r_squared = r_squared_calib,
        rmse = sqrt(mean(residuals_calib^2)),
        fitted_values = smoothed_fitted,
        residuals = residuals_calib,
        convergence = TRUE
      ))

    } else {
      cat("   Ошибка оптимизации. Используем начальные параметры.\n")
      return(list(
        parameters = initial_params,
        log_likelihood = -log_likelihood_taylor(initial_params, calib_data),
        r_squared = NA,
        rmse = NA,
        fitted_values = NA,
        residuals = NA,
        convergence = FALSE
      ))
    }

  }, error = function(e) {
    cat("   Ошибка при калибровке:", e$message, "\n")
    cat("   Используем начальные параметры.\n")
    return(list(
      parameters = initial_params,
      log_likelihood = -log_likelihood_taylor(initial_params, calib_data),
      r_squared = NA,
      rmse = NA,
      fitted_values = NA,
      residuals = NA,
      convergence = FALSE
    ))
  })
}