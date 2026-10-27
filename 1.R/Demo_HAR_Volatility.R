#' Demo_HAR_Volatility.R — Демонстрация HAR модели волатильности ставок
#' Использует точно такой же подход как в Main.R
#' Методология: HAR (Heterogeneous Autoregressive) модель для прогнозирования
#' волатильности процентных ставок с интеграцией в стохастическую модель

#### 0. Подключение библиотек ####
source("1.R/8.SourceAll.R")

#### 1. Константы и параметры ####
set.seed(123)
HORIZON <- 19
TARGET_INFLATION <- 4.0
NEUTRAL_RATE <- 7.0

#### 2. Загрузка данных ####
prep_har_data <- function() {
  # Загружаем данные точно как в Main.R
  inflation <- read_xlsx("3.Data/Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx") %>%
    as.data.table() %>%
    .[,`Дата` := as.Date(paste0("01.", `Дата`), format = "%d.%m.%Y")] %>%
    na.omit() %>%
    setnames("Инфляция, % г/г","cpi_inflation") %>%
    setnames("Ключевая ставка, % годовых", "key_rate") %>%
    .[,-c("Цель по инфляции")]
  
  # Межбанковская реальная ставка
  st_ruonia <- read_xlsx("3.Data/RC_F11_01_2010_T04_07_2025.xlsx") %>%
    as.data.table() %>%
    .[,c("DT","ruo")] %>%
    .[, ym := format(as.Date(DT), "%Y-%m")] %>%
    .[, .(`ruonia_rate` = mean(ruo, na.rm = TRUE)), by = .(ym)] %>%
    .[,`Дата` := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
    .[,-c("ym")]
  
  # Объединение данных
  combined_data <- merge(inflation, st_ruonia, by = "Дата", all.x = TRUE) %>%
    setnames("Дата", "date") %>%
    na.omit() %>%
    setorder(date)
  
  return(combined_data)
}

input_data <- prep_har_data()

#### 3. Построение волатильности ставок точно как в Main.R ####
CreateRateVolatility <- function(data) {
  # Создание волатильности ставок (realized volatility) точно как в Main.R
  data_with_vol <- copy(data) %>%
    .[, rate_volatility := frollapply(ruonia_rate, 6, function(x) {
      if (sum(!is.na(x)) < 3) return(NA)
      sd(x, na.rm = TRUE)
    }, align = "right")] %>%
    .[, rate_volatility := ifelse(is.na(rate_volatility) | rate_volatility == 0,
                                  mean(rate_volatility, na.rm = TRUE), rate_volatility)] %>%
    na.omit()
  
  return(data_with_vol)
}

# Создание волатильности
data_with_volatility <- CreateRateVolatility(input_data)

#### 4. HAR модель точно как в Main.R ####
FitHARModel <- function(train_data, horizon = HORIZON) {
  
  # Подготовка данных для HAR модели точно как в Main.R
  volatility_series <- train_data$rate_volatility[!is.na(train_data$rate_volatility)]
  dates <- seq.Date(as.Date("2016-01-01"), by = "month", length.out = length(volatility_series)-1)
  volatility <- as.xts(volatility_series[-length(volatility_series)], order.by = dates)
  
  # Параметры HAR точно как в Main.R
  periods <- c(5,22)
  
  # Прогноз на horizon месяцев (рекурсивно) точно как в Main.R
  forecast_vals <- numeric(horizon)
  
  cat("Выполняю HAR прогнозирование на", horizon, "периодов...\n")
  
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
    
    if (i %% 5 == 0) cat("  Завершен прогноз для периода", i, "\n")
  }
  
  return(list(
    forecast_volatility = forecast_vals,
    historical_volatility = volatility_series,
    har_parameters = list(periods = periods, horizon = horizon),
    model_info = list(
      type = "HAR",
      estimation_method = "rCov", 
      transform = "none",
      input_type = "RM"
    )
  ))
}

# Обучение HAR модели
har_results <- FitHARModel(data_with_volatility, HORIZON)

#### 5. Интеграция со стохастической компонентой как в Main.R ####
IntegrateHARWithStochastic <- function(har_forecast, n_simulations = 1000) {
  
  # Параметры стохастической волатильности (Хестон процесс)
  vol_params <- list(
    kappa = 0.5,    # Скорость возврата к среднему
    theta = mean(data_with_volatility$rate_volatility, na.rm = TRUE),  # Долгосрочное среднее
    sigma_v = 0.3   # Волатильность волатильности
  )
  
  horizon <- length(har_forecast)
  
  # Симуляция стохастических путей
  stochastic_paths <- matrix(0, nrow = horizon, ncol = n_simulations)
  
  for (sim in 1:n_simulations) {
    current_vol <- tail(data_with_volatility$rate_volatility, 1)
    
    for (t in 1:horizon) {
      # HAR прогноз для периода t
      har_vol_forecast <- har_forecast[t]
      
      # Стохастическая компонента волатильности
      vol_innovation <- rnorm(1, 0, vol_params$sigma_v * sqrt(current_vol))
      stochastic_vol_component <- vol_params$kappa * (vol_params$theta - current_vol) + vol_innovation
      
      # Интеграция точно как в Main.R: 60% HAR + 40% стохастическая эволюция
      current_vol <- 0.6 * har_vol_forecast + 0.4 * (current_vol + stochastic_vol_component)
      
      # Ограничения на положительность
      if (!is.finite(current_vol) || is.na(current_vol) || current_vol <= 0) {
        current_vol <- 1e-4
      }
      
      stochastic_paths[t, sim] <- current_vol
    }
  }
  
  # Статистики интегрированного прогноза
  integrated_forecast <- list(
    har_forecast = har_forecast,
    stochastic_paths = stochastic_paths,
    mean_forecast = apply(stochastic_paths, 1, mean),
    quantiles = apply(stochastic_paths, 1, function(x) {
      quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
    }),
    vol_params = vol_params,
    integration_formula = "vol_t = 0.6 * HAR_forecast_t + 0.4 * stochastic_evolution_t"
  )
  
  return(integrated_forecast)
}

# Интеграция HAR с стохастической компонентой
integrated_har <- IntegrateHARWithStochastic(har_results$forecast_volatility, 1000)

#### 6. Создание прогнозных дат ####
last_date <- max(data_with_volatility$date)
forecast_dates <- seq.Date(last_date + months(1), last_date + months(HORIZON), by = "month")

#### 7. Визуализация результатов HAR модели ####
CreateHARVisualization <- function(historical_data, har_results, integrated_forecast, forecast_dates) {
  
  # Подготовка данных для визуализации
  hist_subset <- tail(historical_data, 36)  # Последние 3 года
  
  # Создание красивого графика
  p <- plot_ly()
  
  # Исторические данные волатильности
  p <- p %>% add_lines(
    data = hist_subset,
    x = ~date, y = ~rate_volatility,
    line = list(color = '#2E86AB', width = 3.5, shape = 'spline'),
    name = 'Историческая волатильность',
    hovertemplate = '<b>%{x}</b><br>Волатильность: %{y:.4f}<extra></extra>'
  )
  
  # HAR прогноз (детерминистический)
  p <- p %>% add_lines(
    x = forecast_dates, y = har_results$forecast_volatility,
    line = list(color = '#E74C3C', width = 3, dash = 'dash'),
    name = 'HAR прогноз (детерминистический)',
    hovertemplate = '<b>%{x}</b><br>HAR прогноз: %{y:.4f}<extra></extra>'
  )
  
  # Интегрированный прогноз (медиана)
  p <- p %>% add_lines(
    x = forecast_dates, y = integrated_forecast$mean_forecast,
    line = list(color = '#F7931E', width = 4, shape = 'spline'),
    name = 'Интегрированный прогноз (HAR + стохастический)',
    hovertemplate = '<b>%{x}</b><br>Интегрированный: %{y:.4f}<extra></extra>'
  )
  
  # Доверительные интервалы
  quantiles_matrix <- t(integrated_forecast$quantiles)
  
  # 90% доверительный интервал
  p <- p %>% add_ribbons(
    x = forecast_dates,
    ymin = quantiles_matrix[, 1], ymax = quantiles_matrix[, 7],
    fillcolor = 'rgba(255, 140, 0, 0.15)',
    line = list(color = 'rgba(255, 140, 0, 0.3)', width = 1),
    name = '90% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>90% ДИ: %{ymin:.4f} - %{ymax:.4f}<extra></extra>'
  )
  
  # 50% доверительный интервал
  p <- p %>% add_ribbons(
    x = forecast_dates,
    ymin = quantiles_matrix[, 3], ymax = quantiles_matrix[, 5],
    fillcolor = 'rgba(255, 180, 0, 0.25)',
    line = list(color = 'rgba(255, 180, 0, 0.4)', width = 1),
    name = '50% доверительный интервал',
    hovertemplate = '<b>%{x}</b><br>50% ДИ: %{ymin:.4f} - %{ymax:.4f}<extra></extra>'
  )
  
  # Точка начала прогноза
  p <- p %>% add_markers(
    x = last_date, y = tail(historical_data$rate_volatility, 1),
    marker = list(
      size = 15,
      color = '#E74C3C',
      symbol = 'diamond',
      line = list(color = '#E74C3C', width = 3)
    ),
    name = 'Начало прогноза',
    hovertemplate = '<b>Начало прогноза</b><br>%{x}<br>Волатильность: %{y:.4f}<extra></extra>'
  )
  
  # Оформление
  p <- p %>% layout(
    title = list(
      text = '<b>HAR модель волатильности процентных ставок ЦБ РФ</b><br><sub>Методология: HAR прогноз + стохастическая интеграция (60% + 40%)</sub>',
      font = list(color = '#2C3E50', size = 20, family = 'Arial'),
      x = 0.5,
      y = 0.98
    ),
    xaxis = list(
      title = list(text = '<b>Временной горизонт</b>', font = list(color = '#2C3E50', size = 16)),
      color = '#2C3E50',
      gridcolor = 'rgba(200,200,200,0.5)',
      tickfont = list(color = '#2C3E50', size = 12)
    ),
    yaxis = list(
      title = list(text = '<b>Волатильность ставок</b>', font = list(color = '#2C3E50', size = 16)),
      color = '#2C3E50',
      gridcolor = 'rgba(200,200,200,0.5)',
      tickfont = list(color = '#2C3E50', size = 12)
    ),
    plot_bgcolor = 'rgba(255,255,255,1)',
    paper_bgcolor = 'rgba(255,255,255,1)',
    font = list(color = '#2C3E50', family = 'Arial'),
    margin = list(l = 100, r = 100, t = 140, b = 100),
    showlegend = TRUE
  )
  
  return(p)
}

# Создание визуализации
har_plot <- CreateHARVisualization(data_with_volatility, har_results, integrated_har, forecast_dates)

#### 8. Диагностические метрики HAR модели ####
har_diagnostics <- data.table(
  Period = 1:HORIZON,
  Date = forecast_dates,
  HAR_Forecast = har_results$forecast_volatility,
  Integrated_Mean = integrated_har$mean_forecast,
  Q05 = t(integrated_har$quantiles)[, 1],
  Q25 = t(integrated_har$quantiles)[, 3],
  Q50 = t(integrated_har$quantiles)[, 4],
  Q75 = t(integrated_har$quantiles)[, 5],
  Q95 = t(integrated_har$quantiles)[, 7],
  Uncertainty_Range = t(integrated_har$quantiles)[, 5] - t(integrated_har$quantiles)[, 3]
)

# Ключевые статистики
current_vol <- tail(data_with_volatility$rate_volatility, 1)
forecast_12m <- har_diagnostics$Integrated_Mean[12]
volatility_trend <- ifelse(forecast_12m > current_vol, "Возрастающий", "Убывающий")

cat("\n=== РЕЗУЛЬТАТЫ HAR МОДЕЛИ ВОЛАТИЛЬНОСТИ ===\n")
cat("Текущая волатильность:", round(current_vol, 4), "\n")
cat("Прогноз через 12 месяцев:", round(forecast_12m, 4), "\n")
cat("Тренд волатильности:", volatility_trend, "\n")
cat("Среднее изменение:", round((forecast_12m/current_vol - 1)*100, 2), "%\n")

#### 9. Сохранение результатов ####
har_model_results <- list(
  historical_data = data_with_volatility,
  har_forecast = har_results,
  integrated_forecast = integrated_har,
  diagnostics = har_diagnostics,
  forecast_dates = forecast_dates,
  methodology = list(
    approach = "HAR (Heterogeneous Autoregressive) модель",
    periods = c(5, 22),
    integration_formula = "vol_t = 0.6 * HAR_forecast_t + 0.4 * stochastic_evolution_t",
    stochastic_process = "Хестон процесс для волатильности волатильности",
    implementation = "Точная копия подхода из Main.R"
  ),
  model_performance = list(
    current_volatility = current_vol,
    forecast_12m = forecast_12m,
    volatility_trend = volatility_trend,
    mean_uncertainty = mean(har_diagnostics$Uncertainty_Range)
  ),
  timestamp = Sys.time()
)

# Сохранение в файл
qsave(har_model_results, "3.Data/har_volatility_model.qs")

#### 10. Отображение результатов ####
print(har_plot)

cat("\n=== HAR ДИАГНОСТИКА ===\n")
print(har_diagnostics[, .(Period, Date, HAR_Forecast, Integrated_Mean, Q25, Q75, Uncertainty_Range)])

cat("\nHAR модель успешно обучена и применена точно по методологии Main.R\n")
cat("Результаты сохранены в: 3.Data/har_volatility_model.qs\n")
