# Demo_Volatility_Analysis.R
# Enhanced HAR Model with Vol-of-Vol and Stochastic Integration
# Author: Nikita Burakov
# Based on: Corsi (2009), Clements & Preve (2018), Lu (2022)

#### 0. Load Libraries ####
library(data.table)
library(plotly)
library(tidyverse)
library(zoo)
library(lubridate)
library(rugarch)
library(xts)
library(highfrequency)

#### 1. Generate Realistic Financial Data ####

set.seed(42)
n_days <- 2500  # ~6 years of daily data
start_date <- as.Date("2015-01-01")
dates <- seq(start_date, by = "day", length.out = n_days)

# Generate realistic financial data with volatility events
generate_realistic_financial_data <- function() {

  # Parameters
  mu <- 0.0003
  base_vol <- 0.012

  # Create event calendar (crises, shocks, recoveries)
  crisis_periods <- list(
    list(start = 200, end = 240, intensity = 3.5/10),
    list(start = 800, end = 830, intensity = 2.8/10),
    list(start = 1500, end = 1520, intensity = 4.2/10),
    list(start = 2000, end = 2150, intensity = 3.0/10)
  )

  # Enhanced volatility regimes
  returns <- numeric(n_days)
  vol_series <- numeric(n_days)
  shock_indicator <- numeric(n_days)

  # Initial conditions
  vol_series[1] <- base_vol
  returns[1] <- rnorm(1, mu, vol_series[1])

  # GARCH parameters
  omega <- 0.0000005  # constant
  alpha <- 0.12/5      # ARCH effect
  beta <- 0.85/5       # GARCH effect

  for(t in 2:n_days) {
    crisis_multiplier <- 1.0
    current_shock <- 0

    for(crisis in crisis_periods) {
      if(t >= crisis$start && t <= crisis$end) {
        crisis_phase <- (t - crisis$start) / (crisis$end - crisis$start)
        if(crisis_phase <= 0.3) {
          intensity <- 1 + (crisis$intensity - 1) * (crisis_phase / 0.3)^2
        } else if(crisis_phase <= 0.7) {
          intensity <- crisis$intensity * (1 + 0.3 * sin(10 * pi * crisis_phase))
        } else {
          intensity <- crisis$intensity * (1 - 0.7 * ((crisis_phase - 0.7) / 0.3)^0.5)
        }
        crisis_multiplier <- max(crisis_multiplier, intensity)
        current_shock <- crisis$intensity
      }
    }

    # Volatile days (random shocks)
    if(runif(1) < 0.02) {
      vol_shock <- rnorm(1, 0, 0.02)
      crisis_multiplier <- crisis_multiplier * (1 + abs(vol_shock) * 5)
      current_shock <- max(current_shock, abs(vol_shock) * 10)
    }

    # Weekly cycles (Monday/Friday effects)
    day_of_week <- (t %% 7) + 1
    week_effect <- ifelse(day_of_week == 2, 1.1,      # Monday
                   ifelse(day_of_week == 6, 1.15, 1.0)) # Friday

    # Monthly effects (end of month, quarter)
    month_day <- (t %% 22) + 1
    month_effect <- ifelse(month_day <= 2 || month_day >= 21, 1.05, 1.0)

    # GARCH with crisis multipliers
    vol_series[t] <- sqrt(omega + alpha * returns[t-1]^2 + beta * vol_series[t-1]^2) *
                     crisis_multiplier * week_effect * month_effect

    # Leverage effect (asymmetry)
    leverage_effect <- ifelse(returns[t-1] < 0, 1.2, 0.9)
    vol_series[t] <- vol_series[t] * leverage_effect

    # Generate returns with changing distributions
    if(crisis_multiplier > 2.0) {
      # Crisis: t-distribution with low degrees of freedom
      returns[t] <- mu + vol_series[t] * rt(1, df = 3)
    } else if(crisis_multiplier > 1.5) {
      # Moderate stress: mixture of normal and t-distribution
      if(runif(1) < 0.7) {
        returns[t] <- mu + vol_series[t] * rnorm(1)
      } else {
        returns[t] <- mu + vol_series[t] * rt(1, df = 5)
      }
    } else {
      # Calm periods: normal distribution
      returns[t] <- mu + vol_series[t] * rnorm(1)
    }

    shock_indicator[t] <- current_shock
  }

  # Add flash crashes (moderate one-day shocks)
  flash_crash_days <- sample(100:2400, 15)
  for(day in flash_crash_days) {
    returns[day] <- returns[day] + sample(c(-1, 1), 1) * 0.03 * (1 + runif(1) * 0.5)  # reduced
    vol_series[day] <- vol_series[day] * 2
    shock_indicator[day] <- 3.0  # reduced from 5.0 to 3.0
  }
  vol_quantiles <- quantile(vol_series, c(0.6, 0.85, 0.95))
  regime <- ifelse(vol_series <= vol_quantiles[1], 1,
            ifelse(vol_series <= vol_quantiles[2], 2,
            ifelse(vol_series <= vol_quantiles[3], 3, 4)))

  data.table(
    date = dates,
    returns = returns,
    true_volatility = vol_series,
    regime = regime,
    shock_indicator = shock_indicator,
    crisis_multiplier = c(1, diff(vol_series) / vol_series[-length(vol_series)] + 1)
  )
}

# Generate data
financial_data <- generate_realistic_financial_data()

#### 2. Volatility Estimation ####

# Realized volatility
calculate_realized_volatility <- function(returns, window = 22) {
  rv <- rollapply(returns^2, window, sum, fill = NA, align = "right") * sqrt(252)
  return(rv)
}

# HAR model
fit_har_model_advanced <- function(volatility_series = train_data$rv_simple, periods = c(1,5, 22, round(365/4),365)) {

  # Data preparation
  dates <- seq.Date(as.Date("2018-01-01"), by = "day", length.out = length(volatility_series)-1)
  volatility <- as.xts(volatility_series[-length(volatility_series)], order.by = dates)

  # Log transformation for stabilization (Clements & Preve 2018)
  log_volatility <- log(pmax(volatility, 1e-6))  # avoid log(0)

  # Create enhanced HAR model with log transformation
  har_model <- highfrequency::HARmodel(
    data = log_volatility,
    periods = periods,
    RVest = c("rCov"),
    type = "HAR",
    h = 1,
    transform = NULL,
    inputType = "RM"
  )

  return(list(
    model = har_model,
    volatility_xts = volatility,
    log_volatility_xts = log_volatility,
    use_log = TRUE
  ))
}

# HAR forecast
forecast_har_with_shock_analysis <- function(har_model_result, train_data, horizon = 300, periods = c(1,5, 22, round(365/4),365)) {

  forecast_vals <- numeric(horizon)
  shock_probabilities <- numeric(horizon)
  regime_probabilities <- matrix(0, nrow = horizon, ncol = 4)

  log_volatility <- har_model_result$log_volatility_xts

  # Analysis of historical shocks
  historical_shocks <- train_data$shock_indicator
  shock_data <- historical_shocks[historical_shocks > 0]

  # Check if enough data for ACF
  if(length(shock_data) > 15) {
    shock_persistence <- acf(shock_data, lag.max = min(10, length(shock_data)/3), plot = FALSE)$acf
  } else {
    # Use simple exponential decay function
    shock_persistence <- c(1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)
  }

  # Shock probability model
  shock_frequency <- mean(historical_shocks > 1.0)  # base shock frequency

  cat("Performing HAR forecast with shock analysis for", horizon, "periods...\n")
  cat("• Historical shock frequency:", round(shock_frequency * 100, 1), "%\n")

  current_shock_state <- 0

  for (i in 1:horizon) {
    # Базовый HAR прогноз
    har_model <- highfrequency::HARmodel(
      data = log_volatility,
      periods = periods,
      RVest = c("rCov"),
      type = "HAR",
      h = 1,
      transform = NULL,
      inputType = "RM"
    )

    log_pred <- predict(har_model)
    base_pred <- exp(as.numeric(log_pred))

    # Анализ вероятности шоков
    # 1. Циклическая компонента (кризисы каждые ~300-400 дней)
    cycle_position <- (nrow(train_data) + i) %% 350
    cycle_shock_prob <- 0.5 * (1 + cos(2 * pi * cycle_position / 350))

    # 2. Персистентность шоков
    if(current_shock_state > 0) {
      persistence_idx <- min(current_shock_state + 1, length(shock_persistence))
      persistence_prob <- ifelse(persistence_idx <= length(shock_persistence),
                                shock_persistence[persistence_idx], 0.01)
      current_shock_state <- max(0, current_shock_state - 1)
    } else {
      persistence_prob <- 0
    }

    # Обрабатываем NA в persistence_prob
    if(is.na(persistence_prob) || !is.finite(persistence_prob)) {
      persistence_prob <- 0
    }

    # 3. Календарные эффекты (конец месяца, квартала)
    day_in_month <- (i %% 22) + 1
    calendar_shock_prob <- ifelse(day_in_month <= 2 || day_in_month >= 20, 0.03, 0.01)

    # 4. Волатильность волатильности
    recent_vol <- mean(as.numeric(tail(log_volatility, 5)), na.rm = TRUE)
    vol_mean <- mean(as.numeric(log_volatility), na.rm = TRUE)
    vol_sd <- sd(as.numeric(log_volatility), na.rm = TRUE)

    if(is.finite(vol_sd) && vol_sd > 0) {
      vol_shock_prob <- abs(exp(recent_vol) - exp(vol_mean)) / exp(vol_sd)
    } else {
      vol_shock_prob <- 0
    }

    # Обрабатываем NA в vol_shock_prob
    if(is.na(vol_shock_prob) || !is.finite(vol_shock_prob)) {
      vol_shock_prob <- 0
    }

    # Общая вероятность шока
    total_shock_prob <- pmin(0.2, shock_frequency + cycle_shock_prob * 0.1 +
                           persistence_prob * 0.08 + calendar_shock_prob + vol_shock_prob)

    # Обрабатываем NA значения
    if(is.na(total_shock_prob) || !is.finite(total_shock_prob)) {
      total_shock_prob <- shock_frequency
    }

    shock_probabilities[i] <- total_shock_prob

    # Моделируем возможный шок (ограничиваем интенсивность)
    if(runif(1) < total_shock_prob) {
      # Генерируем умеренный шок
      shock_intensity <- rgamma(1, shape = 1.5, rate = 2) + 1.9
      shock_direction <- sample(c(0.9, 1.5), 1, prob = c(0.3, 0.7))  # более умеренные изменения

      # Ограничиваем максимальное изменение
      shock_multiplier <- shock_direction * shock_intensity
      base_pred <- base_pred * shock_multiplier
      current_shock_state <- round(shock_intensity)

      cat("  Shock at period", i, ": intensity =", round(shock_intensity, 2), "\n")
    }

    # Анализ режимов волатильности
    vol_level <- base_pred * sqrt(252)  # аннуализируем
    if(vol_level < 0.15) {
      regime_probs <- c(0.7, 0.25, 0.04, 0.01)  # низкая волатильность
    } else if(vol_level < 0.25) {
      regime_probs <- c(0.3, 0.5, 0.15, 0.05)   # умеренная волатильность
    } else if(vol_level < 0.4) {
      regime_probs <- c(0.1, 0.3, 0.45, 0.15)   # высокая волатильность
    } else {
      regime_probs <- c(0.05, 0.15, 0.35, 0.45) # экстремальная волатильность
    }

    regime_probabilities[i, ] <- regime_probs

    # Адаптивная коррекция на основе режима
    regime_adjustment <- ifelse(which.max(regime_probs) == 4, 1.2,    # экстремальный
                         ifelse(which.max(regime_probs) == 3, 1.1,    # высокий
                         ifelse(which.max(regime_probs) == 1, 0.95, 1.0))) # низкий

    base_pred <- base_pred * regime_adjustment


    forecast_vals[i] <- base_pred

    # Добавляем прогноз для следующей итерации
    next_date <- end(log_volatility) + days(1)
    log_volatility <- rbind(log_volatility, xts(log(base_pred), order.by = next_date))

    if (i %% 5 == 0) cat("  Completed period", i, "/", horizon, "\n")
  }

  return(list(
    forecast = forecast_vals,
    shock_probabilities = shock_probabilities,
    regime_probabilities = regime_probabilities,
    shock_analysis = list(
      historical_frequency = shock_frequency,
      persistence = shock_persistence[1:5]
    )
  ))
}

# Enhanced stochastic integration with volatility of volatility
forecast_har_stochastic_advanced <- function(har_forecast_vals =har_forecast , horizon = forecast_horizon, n_simulations = 1000) {

  # Parameters of enhanced Heston process (based on Lu 2022)
  vol_params <- list(
    kappa = 0.5,    # Speed of mean reversion
    theta = mean(har_forecast_vals, na.rm = TRUE),  # Long-term mean
    sigma_v = 0.3,  # Volatility of volatility
    leverage = -0.5 # Leverage effect (asymmetry)
  )

  simulations <- matrix(0, n_simulations, horizon)

  # Constrain base forecasts
  har_forecast_vals_safe <- pmax(1e-4, pmin(har_forecast_vals, mean(har_forecast_vals) * 5))

  for(sim in 1:n_simulations) {
    current_vol <- har_forecast_vals_safe[1]

    for(t in 1:horizon) {
      # HAR forecast for period t
      har_vol_forecast <- har_forecast_vals_safe[t]

      # Moderate stochastic component
      vol_innovation <- rnorm(1, 0, vol_params$sigma_v * sqrt(pmax(1e-4, current_vol)))
      leverage_effect <- vol_params$leverage * sqrt(pmax(1e-4, current_vol)) * rnorm(1)

      stochastic_vol_component <- vol_params$kappa * (vol_params$theta - current_vol) +
                                 vol_innovation + leverage_effect

      # Constrain stochastic component
      stochastic_vol_component <- pmax(-current_vol * 0.5,
                                      pmin(stochastic_vol_component, current_vol * 0.5))

      # Enhanced integration with adaptive weights
      weight_har <- 0.6 + 0.2 * exp(-t/5)
      weight_stoch <- 1 - weight_har

      current_vol <- weight_har * har_vol_forecast + weight_stoch * (current_vol + stochastic_vol_component)

      # Strict constraints
      vol_mean <- mean(har_forecast_vals_safe)
      current_vol <- pmax(vol_mean * 0.1, pmin(current_vol, vol_mean * 3))

      simulations[sim, t] <- current_vol
    }
  }

  # Extended quantiles for better analysis
  list(
    median = apply(simulations, 2, median, na.rm = TRUE),
    q01 = apply(simulations, 2, quantile, 0.01, na.rm = TRUE),
    q05 = apply(simulations, 2, quantile, 0.05, na.rm = TRUE),
    q10 = apply(simulations, 2, quantile, 0.10, na.rm = TRUE),
    q25 = apply(simulations, 2, quantile, 0.25, na.rm = TRUE),
    q75 = apply(simulations, 2, quantile, 0.75, na.rm = TRUE),
    q90 = apply(simulations, 2, quantile, 0.90, na.rm = TRUE),
    q95 = apply(simulations, 2, quantile, 0.95, na.rm = TRUE),
    q99 = apply(simulations, 2, quantile, 0.99, na.rm = TRUE),
    vol_of_vol = apply(simulations, 2, function(x) {
      vol_of_vol_val <- sd(x, na.rm = TRUE)
      if(is.na(vol_of_vol_val) || !is.finite(vol_of_vol_val)) return(0.001)
      return(vol_of_vol_val)  # constrain vol_of_vol
    })
  )
}

#### 3. GARCH for comparison ####
fit_garch_model <- function(returns) {
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(0,0)),
    distribution.model = "norm"
  )

  tryCatch({
    fitted_model <- ugarchfit(spec, returns, solver = "hybrid")
    return(fitted_model)
  }, error = function(e) {
    cat("GARCH did not converge, using simplified model\n")
    return(NULL)
  })
}

#### 4. Apply enhanced methods to data ####
financial_data[, rv_simple := calculate_realized_volatility(returns)]

# Split data
train_size <- 2200
train_data <- financial_data[1:train_size]
test_data <- financial_data[(train_size + 1):nrow(financial_data)]

har_model_result <- fit_har_model_advanced(train_data$rv_simple, periods = c(1,5, 22, round(365/4),365, 365*2))
garch_model <- fit_garch_model(train_data$returns)

# Forecast
forecast_horizon <- nrow(test_data)

har_result_with_shocks <- forecast_har_with_shock_analysis(har_model_result, train_data, forecast_horizon, periods = c(1,5, 22, round(365/4),365, 365*2))
har_forecast <- har_result_with_shocks$forecast
har_stochastic_forecast <- forecast_har_stochastic_advanced(har_forecast, forecast_horizon)

# Create forecast dates
forecast_dates <- test_data$date

# GARCH forecast
if (!is.null(garch_model)) {
  garch_forecast <- ugarchforecast(garch_model, n.ahead = forecast_horizon)
  garch_vol <- as.numeric(sigma(garch_forecast)) * sqrt(252)
} else {
  garch_vol <- rep(mean(train_data$rv_simple, na.rm = TRUE), forecast_horizon)
}

#### 5. Enhanced visualization with volatility of volatility ####
create_advanced_volatility_plot <- function() {

  hist_dates <- train_data$date
  hist_rv <- train_data$rv_simple
  hist_true <- train_data$true_volatility * sqrt(252)

  forecast_dates <- test_data$date
  actual_rv <- test_data$rv_simple
  actual_true <- test_data$true_volatility * sqrt(252)

  p <- plot_ly() %>%

    # Historical data
    add_trace(
      x = hist_dates, y = hist_rv,
      type = 'scatter', mode = 'lines',
      line = list(color = '#2E86C1', width = 2),
      fill = 'tozeroy',
      fillcolor = 'rgba(46, 134, 193, 0.1)',
      name = 'Historical Volatility',
      hovertemplate = '<b>%{x}</b><br>RV: %{y:.3f}<extra></extra>'
    ) %>%

    # True volatility
    add_trace(
      x = hist_dates, y = hist_true,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(128, 128, 128, 0.7)', width = 1, dash = 'dot'),
      name = 'True Volatility',
      hovertemplate = '<b>%{x}</b><br>True: %{y:.3f}<extra></extra>'
    ) %>%

    # HAR forecast
    add_trace(
      x = forecast_dates, y = har_forecast,
      type = 'scatter', mode = 'lines',
      line = list(color = '#E74C3C', width = 3),
      name = 'HAR Forecast (Enhanced)',
      hovertemplate = '<b>%{x}</b><br>HAR: %{y:.3f}<extra></extra>'
    ) %>%

    # Extreme intervals (99%)
    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q99,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(231, 76, 60, 0)', width = 0),
      showlegend = FALSE, hoverinfo = 'skip'
    ) %>%

    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q01,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(231, 76, 60, 0)', width = 0),
      fill = 'tonexty', fillcolor = 'rgba(231, 76, 60, 0.05)',
      name = '99% Extreme Interval',
      hovertemplate = '<b>%{x}</b><br>Extreme CI: %{ymin:.3f} - %{ymax:.3f}<extra></extra>'
    ) %>%

    # 95% interval
    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q95,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(231, 76, 60, 0)', width = 0),
      showlegend = FALSE, hoverinfo = 'skip'
    ) %>%

    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q05,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(231, 76, 60, 0)', width = 0),
      fill = 'tonexty', fillcolor = 'rgba(231, 76, 60, 0.1)',
      name = '90% Confidence Interval',
      hovertemplate = '<b>%{x}</b><br>90% CI: %{ymin:.3f} - %{ymax:.3f}<extra></extra>'
    ) %>%

    # 50% interval
    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q75,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(255, 180, 0, 0)', width = 0),
      showlegend = FALSE, hoverinfo = 'skip'
    ) %>%

    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$q25,
      type = 'scatter', mode = 'lines',
      line = list(color = 'rgba(255, 180, 0, 0)', width = 0),
      fill = 'tonexty', fillcolor = 'rgba(255, 180, 0, 0.25)',
      name = '50% Confidence Interval',
      hovertemplate = '<b>%{x}</b><br>50% CI: %{ymin:.3f} - %{ymax:.3f}<extra></extra>'
    ) %>%

    # Median stochastic forecast
    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$median,
      type = 'scatter', mode = 'lines',
      line = list(color = '#8E44AD', width = 2, dash = 'dash'),
      name = 'HAR Stochastic (Median)',
      hovertemplate = '<b>%{x}</b><br>Median: %{y:.3f}<extra></extra>'
    ) %>%

    # GARCH forecast
    add_trace(
      x = forecast_dates, y = garch_vol,
      type = 'scatter', mode = 'lines',
      line = list(color = '#F39C12', width = 2),
      name = 'GARCH(1,1)',
      hovertemplate = '<b>%{x}</b><br>GARCH: %{y:.3f}<extra></extra>'
    ) %>%

    # Actual volatility
    add_trace(
      x = forecast_dates, y = actual_rv,
      type = 'scatter', mode = 'lines',
      line = list(color = '#27AE60', width = 2),
      fill = 'tozeroy',
      fillcolor = 'rgba(39, 174, 96, 0.1)',
      name = 'Actual Volatility',
      hovertemplate = '<b>%{x}</b><br>Actual: %{y:.3f}<extra></extra>'
    ) %>%

    # Separator line
    add_trace(
      x = c(max(hist_dates), max(hist_dates)),
      y = c(0, max(c(hist_rv, actual_rv), na.rm = TRUE)),
      type = 'scatter', mode = 'lines',
      line = list(color = "gray", width = 2, dash = "dot"),
      showlegend = FALSE,
      name = "Forecast Start",
      hoverinfo = 'skip'
    ) %>%

    # Styling
    layout(
      title = list(
        text = "<b>Enhanced HAR Model: Volatility + Volatility of Volatility</b><br><sub>Log Transformation + Adaptive Integration + Leverage Effect</sub>",
        font = list(size = 18, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "<b>Date</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      yaxis = list(
        title = "<b>Annualized Volatility</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      plot_bgcolor = 'rgba(255, 255, 255, 0)',
      paper_bgcolor = 'rgba(255, 255, 255, 0)',
      legend = list(
        orientation = "v",
        x = 1.02, y = 1,
        font = list(size = 10)
      ),
      hovermode = 'x unified',
      margin = list(l = 60, r = 150, t = 100, b = 60)
    )

  return(p)
}

# Volatility of volatility plot
create_vol_of_vol_plot <- function() {

  p <- plot_ly() %>%
    add_trace(
      x = forecast_dates, y = har_stochastic_forecast$vol_of_vol,
      type = 'scatter', mode = 'lines+markers',
      line = list(color = '#9B59B6', width = 3),
      marker = list(size = 6, color = '#9B59B6'),
      name = 'Volatility of Volatility',
      hovertemplate = '<b>%{x}</b><br>Vol-of-Vol: %{y:.4f}<extra></extra>'
    ) %>%

    layout(
      title = list(
        text = '<b>Volatility of Volatility (Enhanced HAR Model)</b><br><sub>Dynamic Uncertainty Assessment of Forecasts</sub>',
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "<b>Date</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      yaxis = list(
        title = "<b>Volatility of Volatility</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      plot_bgcolor = 'rgba(255, 255, 255, 0)',
      paper_bgcolor = 'rgba(255, 255, 255, 0)',
      margin = list(l = 60, r = 40, t = 100, b = 60)
    )

  return(p)
}

# Shock analysis plot
create_shock_analysis_plot <- function() {

  # Historical shocks
  hist_shock_data <- data.table(
    date = train_data$date,
    shock_intensity = train_data$shock_indicator,
    volatility = train_data$rv_simple
  )

  # Forecasted shocks
  forecast_shock_data <- data.table(
    date = forecast_dates,
    shock_probability = har_result_with_shocks$shock_probabilities * 100,
    volatility_forecast = har_forecast
  )

  p <- plot_ly() %>%

    # Historical shocks (intensity)
    add_trace(
      data = hist_shock_data[shock_intensity > 0],
      x = ~date, y = ~shock_intensity,
      type = 'scatter', mode = 'markers',
      marker = list(
        size = ~pmax(6, shock_intensity * 3),
        color = '#E74C3C',
        opacity = 0.7,
        line = list(color = 'darkred', width = 1)
      ),
      name = 'Historical Shocks',
      hovertemplate = '<b>%{x}</b><br>Intensity: %{y:.2f}<extra></extra>'
    ) %>%

    # Separator line
    add_trace(
      x = c(max(train_data$date), max(train_data$date)),
      y = c(0, 6),
      type = 'scatter', mode = 'lines',
      line = list(color = 'gray', width = 2, dash = 'dot'),
      showlegend = FALSE,
      hoverinfo = 'skip'
    ) %>%

    # Future shock probabilities
    add_trace(
      data = forecast_shock_data,
      x = ~date, y = ~shock_probability,
      type = 'scatter', mode = 'lines+markers',
      line = list(color = '#F39C12', width = 3),
      marker = list(size = 6, color = '#F39C12'),
      name = 'Shock Probability (%)',
      yaxis = 'y2',
      hovertemplate = '<b>%{x}</b><br>Probability: %{y:.1f}%<extra></extra>'
    ) %>%

    layout(
      title = list(
        text = '<b>Shock Analysis of Volatility</b><br><sub>Historical Events and Predicted Probabilities</sub>',
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "<b>Date</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      yaxis = list(
        title = "<b>Historical Shock Intensity</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)',
        side = 'left'
      ),
      yaxis2 = list(
        title = "<b>Future Shock Probability (%)</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        overlaying = 'y',
        side = 'right',
        gridcolor = 'rgba(255, 165, 0, 0.2)'
      ),
      plot_bgcolor = 'rgba(255, 255, 255, 0)',
      paper_bgcolor = 'rgba(255, 255, 255, 0)',
      legend = list(x = 0.02, y = 0.98),
      margin = list(l = 60, r = 60, t = 100, b = 60)
    )

  return(p)
}

# Regime analysis plot
create_regime_analysis_plot <- function() {

  # Historical regimes
  regime_colors <- c('#2ECC71', '#F39C12', '#E67E22', '#E74C3C')
  regime_names <- c('Calm', 'Moderate', 'High', 'Extreme')

  hist_regime_data <- data.table(
    date = train_data$date,
    regime = train_data$regime,
    volatility = train_data$rv_simple * sqrt(252)
  )

  # Predicted regime probabilities
  regime_probs <- har_result_with_shocks$regime_probabilities

  p <- plot_ly() %>%

    # Historical regimes (background)
    add_trace(
      data = hist_regime_data,
      x = ~date, y = ~volatility*100,
      type = 'scatter', mode = 'markers',
      marker = list(
        size = 4,
        color = ~regime,
        colorscale = list(
          c(0, regime_colors[1]), c(0.33, regime_colors[2]),
          c(0.67, regime_colors[3]), c(1, regime_colors[4])
        ),
        opacity = 0.6
      ),
      name = 'Historical Regimes',
      hovertemplate = '<b>%{x}</b><br>Volatility: %{y:.1f}%<br>Regime: %{marker.color}<extra></extra>'
    )

  # Add predicted regime probabilities
  for(regime in 1:4) {
    p <- p %>% add_trace(
      x = forecast_dates,
      y = regime_probs[, regime] * 100,
      type = 'scatter', mode = 'lines',
      line = list(color = regime_colors[regime], width = 2),
      name = paste('Probability:', regime_names[regime]),
      yaxis = 'y2',
      hovertemplate = paste0('<b>%{x}</b><br>', regime_names[regime], ': %{y:.1f}%<extra></extra>')
    )
  }

  p <- p %>%
    # Separator line
    add_trace(
      x = c(max(train_data$date), max(train_data$date)),
      y = c(0, 100),
      type = 'scatter', mode = 'lines',
      line = list(color = 'gray', width = 2, dash = 'dot'),
      showlegend = FALSE,
      hoverinfo = 'skip'
    ) %>%

    layout(
      title = list(
        text = '<b>Regime Analysis</b><br><sub>Historical Regimes and Predicted Probabilities</sub>',
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "<b>Date</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)'
      ),
      yaxis = list(
        title = "<b>Historical Volatility (%)</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        gridcolor = 'rgba(128, 128, 128, 0.2)',
        side = 'left'
      ),
      yaxis2 = list(
        title = "<b>Regime Probabilities (%)</b>",
        titlefont = list(size = 14),
        tickfont = list(size = 12),
        overlaying = 'y',
        side = 'right',
        range = c(0, 100)
      ),
      plot_bgcolor = 'rgba(255, 255, 255, 0)',
      paper_bgcolor = 'rgba(255, 255, 255, 0)',
      legend = list(x = 0.02, y = 0.98),
      margin = list(l = 60, r = 80, t = 100, b = 60)
    )

  return(p)
}

#### 6. Forecast Accuracy Analysis ####
calculate_forecast_accuracy <- function(actual, forecast, model_name) {
  actual_clean <- actual[!is.na(actual)]
  forecast_clean <- forecast[1:length(actual_clean)]

  mse <- mean((actual_clean - forecast_clean)^2)
  mae <- mean(abs(actual_clean - forecast_clean))
  mape <- mean(abs((actual_clean - forecast_clean) / actual_clean)) * 100

  data.table(
    model = model_name,
    mse = mse,
    mae = mae,
    mape = mape
  )
}

# Compare models
actual_test_rv <- test_data$rv_simple[!is.na(test_data$rv_simple)]
forecast_length <- length(actual_test_rv)

accuracy_results <- rbind(
  calculate_forecast_accuracy(actual_test_rv, har_forecast[1:forecast_length], "Enhanced HAR"),
  calculate_forecast_accuracy(actual_test_rv, har_stochastic_forecast$median[1:forecast_length], "Enhanced Stochastic HAR"),
  calculate_forecast_accuracy(actual_test_rv, garch_vol[1:forecast_length], "GARCH(1,1)")
)

#### 7. Display Results ####

cat("FORECAST ACCURACY:\n")
print(accuracy_results)

# Main volatility plot
main_plot <- create_advanced_volatility_plot()
print(main_plot)

# Shock analysis plot
shock_plot <- create_shock_analysis_plot()
print(shock_plot)

# Regime analysis plot
regime_plot <- create_regime_analysis_plot()
print(regime_plot)

# Volatility of volatility plot
vol_of_vol_plot <- create_vol_of_vol_plot()
print(vol_of_vol_plot)
