# Demo_SVD_Forecasting.R
# SVD forecasting of macroeconomic indicators
# Author: Burakov N.S.

source("1.R/1.functions/all_fun.R")

#### Data Generation ####
set.seed(123)
n_months <- 120
start_date <- as.Date("2014-01-01")
dates <- seq(start_date, by = "month", length.out = n_months)
t <- 1:n_months

generate_macro_data <- function() {
  # Inflation with strong annual seasonality
  inflation_base <- 4 +
    1.8 * sin(2 * pi * t / 12) +                    # enhanced annual seasonality
    0.3 * sin(2 * pi * t / 6) -                    # reduced semi-annual
    0.4 * sin(2 * pi * t / 60) +
    0.3 * sin(2 * pi * t / 36) +
    0.1 * sin(2 * pi * t / 3) +                     # weak quarterly
    cumsum(rnorm(n_months, 0, 0.3)) +               # random walk
    rnorm(n_months, 0, 0.4)                         # additional noise

  # GDP with annual seasonality
  gdp_growth <- 8 +
    1.4 * sin(2 * pi * t / 12 + pi/3) +            # strong annual seasonality
    0.2 * sin(2 * pi * t / 4) -                     # weak quarterly
    0.8 * sin(2 * pi * t / 48) +
    0.3 * cos(2 * pi * t / 84) +
    cumsum(rnorm(n_months, 0, 0.4)) +               # random walk
    rnorm(n_months, 0, 0.3)                         # noise

  # Key rate (Taylor rule)
  key_rate <- 10 +
    0.8 * (inflation_base - 4) +
    -0.4 * (gdp_growth - 2.5) +
    0.8 * sin(2 * pi * t / 12) +                    # annual CB seasonality
    0.6 * sin(2 * pi * t / 24) +
    cumsum(rnorm(n_months, 0, 0.2)) +               # random walk
    rnorm(n_months, 0, 0.25)                        # decision noise

  # Exchange rate
  exchange_rate <- 65 +
    3.5 * sin(2 * pi * t / 18) +                    # oil cycle
    2 * cos(2 * pi * t / 42) +
    2.0 * sin(2 * pi * t / 12) +                    # annual export seasonality
    0.6 * (inflation_base - 4) +
    cumsum(rnorm(n_months, 0, 1.0)) +               # random walk
    rnorm(n_months, 0, 1.2)                         # volatility

  # Credit growth
  credit_growth <- 12 +
    -0.5 * (key_rate - 6) +
    0.4 * (gdp_growth - 2.5) +
    1.8 * sin(2 * pi * t / 12 + pi/6) -            # strong annual seasonality
    0.2 * sin(2 * pi * t / 4) +                     # weak quarterly
    cumsum(rnorm(n_months, 0, 0.5)) +               # random walk
    rnorm(n_months, 0, 0.6)                         # noise

  data.table(
    date = dates,
    inflation = inflation_base,
    gdp_growth = gdp_growth,
    key_rate = key_rate,
    exchange_rate = exchange_rate,
    credit_growth = credit_growth
  )
}

macro_data <- generate_macro_data()

#### SVD Forecasting ####
train_size <- 96
train_data <- macro_data[1:train_size]
test_data <- macro_data[(train_size + 1):nrow(macro_data)]

svd_result <- PredictTimeSeriesWithSVD(
  input_dataset = train_data,
  prediction_steps = nrow(test_data),
  temporal_var = "date",
  seasonal_period = 12,
  enable_seasonality = TRUE,
  enable_cyclical = TRUE
)

forecast_dt <- svd_result$forecast %>%
  .[,date := test_data$date]
historical_dt <- svd_result$historical %>%
  .[,date := train_data$date]
var_names <- c("inflation", "gdp_growth", "key_rate", "exchange_rate", "credit_growth")

#### Visualization ####
create_individual_plots <- function(historical = train_data, forecast = forecast_dt, test_actual = test_data, i = 1) {
  plot_titles <- c(
    "Inflation, % y/y",
    "GDP Growth, % y/y",
    "Key Rate, %",
    "USD/RUB",
    "Credit Growth, % y/y"
  )

  var <- var_names[i]
  title <- plot_titles[i]

  p <- plot_ly() %>%
    add_trace(
      data = historical,
      x = ~date, y = ~get(var),
      type = 'scatter', mode = 'lines',
      line = list(color = '#2E86C1', width = 3),
      fill = 'tozeroy',
      fillcolor = 'rgba(46, 134, 193, 0.15)',
      name = 'Historical',
      hovertemplate = '<b>%{x}</b><br>%{y:.2f}<extra></extra>'
    ) %>%
    add_trace(
      data = forecast,
      x = ~date, y = ~get(var),
      type = 'scatter', mode = 'lines',
      line = list(color = '#E74C3C', width = 3, dash = 'dash'),
      fill = 'tozeroy',
      fillcolor = 'rgba(231, 76, 60, 0.15)',
      name = 'Forecast',
      hovertemplate = '<b>%{x}</b><br>%{y:.2f}<extra></extra>'
    ) %>%
    add_trace(
      data = test_actual,
      x = ~date, y = ~get(var),
      type = 'scatter', mode = 'lines',
      line = list(color = '#27AE60', width = 3, dash = 'dot'),
      name = 'Actual',
      hovertemplate = '<b>%{x}</b><br>%{y:.2f}<extra></extra>'
    ) %>%
    layout(
      title = list(text = paste0("<b>", title, "</b>"), font = list(size = 16)),
      xaxis = list(title = "", gridcolor = 'rgba(128, 128, 128, 0.2)'),
      yaxis = list(title = "", gridcolor = 'rgba(128, 128, 128, 0.2)'),
      plot_bgcolor = 'white',
      paper_bgcolor = 'white',
      showlegend = FALSE,
      margin = list(l = 50, r = 40, t = 50, b = 60))


  return(p)
}


subplot(create_individual_plots(train_data, forecast_dt, test_data,1),create_individual_plots(train_data, forecast_dt, test_data,2),
        create_individual_plots(train_data, forecast_dt, test_data,3),create_individual_plots(train_data, forecast_dt, test_data,4),
        create_individual_plots(train_data, forecast_dt, test_data,5), nrows = 3)

