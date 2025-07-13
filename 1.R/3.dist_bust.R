# Author: Nikita Burakov (R adaptation)
# Date: 2025-06-25
# Version: 1.5

#### 0. Load Required Packages ####

# Подключение всех библиотек проекта
source("1.R/9.libraries.r")

#### 1. Parameter Setup and Data Generation ####

set.seed(42)  # fix seed for reproducibility

N_SAMPLES <- 100000
BASE_VOL <- 0.3
VOL_VOL <- 0.6
CRASH_PROB <- 0.01
CRASH_MULT <- 5

cat("Generating data...\n")

dt <- data.table(
  normal = rnorm(N_SAMPLES, mean = 0, sd = BASE_VOL),

  stochastic = replicate(N_SAMPLES, {
    stochastic_vol <- BASE_VOL * exp(rnorm(1, mean = 0, sd = VOL_VOL))
    rnorm(1, mean = 0, sd = stochastic_vol)
  }),

  mixture = replicate(N_SAMPLES, {
    if (runif(1) < CRASH_PROB) {
      rnorm(1, mean = 0, sd = BASE_VOL * CRASH_MULT)
    } else {
      stochastic_vol <- BASE_VOL * exp(rnorm(1, mean = 0, sd = VOL_VOL))
      rnorm(1, mean = 0, sd = stochastic_vol)
    }
  })
)

cat("Done! Sample data:\n")
print(head(dt, 5))

#### 2. Prepare Data for Plots ####

log_dt <- melt(dt[, .(
  normal = log(abs(normal)),
  stochastic = log(abs(stochastic)),
  mixture = log(abs(mixture)))
], measure.vars = c("normal", "stochastic", "mixture"),
.variable.name = "distribution", value.name = "log_return")

density_dt <- melt(dt, measure.vars = c("normal", "stochastic", "mixture"),
                   variable.name = "distribution", value.name = "return")

#### 3. Plot Density Graphs - Logarithmic and Linear ####

# A. Logarithms of absolute values
plot_log_density <- plot_ly(log_dt, x = ~log_return, color = ~variable, type = 'histogram', histnorm = "probability density", opacity = 0.6) %>%
  layout(
    title = "Logarithms of Absolute Returns",
    xaxis = list(title = "log(|return|)"),
    yaxis = list(title = "Density"),
    barmode = "overlay"
  )

# B. Central parts of distributions
plot_density <- plot_ly(density_dt, x = ~return, color = ~distribution, type = 'histogram', histnorm = "probability density", opacity = 0.6) %>%
  layout(
    title = "Central Parts of Distributions",
    xaxis = list(title = "Return", range = c(-2, 2)),
    yaxis = list(title = "Density"),
    barmode = "overlay"
  )

#### 4. Create QQ-Plots Manually via Quantile ####

create_qqplot_plotly <- function(data, title) {
  theory_quantiles <- qnorm(ppoints(100), mean = mean(data), sd = sd(data))
  sample_quantiles <- quantile(data, probs = ppoints(100))

  res <- plot_ly() %>%
    add_trace(
    x = theory_quantiles, y = sample_quantiles,
    type = "scatter", mode = "markers",
    marker = list(color = 'steelblue', size = 5, opacity = 0.6)) %>%
    add_lines(
      x = theory_quantiles, y = theory_quantiles, name = "Perfect Line"
    ) %>%
    layout(
      title = title,
      xaxis = list(title = "Theoretical Quantiles"),
      yaxis = list(title = "Sample Quantiles")
    )
}

qq1 <- create_qqplot_plotly(dt$normal, "QQ-Plot: Normal")
qq2 <- create_qqplot_plotly(dt$stochastic, "QQ-Plot: Stochastic Volatility")
qq3 <- create_qqplot_plotly(dt$mixture, "QQ-Plot: Mixture Distribution")

#### 5. Display All in Single Dashboard ####

cat("Displaying plots...\n")
subplot(
  plot_log_density,
  subplot(qq1, qq2, qq3, nrows = 1, margin = 0.03, titleX = TRUE, titleY = TRUE),
  plot_density,
  nrows = 3,
  margin = 0.07,
  titleX = TRUE,
  titleY = TRUE
) %>% layout(title = "Comparison of Distributions with Different Tails (Taleb-style)")

#### 6. Tail Analysis ####

tail_analysis <- function(data, name) {
  abs_data <- abs(data)
  data.table(
    Distribution = name,
    `P(|X|>2σ)` = mean(abs_data > 2 * BASE_VOL),
    `P(|X|>3σ)` = mean(abs_data > 3 * BASE_VOL),
    `P(|X|>4σ)` = mean(abs_data > 4 * BASE_VOL),
    `P(|X|>5σ)` = mean(abs_data > 5 * BASE_VOL),
    Kurtosis = kurtosis(data),
    Mean = mean(data),
    Variance = var(data)
  )
}

results <- rbind(
  tail_analysis(dt$normal, "Normal"),
  tail_analysis(dt$stochastic, "Stochastic"),
  tail_analysis(dt$mixture, "Mixture")
)

cat("Tail analysis results:\n")
print(results)

theoretical <- data.table(
  Event = c(">2σ", ">3σ", ">4σ", ">5σ"),
  Probability = c(
    2 * (1 - pnorm(2)),
    2 * (1 - pnorm(3)),
    2 * (1 - pnorm(4)),
    2 * (1 - pnorm(5)))
)

tail_dt <- melt(results, id.vars = "Distribution",
                measure.vars = c("P(|X|>2σ)", "P(|X|>3σ)", "P(|X|>4σ)", "P(|X|>5σ)"),
                variable.name = "Event", value.name = "Probability")

# Convert X values to text
tail_dt$Event <- gsub("P\\(\\|X\\|>", ">", tail_dt$Event)
tail_dt$Event <- gsub("σ\\)", "σ", tail_dt$Event)

# Plot with plotly
plot_ly(tail_dt, x = ~`Event`, y = ~`Probability`,
        color = ~`Distribution`, type = 'bar',
        barmode = 'group') %>%
  layout(
    title = "Tail Probabilities: Experimental Results",
    yaxis = list(type = "log", title = "Probability (log scale)"),
    xaxis = list(title = "Event"),
    legend = list(orientation = 'h')
  )
