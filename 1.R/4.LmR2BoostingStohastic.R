# Author: Nikita Burakov
# Date: 2025-07-11
# Version: 2.1
# Fat-tailed residuals in linear regression and R² correction

#### 0. Load Libraries ####

library(data.table)
library(plotly)
library(extraDistr)
library(moments)
library(magrittr)
library(stats)  # For t-distribution and F-distribution functions

#### 1. Function Definitions ####

# Generate different types of residuals for simulation
GenerateResiduals <- function(n, residual_type = "normal") {
  switch(residual_type,
    "normal" = rnorm(n, 0, 0.5),
    "fat_tail" = {
      base_error <- rt(n, df = ALPHA_TAIL) * 0.45
      extreme_mask <- runif(n) < 0.04  # 4% outliers
      base_error[extreme_mask] <- base_error[extreme_mask] * 3.5
      base_error
    },
    "mixture" = {
      is_outlier <- rbinom(n, 1, 0.07)  # 7% outliers
      normal_part <- rnorm(n, 0, 0.5)
      outlier_part <- rt(n, df = 3) * 1.5  # Moderately heavy t-distribution
      (1 - is_outlier) * normal_part + is_outlier * outlier_part
    }
  )
}

# Fit standard OLS model and extract key metrics
FitOlsModel <- function(dt_subset) {
  model <- lm(y ~ x1 + x2, data = dt_subset)

  list(
    coefficients = coef(model),
    r_squared = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared,
    residuals = residuals(model),
    fitted = fitted(model),
    sigma = summary(model)$sigma
  )
}

# Estimation of α parameter for power law (Hill method)
EstimateAlpha <- function(data) {
  data_positive <- data[data > 0]
  if(length(data_positive) < 10) return(NA)

  # Hill method for α estimation
  n <- length(data_positive)
  k <- min(100, n %/% 4)  # use top 25% or 100 observations

  data_sorted <- sort(data_positive, decreasing = TRUE)
  alpha_hill <- 1 / mean(log(data_sorted[1:k]) - log(data_sorted[k+1]))

  return(alpha_hill)
}

# R² correction for fat tails
CalculateRobustR2 <- function(y_actual, y_fitted, alpha_tail = 3.5) {
  residuals <- y_actual - y_fitted
  residuals_abs <- abs(residuals)
  n <- length(residuals)

  # Standard R² = 1 - SS_res/SS_tot
  ss_res <- sum(residuals^2)
  ss_tot <- sum((y_actual - mean(y_actual))^2)
  r2_standard <- 1 - ss_res/ss_tot
  if(alpha_tail > 2) {  # Fat tails with finite second moment

    # Method 1: Taleb's exact R² formula (6.14) for t-distributed errors
    # When ε ~ t(α), then E[ε²] = α/(α-2) for α > 2
    # Formula: R² = a²(α-2)/(a²(α-2) + α)

    # Estimate regression coefficients variance from data
    model_temp <- lm(y_actual ~ y_fitted)
    coef_variance <- var(coef(model_temp)[2])  # variance of slope coefficient

    # Apply Taleb's exact formula (6.14): R² = a²(α-2)/(a²(α-2) + α)
    if(alpha_tail > 2.1) {
      # Calculate true signal to noise ratio from regression
      signal_variance <- var(y_fitted)
      observed_noise_variance <- var(residuals)

      # Theoretical correction: E[ε²] = α/(α-2) for t-distribution
      theoretical_noise_multiplier <- alpha_tail / (alpha_tail - 2)

      # Apply Taleb's exact formula
      # R² = signal_var / (signal_var + theoretical_noise_var)
      theoretical_noise_variance <- observed_noise_variance * theoretical_noise_multiplier
      r2_taleb_exact <- signal_variance / (signal_variance + theoretical_noise_variance)
      r2_taleb_exact <- max(0, min(1, r2_taleb_exact))
    } else {
      r2_taleb_exact <- 0  # As α approaches 2, E[R²] approaches 0
    }

    # Method 2: MLE correction accounting for E[ε²] = α/(α-2)
    if(alpha_tail > 2.1) {
      # Correct the sum of squares using theoretical expectation
      expected_ss_multiplier <- alpha_tail / (alpha_tail - 2)
      corrected_ss_res <- ss_res * expected_ss_multiplier
      r2_mle_corrected <- max(0, 1 - corrected_ss_res / ss_tot)
    } else {
      r2_mle_corrected <- 0
    }

  } else {
    # For α ≤ 2, infinite variance case
    r2_taleb_exact <- 0
    r2_mle_corrected <- 0
  }

  list(
    r2_standard = r2_standard,
    r2_taleb_exact = r2_taleb_exact,
    r2_mle_corrected = r2_mle_corrected,
    alpha_estimated = EstimateAlpha(residuals_abs),
    degrees_freedom = alpha_tail,
    expected_error_variance = if(alpha_tail > 2) alpha_tail / (alpha_tail - 2) else Inf
  )
}

#### 2. Parameter Setup and Base Data Generation ####

set.seed(123)  # Standard seed for reproducible demo
N_OBS <- 10000  # Large sample for stable results
ALPHA_TAIL <- 2.8  # Optimal fat tail parameter for demonstration
# Create true relationship Y = β₀ + β₁X₁ + β₂X₂ + ε
true_beta0 <- 2.0
true_beta1 <- 1.5
true_beta2 <- -0.8

dt_base <- data.table(
  x1 = rnorm(N_OBS, mean = 0, sd = 1),
  x2 = rnorm(N_OBS, mean = 0, sd = 1.5)
) %>%
  .[, y_true := true_beta0 + true_beta1 * x1 + true_beta2 * x2]

#### 3. Generate Models with Different Residual Types ####
# Create comprehensive dataset using melt approach
dt_models <- copy(dt_base) %>%
  .[,`Normal Residuals` := GenerateResiduals(N_OBS, "normal")] %>%
  .[,`Fat Tails` := GenerateResiduals(N_OBS, "fat_tail")] %>%
  .[,`Mixture` := GenerateResiduals(N_OBS, "mixture")] %>%
  melt(id.vars = c("x1", "x2", "y_true"), variable.name = "model_type", value.name = "residuals") %>%
  .[, y := y_true + residuals]


#### 4. Fit Standard OLS Models ####

# Fit models for each type
models_results <- dt_models[, {
  model_fit <- FitOlsModel(.SD)
  list(
    beta0 = model_fit$coefficients[1],
    beta1 = model_fit$coefficients[2],
    beta2 = model_fit$coefficients[3],
    r2_standard = model_fit$r_squared,
    adj_r2 = model_fit$adj_r_squared,
    sigma_residuals = model_fit$sigma,
    residuals = list(model_fit$residuals),
    fitted = list(model_fit$fitted)
  )
}, by = model_type]

print(models_results[, .(model_type, beta0, beta1, beta2, r2_standard, sigma_residuals)])

#### 5. Apply R² Corrections ####

# Apply robust estimates with correct α for each model type
robust_results <- dt_models[, {
  subset_data <- .SD
  model_fit <- FitOlsModel(subset_data)

  # Use appropriate α based on model type
  model_type_char <- as.character(unique(.BY$model_type))
  alpha_for_model <- switch(model_type_char,
    "Normal Residuals" = 30,  # Large α for normal distribution
    "Fat Tails" = ALPHA_TAIL,  # Use defined fat tail parameter
    "Mixture" = max(4, ALPHA_TAIL * 0.8)  # Intermediate value for mixture
  )

  robust_metrics <- CalculateRobustR2(subset_data$y, model_fit$fitted, alpha_for_model)

  list(
    r2_standard = robust_metrics$r2_standard,
    r2_taleb_exact = robust_metrics$r2_taleb_exact,
    r2_mle_corrected = robust_metrics$r2_mle_corrected,
    alpha_estimated = robust_metrics$alpha_estimated,
    degrees_freedom = robust_metrics$degrees_freedom,
    expected_error_variance = robust_metrics$expected_error_variance,
    alpha_used = alpha_for_model
  )
}, by = model_type]

final_comparison <- merge(models_results[, .(model_type, sigma_residuals)],
                         robust_results, by = "model_type")

print(final_comparison)

#### 6. Visualization of Results ####

# Prepare data for plots
residuals_data <- dt_models[, {
  model_fit <- FitOlsModel(.SD)
  data.table(
    residuals = model_fit$residuals,
    fitted = model_fit$fitted,
    abs_residuals = abs(model_fit$residuals)
  )
}, by = model_type]

# Plot 1: Q-Q plot of residuals
plot_qq_residuals <- residuals_data %>%
  split(by = "model_type") %>%
  lapply(function(dt) {
    sample_q <- quantile(dt$residuals, ppoints(100))
    theory_q <- qnorm(ppoints(100), 0, sd(dt$residuals))

    plot_ly(x = ~theory_q, y = ~sample_q, type = "scatter", mode = "markers",
            name = unique(dt$model_type), opacity = 0.7) %>%
      add_lines(x = theory_q, y = theory_q, name = "Perfect Line",
               line = list(dash = "dash", color = "red"))
  }) %>%
  subplot(nrows = 1, margin = 0.05, titleX = TRUE, titleY = TRUE) %>%
  layout(title = "Q-Q Plots of Residuals vs Normal Distribution")

plot_residuals_dist <- plot_ly() %>%
  add_histogram(data = residuals_data[model_type == "Normal Residuals"], x = ~residuals,
               name = "Normal Residuals", opacity = 0.6, nbinsx = 60,
               histnorm = "probability density",
               marker = list(color = "#1f77b4", line = list(color = "#0d47a1", width = 1))) %>%
  add_histogram(data = residuals_data[model_type == "Fat Tails"], x = ~residuals,
               name = "Fat Tails", opacity = 0.6, nbinsx = 60,
               histnorm = "probability density",
               marker = list(color = "#d62728", line = list(color = "#b71c1c", width = 1))) %>%
  add_histogram(data = residuals_data[model_type == "Mixture"], x = ~residuals,
               name = "Mixture Model", opacity = 0.6, nbinsx = 60,
               histnorm = "probability density",
               marker = list(color = "#2ca02c", line = list(color = "#1b5e20", width = 1))) %>%
    layout(
    title = "Residual Distributions by Model Type",
    xaxis = list(title = "Residuals", range = c(-10, 10)),
    yaxis = list(title = "Probability Density"),
    barmode = "overlay",
    showlegend = TRUE
  )

# Plot 3: Tail analysis
tail_probs <- residuals_data[, {
  abs_res <- abs(residuals)
  sigma_est <- sd(residuals)

  data.table(
    sigma_2 = mean(abs_res > 2 * sigma_est),
    sigma_3 = mean(abs_res > 3 * sigma_est),
    sigma_4 = mean(abs_res > 4 * sigma_est),
    sigma_5 = mean(abs_res > 5 * sigma_est)
  )
}, by = model_type] %>%
  melt(id.vars = "model_type", variable.name = "threshold", value.name = "probability")

plot_tail_probs <- plot_ly(tail_probs, x = ~threshold, y = ~probability,
                          color = ~model_type, type = 'bar', barmode = 'group') %>%
  layout(
    title = "Tail Probabilities of Residuals",
    yaxis = list(title = "Probability"),
    xaxis = list(title = "Threshold (in σ)")
  )

# Plot 4: Dramatic R² comparison showing overestimation
r2_comparison <- final_comparison %>%
  melt(id.vars = "model_type",
       measure.vars = c("r2_standard", "r2_taleb_exact", "r2_mle_corrected"),
       variable.name = "r2_type", value.name = "r2_value")

# Add better labels for methods
r2_comparison[, r2_type_label := factor(r2_type,
  levels = c("r2_standard", "r2_taleb_exact", "r2_mle_corrected"),
  labels = c("Standard R²", "Taleb Exact R²", "MLE Corrected R²"))]

# Remove NA values for plotting
r2_comparison <- r2_comparison[!is.na(r2_value)]

plot_r2_comparison <- plot_ly(r2_comparison, x = ~model_type, y = ~r2_value,
                             color = ~r2_type_label, type = 'bar', barmode = 'group',
                             opacity = 0.8) %>%
  layout(
    title = "Conservative R² Correction Methods vs Standard R²",
    yaxis = list(title = "R² Value"),
    xaxis = list(title = "Model Type"),
    legend = list(title = "Method")
  )

# Plot 5: Method comparison
effectiveness_data <- final_comparison %>%
  .[!is.na(r2_taleb_exact) & !is.na(r2_mle_corrected)] %>%
  .[, .(
    model_type,
    difference_exact = r2_standard - r2_taleb_exact,
    difference_mle = r2_standard - r2_mle_corrected
  )] %>%
  melt(id.vars = "model_type",
       measure.vars = c("difference_exact", "difference_mle"),
       variable.name = "method", value.name = "difference")

plot_effectiveness <- plot_ly(effectiveness_data, x = ~model_type, y = ~difference,
                             color = ~method, type = 'bar', barmode = 'group',
                             opacity = 0.7) %>%
  layout(
    title = "Difference between Standard and Taleb Methods",
    yaxis = list(title = "R² Difference"),
    xaxis = list(title = "Model Type"),
    legend = list(title = "Method")
  )

#### 7. Final Dashboard ####
subplot(
  subplot(plot_residuals_dist, plot_tail_probs, nrows = 1, margin = 0.03),
  subplot(plot_r2_comparison, plot_effectiveness, nrows = 1, margin = 0.03),
  subplot(plot_qq_residuals, nrows = 1, margin = 0.03),
  nrows = 3, margin = 0.05, titleX = TRUE, titleY = TRUE
) %>%
  layout(title = "Fat-Tailed Residuals Analysis: Taleb's Exact R² Correction Methods")

#### 8. Final Results ####
# Calculate differences for interpretation
differences_summary <- final_comparison[, .(
  model_type,
  alpha_used,
  standard_r2 = round(r2_standard, 4),
  taleb_exact_r2 = round(r2_taleb_exact, 4),
  mle_corrected_r2 = round(r2_mle_corrected, 4),
  bias_exact = round(r2_standard - r2_taleb_exact, 4),
  bias_mle = round(r2_standard - r2_mle_corrected, 4),
  theoretical_e_eps2 = round(expected_error_variance, 2)
)]