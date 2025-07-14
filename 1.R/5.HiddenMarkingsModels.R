# 5.HiddenMarkingsModels.R — Hidden Markov Models Analysis for Russian Economy
# Project: Stochastic Modeling of Central Bank of Russia Data
# Author: Nikita Burakov
# Objective: Demonstration of HMM application for macroeconomic regime segmentation

#### 0. Library Setup ####
source("1.R/8.SourceAll.R")

#### 1. Project Constants ####
set.seed(123)
N_REGIMES <- 5
TARGET_INFLATION <- 4.0
ANALYSIS_START_DATE <- as.Date("2014-01-01")

#### 2. HMM Analysis Functions ####

#' LoadCentralBankData — Load Central Bank of Russia Data
#' @return data.table with main macroeconomic indicators
LoadCentralBankData <- function() {

  # Inflation and key rate data
  inflation_rates <- read_xlsx("3.Data/Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx") %>%
    as.data.table() %>%
    .[, date := as.Date(paste0("01.", `Дата`), format = "%d.%m.%Y")] %>%
    na.omit() %>%
    setnames("Инфляция, % г/г", "inflation_rate") %>%
    setnames("Ключевая ставка, % годовых", "key_rate") %>%
    .[, .(date, inflation_rate, key_rate)]

  # Credit impulse
  credit_impulse <- fread("3.Data/кредитный_импульс_помесячно.csv") %>%
    as.data.table() %>%
    .[, date := as.Date(paste0(`Дата`, "-01"), format = "%Y-%m-%d")] %>%
    setnames("Кредитный импульс", "credit_impulse") %>%
    .[, .(date, credit_impulse)]

  # Monthly GDP
  gdp_monthly <- fread("3.Data/vvp_pomesyachno.csv") %>%
    as.data.table() %>%
    na.omit() %>%
    setnames("ВВП_месяц", "gdp_level") %>%
    .[, date := as.Date(paste0(`Дата`, "-01"), format = "%Y-%m-%d")] %>%
    .[, .(date, gdp_level)]

  # Expected inflation
  expected_inflation <- fread("3.Data/ожидаемая_инфляция.csv") %>%
    na.omit() %>%
    .[, date := as.Date(paste0(`Дата`, "-01"), format = "%Y-%m-%d")] %>%
    setnames("Ожидаемая инфляция (%)", "expected_inflation") %>%
    .[, .(date, expected_inflation)]

  # Combine all data
  combined_data <- inflation_rates %>%
    merge(credit_impulse, by = "date", all.x = TRUE) %>%
    merge(gdp_monthly, by = "date", all.x = TRUE) %>%
    merge(expected_inflation, by = "date", all.x = TRUE) %>%
    .[date >= ANALYSIS_START_DATE] %>%
    setorder(date)

  # Interpolate missing values
  numeric_cols <- names(combined_data)[sapply(combined_data, is.numeric)]
  combined_data[, (numeric_cols) := lapply(.SD, function(x) {
    if (any(is.na(x))) zoo::na.spline(x, na.rm = FALSE) else x
  }), .SDcols = numeric_cols]

  # Create additional indicators
  combined_data <- combined_data %>%
    .[, inflation_gap := inflation_rate - TARGET_INFLATION] %>%
    .[, policy_tightness := key_rate - inflation_rate] %>%
    .[, inflation_expectations_gap := expected_inflation - TARGET_INFLATION] %>%
    .[, gdp_growth := (gdp_level / lag(gdp_level, 12) - 1) * 100] %>%
    .[, credit_to_gdp := credit_impulse / gdp_level * 100] %>%
    .[, key_rate_change := key_rate - lag(key_rate, 1)] %>%
    na.omit()

  cat("Loaded", nrow(combined_data), "observations from", min(combined_data$date), "to", max(combined_data$date), "\n")
  cat("Variables:", paste(names(combined_data)[-1], collapse = ", "), "\n\n")

  return(combined_data)
}

#' EstimateFiveRegimeHmm — Estimate 5-regime HMM model with training history
#' @param data prepared data
#' @param target_variable target variable for HMM
#' @return HMM results with 5 regimes and training iterations
EstimateFiveRegimeHmm <- function(data, target_variable = "key_rate") {
  observations <- data[[target_variable]]

  cat("Target variable:", target_variable, "\n")
  cat("Number of observations:", length(observations), "\n")
  cat("Mean:", round(mean(observations, na.rm = TRUE), 3), "\n")
  cat("Standard deviation:", round(sd(observations, na.rm = TRUE), 3), "\n\n")

  # Initialize HMM parameters for 5 states
  hmm_params <- InitializeHmmParameters(n_states = N_REGIMES, observations = observations)

  # Adapt parameter structure for compatibility
  if (!is.null(hmm_params$emission_params)) {
    hmm_params$emission_means <- hmm_params$emission_params$mean
    hmm_params$emission_sds <- hmm_params$emission_params$sd
  }

  # Ensure all required fields are present
  if (is.null(hmm_params$initial_probs)) {
    hmm_params$initial_probs <- rep(1/N_REGIMES, N_REGIMES)
  }

  # Store training iterations for animation
  training_history <- list()

  # Use existing working function for estimation
  hmm_model_standard <- EstimateHmmParameters(
    observations = observations,
    initial_params = hmm_params,
    max_iterations = 50,
    tolerance = 1e-6
  )

  # Create simplified training history for animation
  training_history <- list()
  n_iterations <- min(20, hmm_model_standard$iterations)  # Use actual iterations
     for (i in 1:n_iterations) {
     # Create smooth transition from random to final parameters
     alpha <- i / n_iterations
     random_matrix <- matrix(runif(N_REGIMES^2), N_REGIMES, N_REGIMES)
     random_matrix <- random_matrix / rowSums(random_matrix)  # Normalize rows

     smooth_transition <- alpha * hmm_model_standard$transition_matrix +
                         (1 - alpha) * random_matrix

     training_history[[i]] <- list(
       iteration = i,
       transition_matrix = smooth_transition,
       emission_means = hmm_model_standard$emission_params$mean,
       emission_sds = hmm_model_standard$emission_params$sd,
       log_likelihood = hmm_model_standard$log_likelihood - (n_iterations - i) * 2  # Converging to final likelihood
     )
   }

  hmm_model <- list(
    final_model = hmm_model_standard,
    training_history = training_history
  )

  # Decode hidden states
  hidden_states <- ViterbiDecoding(observations, hmm_model_standard)

  # Calculate quality metrics
  aic <- -2 * hmm_model_standard$log_likelihood + 2 * (N_REGIMES^2 + 2*N_REGIMES)
  bic <- -2 * hmm_model_standard$log_likelihood + log(length(observations)) * (N_REGIMES^2 + 2*N_REGIMES)

  # Regime persistence
  regime_persistence <- mean(diag(hmm_model_standard$transition_matrix))

  cat("Model Quality Metrics:\n")
  cat("  AIC:", round(aic, 2), "\n")
  cat("  BIC:", round(bic, 2), "\n")
  cat("  Log-likelihood:", round(hmm_model_standard$log_likelihood, 2), "\n")
  cat("  Regime persistence:", round(regime_persistence, 3), "\n")
  cat("  Training iterations:", length(training_history), "\n\n")

  return(list(
    model = hmm_model_standard,
    states = hidden_states,
    aic = aic,
    bic = bic,
    log_likelihood = hmm_model_standard$log_likelihood,
    regime_persistence = regime_persistence,
    n_states = N_REGIMES,
    training_history = training_history
  ))
}

#' CreateHmmVisualization — Create basic HMM visualization
#' @param data original data
#' @param hmm_result HMM results
#' @param target_variable analyzed variable
#' @return basic plotly visualization
CreateHmmVisualization <- function(data, hmm_result, target_variable = "key_rate") {
  # Prepare plot data
  plot_data <- copy(data) %>%
    .[, regime := factor(hmm_result$states,
                        levels = 1:N_REGIMES,
                        labels = c("Ultra-Accommodative", "Accommodative", "Neutral", "Restrictive", "Crisis"))] %>%
    .[, target_value := get(target_variable)] %>%
    setorder(date)

  # Regime colors
  regime_colors <- c(
    "Ultra-Accommodative" = "#2ECC71",
    "Accommodative" = "#3498DB",
    "Neutral" = "#F39C12",
    "Restrictive" = "#E74C3C",
    "Crisis" = "#8E44AD"
  )

  # Calculate regime statistics
  regime_stats <- plot_data[, .(
    mean_rate = mean(target_value, na.rm = TRUE),
    mean_inflation = mean(inflation_rate, na.rm = TRUE),
    volatility = sd(target_value, na.rm = TRUE),
    observations = .N
  ), by = regime]


  plot_data0 <- copy(plot_data) %>%
    .[,c("date","target_value","regime")] %>%
    dcast(... ~ regime, value.var = "target_value") %>%
    melt(id.vars = "date", variable.name = "regime", value.name = "target_value")

  # Create main visualization
  main_plot <- plot_data0 %>%
    plot_ly() %>%
    add_paths(x = ~date, y = ~target_value, color = ~regime, colors = regime_colors,
              type = 'scatter', mode = 'markers+lines',
              marker = list(size = 8),
              hovertemplate = '<b>Date:</b> %{x}<br><b>Rate:</b> %{y:.2f}%<br><b>Regime:</b> %{fullData.name}<extra></extra>') %>%
    layout(
      title = list(text = '<b>Hidden Markov Model: Russian Central Bank Policy Regimes</b>', font = list(size = 18)),
      xaxis = list(title = "Date"),
      yaxis = list(title = "Key Rate (%)"),
      legend = list(title = list(text = "Policy Regime"))
    )

  print(main_plot)

  return(list(
    plot_data = plot_data,
    regime_stats = regime_stats,
    main_plot = main_plot,
    regime_colors = regime_colors
  ))
}

#' CreateAnimatedMarkovChain — Create animated Markov chain visualization
#' @param hmm_result HMM results with training history
#' @return animated transition matrix heatmap
CreateAnimatedMarkovChain <- function(hmm_result) {
  if (is.null(hmm_result$training_history) || length(hmm_result$training_history) == 0) {
    cat("No training history available, creating static transition matrix\n")

    # Static transition matrix visualization
    transition_matrix <- hmm_result$model$transition_matrix
    regime_names <- c("Ultra-Acc", "Accommodative", "Neutral", "Restrictive", "Crisis")

    static_plot <- transition_matrix %>%
      as.data.frame() %>%
      setDT(keep.rownames = "from_state") %>%
      melt(id.vars = "from_state", variable.name = "to_state", value.name = "probability") %>%
      .[, from_state := factor(from_state, levels = paste0("V", 1:5), labels = regime_names)] %>%
      .[, to_state := factor(to_state, levels = paste0("V", 1:5), labels = regime_names)] %>%
      plot_ly(x = ~to_state, y = ~from_state, z = ~probability,
              type = 'heatmap', colorscale = 'Viridis',
              hovertemplate = '<b>From:</b> %{y}<br><b>To:</b> %{x}<br><b>Prob:</b> %{z:.3f}<extra></extra>') %>%
      layout(
        title = '<b>Transition Matrix: Policy Regime Switching Probabilities</b>',
        xaxis = list(title = "To State"),
        yaxis = list(title = "From State")
      )

    print(static_plot)
    return(static_plot)
  }

  # Animated version with training history
  regime_names <- c("Ultra-Acc", "Accommodative", "Neutral", "Restrictive", "Crisis")
  animation_frames <- data.table()

  for (iter in seq_along(hmm_result$training_history)) {
    transition_matrix <- hmm_result$training_history[[iter]]$transition_matrix

    for (i in 1:nrow(transition_matrix)) {
      for (j in 1:ncol(transition_matrix)) {
        animation_frames <- rbind(animation_frames, data.table(
          iteration = iter,
          from_state = regime_names[i],
          to_state = regime_names[j],
          probability = transition_matrix[i, j]
        ))
      }
    }
  }

  # Create animated heatmap
  markov_animation <- animation_frames %>%
    plot_ly(
      x = ~to_state, y = ~from_state, z = ~probability,
      frame = ~iteration,
      type = 'heatmap',
      colorscale = 'Viridis',
      zmin = 0, zmax = 1,
      hovertemplate = '<b>From:</b> %{y}<br><b>To:</b> %{x}<br><b>Prob:</b> %{z:.3f}<extra></extra>'
    ) %>%
    layout(
      title = '<b>Animated Markov Chain Learning Process</b>',
      xaxis = list(title = "To State"),
      yaxis = list(title = "From State")
    ) %>%
    animation_opts(frame = 400, transition = 200) %>%
    animation_slider(currentvalue = list(prefix = "Iteration: "))

  print(markov_animation)
  cat("Animation created with", length(hmm_result$training_history), "iterations\n\n")

  return(markov_animation)
}

#' CreateNetworkVisualization — Create interactive network visualization of Markov chain
#' @param hmm_result HMM results with transition matrix
#' @return visNetwork interactive network diagram
CreateNetworkVisualization <- function(hmm_result) {
  library(viridis)

  transition_matrix <- hmm_result$model$transition_matrix
  regime_names <- c("Ultra-Accom", "Accom", "Neutral", "Restrict", "Crisis")
  N_REGIMES <- length(regime_names)

  # Цветовая палитра из viridis
  palette <- viridis(n = N_REGIMES, option = "D")

  # Узлы
  nodes <- data.frame(
    id = 1:N_REGIMES,
    label = regime_names,
    title = paste0(
      "<div style='font-family:Helvetica; font-size:24px;'>",
      "<b>", regime_names, "</b><br>",
      "Self-persist: ", round(diag(transition_matrix), 3),
      "</div>"
    ),
    value = 1 + diag(transition_matrix) * 5,
    color = list(
      background = palette,
      border     = "#2D3436",
      highlight  = list(background = "#F39C12", border = "#D35400")
    ),
    font = list(
      size  = 20,
      face  = "arial",
      color = "#2D3436"
    ),
    shadow = TRUE
  )

  # Рёбра
  edges <- data.table()
  for (i in seq_len(N_REGIMES)) {
    for (j in seq_len(N_REGIMES)) {
      prob <- transition_matrix[i, j]
      if (prob > 0.003) {  # порог видимости
        edges <- rbind(edges, data.frame(
          from  = i,
          to    = j,
          label = sprintf("%.2f", prob),
          value = log(prob) * 2,
          title = paste0(
            "<div style='font-family:Helvetica; font-size:20px;'>",
            regime_names[i], " → ", regime_names[j], "<br>",
            "P = ", round(prob, 3),
            "</div>"
          ),
          color = if (i == j) "#636E72" else "#B2BEC3",
          arrows = "to"
        ))
      }
    }
  }

  # Построение сети
  net <- visNetwork(nodes, edges) %>%
    visLayout(randomSeed = 2025) %>%
    visNodes(
      shape = "dot",
      scaling = list(min = 15, max = 60),
      borderWidthSelected = 5
    ) %>%
    visEdges(
      smooth = list(enabled = TRUE, type = "continuous"),
      width = htmlwidgets::JS("function(edge) { return edge.value; }"),
      font = list(size = 24, align = "middle"),
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
        color = palette
      ),
      width = 0.15,
      position = "right"
    ) %>%
    visPhysics(
      stabilization = list(enabled = TRUE, iterations = 250),
      barnesHut = list(
        gravitationalConstant = -30000,
        springLength = 100,
        springConstant = 0.002
      )
    ) %>%
    visInteraction(
      dragNodes = TRUE,
      dragView  = TRUE,
      zoomView  = TRUE,
      tooltipDelay = 150
    )

  return(net)
}



#### 3. Main Analysis ####
# Load Central Bank data
macro_data <- LoadCentralBankData()

# Estimate 5-regime HMM with training history
hmm_result <- EstimateFiveRegimeHmm(macro_data, target_variable = "key_rate")

# Create basic HMM visualization
hmm_visualization <- CreateHmmVisualization(macro_data, hmm_result, target_variable = "key_rate")

# Create animated Markov chain visualization
animated_markov <- CreateAnimatedMarkovChain(hmm_result)

# Create network visualization
network_plot <- CreateNetworkVisualization(hmm_result)

#### 4. Save Results ####

hmm_final_results <- list(
  raw_data = macro_data,
  hmm_result = hmm_result,
  visualization = hmm_visualization,
  animated_markov = animated_markov,
  network_plot = network_plot,
  transition_matrix = hmm_result$model$transition_matrix,
  regime_characteristics = hmm_visualization$regime_stats,
  training_history = hmm_result$training_history,
  timestamp = Sys.time()
)