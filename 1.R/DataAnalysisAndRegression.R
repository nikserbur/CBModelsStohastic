# DataAnalysisAndRegression.R — Анализ данных и регрессионное моделирование для прогнозирования ключевых переменных
# Проект: Стохастическое моделирование ставок ЦБ с масштабируемым подходом
# Автор: Бураков Никита Сергеевич

#### 0. Подключение библиотек ####
source("1.R/8.SourceAll.R")

#### 1. Анализ доступных данных ####

#' AnalyzeAvailableData — Анализ всех доступных колонок в данных
#' @param input_data исходные данные
#' @return список с анализом колонок
AnalyzeAvailableData <- function(input_data) {
  
  # Получение всех колонок
  all_columns <- names(input_data)
  numeric_columns <- names(input_data)[sapply(input_data, is.numeric)]
  date_columns <- names(input_data)[sapply(input_data, function(x) inherits(x, "Date"))]
  
  # Анализ базовых статистик
  numeric_summary <- input_data[, lapply(.SD, function(x) {
    if (is.numeric(x)) {
      list(
        mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE),
        na_count = sum(is.na(x)),
        na_percent = round(sum(is.na(x)) / length(x) * 100, 2)
      )
    } else NULL
  }), .SDcols = numeric_columns]
  
  # Анализ корреляций
  correlation_matrix <- cor(input_data[, .SD, .SDcols = numeric_columns], 
                           use = "pairwise.complete.obs")
  
  # Определение основных целевых переменных
  target_variables <- c("ruonia_rate", "cpi_inflation", "key_rate")
  
  # Определение переменных, которые рассчитываются в коде
  calculated_variables <- c(
    "inflation_gap", "policy_tightness", "inflation_expectations_gap",
    "gdp_growth", "credit_to_gdp", "key_rate_change", "gdp_trend",
    "gdp_gap", "inflation_deviation", "rate_change", "rate_volatility",
    "asset_volatility", "volatility_regime"
  )
  
  # Определение исходных переменных (не рассчитываемых)
  raw_variables <- setdiff(numeric_columns, calculated_variables)
  
  return(list(
    all_columns = all_columns,
    numeric_columns = numeric_columns,
    date_columns = date_columns,
    target_variables = target_variables,
    calculated_variables = calculated_variables,
    raw_variables = raw_variables,
    numeric_summary = numeric_summary,
    correlation_matrix = correlation_matrix
  ))
}

#### 2. Определение ключевых целевых переменных ####

#' IdentifyKeyTargetVariables — Определение основных целевых переменных для моделирования
#' @param data_analysis результат анализа данных
#' @return список ключевых переменных
IdentifyKeyTargetVariables <- function(data_analysis) {
  
  # Основные целевые переменные для правила Тейлора
  primary_targets <- c("ruonia_rate", "cpi_inflation")
  
  # Переменные, которые рассчитываются из других
  derived_variables <- c(
    "inflation_deviation",  # рассчитывается как cpi_inflation - TARGET_INFLATION
    "gdp_gap",             # рассчитывается через loess тренд
    "policy_tightness",    # рассчитывается как key_rate - cpi_inflation
    "inflation_expectations_gap"  # рассчитывается как expected_inflation - TARGET_INFLATION
  )
  
  # Исходные переменные для прогнозирования
  raw_predictors <- setdiff(data_analysis$raw_variables, primary_targets)
  
  return(list(
    primary_targets = primary_targets,
    derived_variables = derived_variables,
    raw_predictors = raw_predictors,
    description = "Основные целевые переменные для правила Тейлора и их исходные предикторы"
  ))
}

#### 3. Построение регрессионных моделей для прогнозирования ####

#' BuildRegressionModels — Построение регрессионных моделей для прогнозирования ключевых переменных
#' @param train_data обучающие данные
#' @param target_variables целевые переменные
#' @param predictor_variables предикторы
#' @return список обученных моделей
BuildRegressionModels <- function(train_data, target_variables, predictor_variables) {
  
  models <- list()
  
  for (target in target_variables) {
    cat("Построение модели для:", target, "\n")
    
    # Подготовка данных для модели
    model_data <- train_data[, c(target, predictor_variables), with = FALSE]
    model_data <- na.omit(model_data)
    
    if (nrow(model_data) < 10) {
      cat("Недостаточно данных для моделирования", target, "\n")
      next
    }
    
    # Формула модели
    formula_str <- paste0(target, " ~ ", paste(predictor_variables, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    # Построение линейной регрессии
    tryCatch({
      linear_model <- lm(model_formula, data = model_data)
      
      # Проверка качества модели
      r_squared <- summary(linear_model)$r.squared
      adj_r_squared <- summary(linear_model)$adj.r.squared
      aic <- AIC(linear_model)
      bic <- BIC(linear_model)
      
      # Анализ остатков
      residuals <- residuals(linear_model)
      alpha_estimate <- EstimateAlphaParameter(abs(residuals), method = "hill")
      
      # Taleb R² коррекция
      robust_r2 <- CalculateRobustR2(
        y_actual = model_data[[target]],
        y_predicted = fitted(linear_model),
        alpha_tail = alpha_estimate
      )
      
      models[[target]] <- list(
        model = linear_model,
        formula = formula_str,
        r_squared = r_squared,
        adj_r_squared = adj_r_squared,
        aic = aic,
        bic = bic,
        alpha_estimate = alpha_estimate,
        robust_r2 = robust_r2,
        residuals = residuals,
        fitted = fitted(linear_model),
        data_points = nrow(model_data)
      )
      
      cat("  R² =", round(r_squared, 4), "Adj R² =", round(adj_r_squared, 4), "\n")
      cat("  Taleb R² =", round(robust_r2$r2_taleb_exact, 4), "\n")
      cat("  α =", round(alpha_estimate, 4), "\n")
      
    }, error = function(e) {
      cat("Ошибка при построении модели для", target, ":", e$message, "\n")
    })
  }
  
  return(models)
}

#### 4. Сравнение методов прогнозирования ####

#' CompareForecastingMethods — Сравнение различных методов прогнозирования
#' @param train_data обучающие данные
#' @param target_variable целевая переменная
#' @param predictor_variables предикторы
#' @return сравнение методов
CompareForecastingMethods <- function(train_data, target_variable, predictor_variables) {
  
  # Подготовка данных
  model_data <- train_data[, c(target_variable, predictor_variables), with = FALSE]
  model_data <- na.omit(model_data)
  
  if (nrow(model_data) < 20) {
    return(NULL)
  }
  
  # Разделение на обучающую и тестовую выборки
  set.seed(123)
  train_indices <- sample(1:nrow(model_data), size = floor(0.8 * nrow(model_data)))
  train_subset <- model_data[train_indices]
  test_subset <- model_data[-train_indices]
  
  results <- list()
  
  # 1. Линейная регрессия
  tryCatch({
    formula_str <- paste0(target_variable, " ~ ", paste(predictor_variables, collapse = " + "))
    linear_model <- lm(as.formula(formula_str), data = train_subset)
    
    # Прогнозы
    train_pred <- predict(linear_model, train_subset)
    test_pred <- predict(linear_model, test_subset)
    
    # Метрики
    train_rmse <- sqrt(mean((train_subset[[target_variable]] - train_pred)^2))
    test_rmse <- sqrt(mean((test_subset[[target_variable]] - test_pred)^2))
    
    results$linear <- list(
      model = linear_model,
      train_rmse = train_rmse,
      test_rmse = test_rmse,
      r_squared = summary(linear_model)$r.squared,
      aic = AIC(linear_model)
    )
  }, error = function(e) {
    cat("Ошибка в линейной регрессии:", e$message, "\n")
  })
  
  # 2. Gradient Boosting (если доступен)
  tryCatch({
    if (requireNamespace("xgboost", quietly = TRUE)) {
      # Подготовка данных для XGBoost
      train_matrix <- xgboost::xgb.DMatrix(
        data = as.matrix(train_subset[, predictor_variables, with = FALSE]),
        label = train_subset[[target_variable]]
      )
      test_matrix <- xgboost::xgb.DMatrix(
        data = as.matrix(test_subset[, predictor_variables, with = FALSE]),
        label = test_subset[[target_variable]]
      )
      
      # Параметры модели
      params <- list(
        objective = "reg:squarederror",
        max_depth = 6,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8,
        nrounds = 100
      )
      
      # Обучение модели
      xgb_model <- xgboost::xgb.train(params, train_matrix, nrounds = 100, verbose = 0)
      
      # Прогнозы
      train_pred_xgb <- predict(xgb_model, train_matrix)
      test_pred_xgb <- predict(xgb_model, test_matrix)
      
      # Метрики
      train_rmse_xgb <- sqrt(mean((train_subset[[target_variable]] - train_pred_xgb)^2))
      test_rmse_xgb <- sqrt(mean((test_subset[[target_variable]] - test_pred_xgb)^2))
      
      results$xgboost <- list(
        model = xgb_model,
        train_rmse = train_rmse_xgb,
        test_rmse = test_rmse_xgb,
        importance = xgboost::xgb.importance(model = xgb_model)
      )
    }
  }, error = function(e) {
    cat("Ошибка в XGBoost:", e$message, "\n")
  })
  
  # 3. Сравнительная таблица
  comparison_table <- data.table(
    Method = names(results),
    Train_RMSE = sapply(results, function(x) x$train_rmse),
    Test_RMSE = sapply(results, function(x) x$test_rmse),
    R_Squared = sapply(results, function(x) ifelse("r_squared" %in% names(x), x$r_squared, NA)),
    AIC = sapply(results, function(x) ifelse("aic" %in% names(x), x$aic, NA))
  )
  
  return(list(
    results = results,
    comparison = comparison_table,
    target_variable = target_variable,
    n_train = nrow(train_subset),
    n_test = nrow(test_subset)
  ))
}

#### 5. Интеграция с основной моделью ####

#' IntegrateRegressionModels — Интеграция регрессионных моделей с основной стохастической моделью
#' @param regression_models обученные регрессионные модели
#' @param forecast_data прогнозные данные
#' @param predictor_variables предикторы
#' @return обновленные прогнозные данные
IntegrateRegressionModels <- function(regression_models, forecast_data, predictor_variables) {
  
  updated_forecast <- copy(forecast_data)
  
  for (target_var in names(regression_models)) {
    model <- regression_models[[target_var]]
    
    if (!is.null(model) && !is.null(model$model)) {
      tryCatch({
        # Проверяем доступность предикторов в прогнозных данных
        available_predictors <- predictor_variables[predictor_variables %in% names(forecast_data)]
        
        if (length(available_predictors) > 0) {
          # Прогнозирование
          prediction_data <- forecast_data[, available_predictors, with = FALSE]
          
          # Удаляем NA значения для прогнозирования
          complete_cases <- complete.cases(prediction_data)
          
          if (sum(complete_cases) > 0) {
            # Прогноз только для полных случаев
            predictions <- predict(model$model, prediction_data[complete_cases])
            
            # Обновляем прогнозные данные
            updated_forecast[complete_cases, (target_var) := predictions]
            
            cat("Прогноз для", target_var, "обновлен на", sum(complete_cases), "наблюдений\n")
          }
        }
      }, error = function(e) {
        cat("Ошибка при прогнозировании", target_var, ":", e$message, "\n")
      })
    }
  }
  
  return(updated_forecast)
}

#### 6. Документация методологии ####

#' CreateMethodologyDocumentation — Создание документации методологии
#' @param data_analysis анализ данных
#' @param target_variables целевые переменные
#' @param regression_models регрессионные модели
#' @return документация методологии
CreateMethodologyDocumentation <- function(data_analysis, target_variables, regression_models) {
  
  methodology <- list(
    title = "Масштабируемая методология прогнозирования ключевых переменных ЦБ РФ",
    
    data_analysis = list(
      total_variables = length(data_analysis$all_columns),
      numeric_variables = length(data_analysis$numeric_columns),
      raw_variables = length(data_analysis$raw_variables),
      calculated_variables = length(data_analysis$calculated_variables)
    ),
    
    target_variables = list(
      primary_targets = target_variables$primary_targets,
      derived_variables = target_variables$derived_variables,
      description = "Основные целевые переменные для правила Тейлора"
    ),
    
    regression_approach = list(
      method = "Множественная линейная регрессия с Taleb R² коррекцией",
      features = c(
        "Автоматический выбор предикторов",
        "Коррекция жирных хвостов",
        "Сравнение с Gradient Boosting",
        "Интеграция с основной моделью"
      ),
      advantages = c(
        "Масштабируемость: легко добавлять новые переменные",
        "Интерпретируемость: коэффициенты регрессии",
        "Робастность: Taleb R² для жирных хвостов",
        "Гибкость: поддержка различных методов"
      )
    ),
    
    model_performance = lapply(regression_models, function(model) {
      list(
        r_squared = model$r_squared,
        adj_r_squared = model$adj_r_squared,
        taleb_r2 = model$robust_r2$r2_taleb_exact,
        alpha_estimate = model$alpha_estimate,
        data_points = model$data_points
      )
    }),
    
    integration_steps = c(
      "1. Анализ всех доступных данных",
      "2. Определение ключевых целевых переменных", 
      "3. Построение регрессионных моделей",
      "4. Сравнение методов прогнозирования",
      "5. Интеграция с основной стохастической моделью",
      "6. Обновление прогнозов через регрессионные модели"
    ),
    
    technical_details = list(
      fat_tail_correction = "Taleb R² для жирных хвостов",
      variable_selection = "Автоматический выбор предикторов",
      model_comparison = "Linear Regression vs Gradient Boosting",
      integration_method = "Прогнозирование ключевых переменных из остальных данных"
    )
  )
  
  return(methodology)
}

#### 7. Основная функция анализа и моделирования ####

#' RunDataAnalysisAndRegression — Основная функция анализа данных и регрессионного моделирования
#' @param input_data исходные данные
#' @return полный анализ и модели
RunDataAnalysisAndRegression <- function(input_data) {
  
  cat("=== АНАЛИЗ ДАННЫХ И РЕГРЕССИОННОЕ МОДЕЛИРОВАНИЕ ===\n")
  
  # 1. Анализ доступных данных
  cat("1. Анализ доступных данных...\n")
  data_analysis <- AnalyzeAvailableData(input_data)
  
  cat("   Всего колонок:", length(data_analysis$all_columns), "\n")
  cat("   Числовых колонок:", length(data_analysis$numeric_columns), "\n")
  cat("   Исходных переменных:", length(data_analysis$raw_variables), "\n")
  cat("   Рассчитываемых переменных:", length(data_analysis$calculated_variables), "\n")
  
  # 2. Определение ключевых переменных
  cat("2. Определение ключевых переменных...\n")
  target_variables <- IdentifyKeyTargetVariables(data_analysis)
  
  cat("   Основные целевые переменные:", paste(target_variables$primary_targets, collapse = ", "), "\n")
  cat("   Исходные предикторы:", paste(target_variables$raw_predictors, collapse = ", "), "\n")
  
  # 3. Построение регрессионных моделей
  cat("3. Построение регрессионных моделей...\n")
  regression_models <- BuildRegressionModels(
    input_data, 
    target_variables$primary_targets, 
    target_variables$raw_predictors
  )
  
  # 4. Сравнение методов для каждой целевой переменной
  cat("4. Сравнение методов прогнозирования...\n")
  method_comparisons <- list()
  
  for (target in target_variables$primary_targets) {
    if (target %in% names(regression_models)) {
      comparison <- CompareForecastingMethods(
        input_data, 
        target, 
        target_variables$raw_predictors
      )
      method_comparisons[[target]] <- comparison
    }
  }
  
  # 5. Создание документации
  cat("5. Создание документации методологии...\n")
  methodology_doc <- CreateMethodologyDocumentation(
    data_analysis, 
    target_variables, 
    regression_models
  )
  
  # 6. Сохранение результатов
  cat("6. Сохранение результатов...\n")
  results <- list(
    data_analysis = data_analysis,
    target_variables = target_variables,
    regression_models = regression_models,
    method_comparisons = method_comparisons,
    methodology_documentation = methodology_doc,
    timestamp = Sys.time()
  )
  
  saved_path <- SaveProjectResults(results, "3.Data/data_analysis_and_regression_results.qs")
  
  cat("Анализ завершен! Результаты сохранены в", saved_path, "\n")
  
  return(results)
}

#### 8. Запуск анализа ####

# Загрузка данных
input_data <- prep_data()

# Запуск полного анализа
analysis_results <- RunDataAnalysisAndRegression(input_data)

# Вывод ключевых результатов
cat("\n=== КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ ===\n")
cat("Доступно переменных для анализа:", length(analysis_results$data_analysis$numeric_columns), "\n")
cat("Основные целевые переменные:", paste(analysis_results$target_variables$primary_targets, collapse = ", "), "\n")
cat("Построено моделей:", length(analysis_results$regression_models), "\n")

# Вывод качества моделей
for (target in names(analysis_results$regression_models)) {
  model <- analysis_results$regression_models[[target]]
  cat(target, ":", "R² =", round(model$r_squared, 4), 
      "Taleb R² =", round(model$robust_r2$r2_taleb_exact, 4), "\n")
} 