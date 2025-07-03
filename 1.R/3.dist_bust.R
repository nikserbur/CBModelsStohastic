# Автор:  Никита Бураков (адаптация для R)
# Дата: 2025-06-25
# Версия: 1.4
# =============================================================================

#### 0. Подключаем всё нужное ####

# Если какие-то пакеты не стоят — поставьте вручную через install.packages()
library(data.table)
library(plotly)
library(extraDistr)
library(moments)

#### 1. Настройка параметров и генерация данных ####

set.seed(42)  # фиксируем сид — для повторяемости

N_SAMPLES <- 100000
BASE_VOL <- 0.3
VOL_VOL <- 0.6
CRASH_PROB <- 0.01
CRASH_MULT <- 5

cat("Генерим данные...\n")

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

cat("Готово! Пример данных:\n")
print(head(dt, 5))

#### 2. Готовим данные для графиков ####

log_dt <- melt(dt[, .(
  normal = log(abs(normal)),
  stochastic = log(abs(stochastic)),
  mixture = log(abs(mixture)))
], measure.vars = c("normal", "stochastic", "mixture"),
.variable.name = "distribution", value.name = "log_return")

density_dt <- melt(dt, measure.vars = c("normal", "stochastic", "mixture"),
                   variable.name = "distribution", value.name = "return")

#### 3. Рисуем графики плотности — в логе и без ####

# A. Логарифмы абсолютных значений
plot_log_density <- plot_ly(log_dt, x = ~log_return, color = ~variable, type = 'histogram', histnorm = "probability density", opacity = 0.6) %>%
  layout(
    title = "Логарифмы абсолютных доходностей",
    xaxis = list(title = "log(|return|)"),
    yaxis = list(title = "Плотность"),
    barmode = "overlay"
  )

# B. Центр распределений
plot_density <- plot_ly(density_dt, x = ~return, color = ~distribution, type = 'histogram', histnorm = "probability density", opacity = 0.6) %>%
  layout(
    title = "Центральные части распределений",
    xaxis = list(title = "Доходность", range = c(-2, 2)),
    yaxis = list(title = "Плотность"),
    barmode = "overlay"
  )

#### 4. Делаем QQ-графики — вручную через quantile ####

create_qqplot_plotly <- function(data, title) {
  theory_quantiles <- qnorm(ppoints(100), mean = mean(data), sd = sd(data))
  sample_quantiles <- quantile(data, probs = ppoints(100))

  res <- plot_ly() %>%
    add_trace(
    x = theory_quantiles, y = sample_quantiles,
    type = "scatter", mode = "markers",
    marker = list(color = 'steelblue', size = 5, opacity = 0.6)) %>%
    add_lines(
      x = theory_quantiles, y = theory_quantiles, name = "Идеал"
    ) %>%
    layout(
      title = title,
      xaxis = list(title = "Теоретические квантили"),
      yaxis = list(title = "Выборочные квантили")
    )
}

qq1 <- create_qqplot_plotly(dt$normal, "QQ-график: Нормальное")
qq2 <- create_qqplot_plotly(dt$stochastic, "QQ-график: Стохастическая волатильность")
qq3 <- create_qqplot_plotly(dt$mixture, "QQ-график: Смесь распределений")

#### 5. Показываем всё в одной дашке ####

cat("Показываем графики...\n")
subplot(
  plot_log_density,
  subplot(qq1, qq2, qq3, nrows = 1, margin = 0.03, titleX = TRUE, titleY = TRUE),
  plot_density,
  nrows = 3,
  margin = 0.07,
  titleX = TRUE,
  titleY = TRUE
) %>% layout(title = "Сравнение распределений с разными хвостами (Талеб-style)")

#### 6. Хвостовой анализ ####

tail_analysis <- function(data, name) {
  abs_data <- abs(data)
  data.table(
    Распределение = name,
    `P(|X|>2σ)` = mean(abs_data > 2 * BASE_VOL),
    `P(|X|>3σ)` = mean(abs_data > 3 * BASE_VOL),
    `P(|X|>4σ)` = mean(abs_data > 4 * BASE_VOL),
    `P(|X|>5σ)` = mean(abs_data > 5 * BASE_VOL),
    Эксцесс = kurtosis(data),
    Среднее = mean(data),
    Дисперсия = var(data)
  )
}

results <- rbind(
  tail_analysis(dt$normal, "Нормальное"),
  tail_analysis(dt$stochastic, "Стохастическая"),
  tail_analysis(dt$mixture, "Смесь")
)

cat("Результаты хвостов:\n")
print(results)

theoretical <- data.table(
  Событие = c(">2σ", ">3σ", ">4σ", ">5σ"),
  Вероятность = c(
    2 * (1 - pnorm(2)),
    2 * (1 - pnorm(3)),
    2 * (1 - pnorm(4)),
    2 * (1 - pnorm(5)))
)

tail_dt <- melt(results, id.vars = "Распределение",
                measure.vars = c("P(|X|>2σ)", "P(|X|>3σ)", "P(|X|>4σ)", "P(|X|>5σ)"),
                variable.name = "Событие", value.name = "Вероятность")

# Переводим значения по X в текст
tail_dt$Событие <- gsub("P\\(\\|X\\|>", ">", tail_dt$Событие)
tail_dt$Событие <- gsub("σ\\)", "σ", tail_dt$Событие)

# Отрисовка в plotly
plot_ly(tail_dt, x = ~`Событие`, y = ~`Вероятность`,
        color = ~`Распределение`, type = 'bar',
        barmode = 'group') %>%
  layout(
    title = "Хвостовые вероятности: эксперимент",
    yaxis = list(type = "log", title = "Вероятность (лог масштаб)"),
    xaxis = list(title = "Событие"),
    legend = list(orientation = 'h')
  )
