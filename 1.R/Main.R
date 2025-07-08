# Main.R — Основной файл для расчётов и моделирования
# Проект: Стохастическое моделирование ставок
# Разработчик: Бураков Никита Сергеевич

# 0. Подключение всех библиотек ####
source("1.R/8.SourceAll.R")

# 1. Данные (загрузка, очистка, создание лагов/гепов) ####
# исходные данные тут https://www.cbr.ru/dkp/statistics/

# Данные и инфляции
inflation <- read_xlsx("3.Data/Инфляция и ключевая ставка Банка России_F01_12_2013_T18_06_2025.xlsx") %>%
  as.data.table() %>%
  .[,`Дата` := as.Date(paste0("01.", `Дата`), format = "%d.%m.%Y")] %>%
  na.omit() %>%
  setnames("Инфляция, % г/г","cpi_inflation")

# Кредитный импульс
credit_impl <- fread("3.Data/кредитный_импульс_помесячно.csv") %>%
  as.data.table() %>%
  .[,`Дата` := as.Date(paste0(`Дата`,"-01"), format = "%Y-%m-%d")]

# Межбанковская реальная ставка
st_ruonia <- read_xlsx("3.Data/RC_F11_01_2010_T04_07_2025.xlsx") %>%
  as.data.table() %>%
  .[,c("DT","ruo")] %>%
  .[, ym := format(as.Date(DT), "%Y-%m")] %>%
  .[, .(`ruonia_rate` = mean(ruo, na.rm = TRUE)), by = .(ym)] %>%
  .[,`Дата` := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")]

#ВВП
gdp <- fread("3.Data/vvp_pomesyachno.csv") %>%
  as.data.table() %>%
  na.omit() %>%
  setnames("ВВП_месяц","gdp_level")

#Бюджетный импулс
budget_impl <- read_xlsx("3.Data/budget.xlsx") %>%
  as.data.table() %>%
  setnames("Date","Дата") %>%
  .[,`Дата` := as.Date(paste0(`Дата` ,"-01"), format = "%Y-%m-%d")] %>%
  setnames("Value","Бюджетный импулс")

#Скор санкций
sunctions <- fread("3.Data/sunctions.csv") %>%
  setnames(c("Дата","Скор санкций")) %>%
  .[, ym := format(as.Date(`Дата`), "%Y-%m")] %>%
  .[, .(`Скор_санкций` = mean(`Скор санкций`, na.rm = TRUE)), by = .(ym)] %>%
  .[,`Дата` := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")]

respectation_inf <- fread("3.Data/ожидаемая_инфляция.csv") %>%
  na.omit() %>%
  .[,`Дата` := as.Date(paste0(`Дата`,"-01"), format = "%Y-%m-%d")]

#Валатильность активов
wind <- 30
all_input_data <- read_xlsx("3.Data/assets.xlsx") %>%
  as.data.table() %>%
  na.omit() %>%
  .[order(`Дата`)] %>%
  .[, all := MOEX + RTS] %>%
  .[, vol := frollapply(all, wind, function(all){
    mean <- mean(all)
    res <- sum(abs(mean - all))/length(all)
    })*sqrt(wind)] %>%
  .[,c("Дата","vol")] %>%
  na.omit() %>%
  .[, ym := format(as.Date(`Дата`), "%Y-%m")] %>%
  .[, .(vol = mean(vol, na.rm = TRUE)), by = .(ym)] %>%
  .[,`Дата` := as.Date(paste0(ym,"-01"), format = "%Y-%m-%d")] %>%
  .[,-c("ym")] %>%
  merge(inflation, by = "Дата", all.x = TRUE) %>%
  merge(credit_impl, by = "Дата", all.x = TRUE) %>%
  merge(st_ruonia, by = "Дата", all.x = TRUE) %>%
  merge(gdp, by = "Дата", all.x = TRUE) %>%
  merge(budget_impl, by = "Дата", all.x = TRUE) %>%
  # merge(sunctions, by = "Дата", all.x = TRUE) %>%
  merge(respectation_inf, by = "Дата", all.x = TRUE) %>%
  # .[1,`Скор_санкций`:=2037] %>%
  .[, (names(.)[sapply(., is.numeric) & names(.) != "Дата"]) := lapply(.SD, function(x) {
    if (any(is.na(x))) zoo::na.spline(x, na.rm = FALSE) else x}), .SDcols = names(.)[sapply(., is.numeric) & names(.) != "Дата"]] %>%
  setnames("Дата", "date")



prep_input_data <- copy(all_input_data) %>%
  setorder(date) %>%
  .[,gdp_trend := loess(gdp_level ~ as.numeric(date), span=0.3)$fitted] %>%
  .[,gdp_gap_lag1  := lag(gdp_level, 1) - lag(gdp_trend, 1)] %>%
  .[,cpi_infl_lag1 := lag(cpi_inflation, 1)] %>%
  .[,gdp_gap2  := (gdp_level - gdp_trend) / gdp_trend * 100] %>%
  .[,rate_diff := ruonia_rate - lag(ruonia_rate, 1)]

# 2. Хвостовые распределения (оценка α‑stable, Poisson, Pareto) ####
resid0 <- lm(ruonia_rate ~ cpi_inflation + gdp_gap2, data = prep_input_data) %>%
  residuals()

stable_fit <- resid0 %>% stable_fit()
tail95     <- resid0[resid0 > quantile(resid0, .95, na.rm = TRUE)]
pareto_fit <- fitdist(tail95, "pareto", start = list(shape = 1, scale = min(tail95, na.rm = TRUE)))
lambda_hat <- mean(resid0 > quantile(resid0, .95, na.rm = TRUE))

# 3. Оценка ожидаемой инфляции (RNN/Transformer для E_t[π]) ####
dt <- dt %>%
  mutate(E_pi = estimate_expectations(.))  # функция из 3_ml_expectations.R

# 4. Регрессия дрифта (расш. Тейлора) ####
dt <- dt %>%
  mutate(rate_diff = ruonia_rate - lag(ruonia_rate)) %>%
  lm(rate_diff ~ I(cpi_inflation - pi_target) + gdp_gap2 +
       volatility_index + sanctions_score + E_pi, data = .) %>%
  { dt %>% mutate(drift = predict(., newdata = dt)) }

# 5. Подгонка σ, γ, λ (MLE / variational) ####
negLogLik <- function(log_sigma, log_gamma, alpha, lambda) {
  sigma <- exp(log_sigma); gamma <- exp(log_gamma)
  ll_norm   <- dnorm(resid0,   sd = sigma,    log = TRUE)
  ll_stable <- dstable(resid0, alpha = alpha, beta = 0, scale = gamma, log = TRUE)
  ll_pois   <- dpois(resid0 > quantile(resid0, .95, na.rm = TRUE), lambda, log = TRUE)
  -sum(ll_norm + ll_stable + ll_pois, na.rm = TRUE)
}
mle_fit   <- mle(negLogLik,
                 start = list(log_sigma=log(sd(resid0)), log_gamma=log(0.2),
                              alpha=1.7, lambda=lambda_hat),
                 method = "L-BFGS-B",
                 lower  = c(-10, -10, 0.5, 0))
params    <- coef(mle_fit)
sigma_hat <- exp(params["log_sigma"])
gamma_hat <- exp(params["log_gamma"])
alpha_hat <- params["alpha"]
lambda_hat<- params["lambda"]

# 6. Симуляции (Monte Carlo траектории) ####
simulate_path <- function(drift, r0, sigma, gamma, alpha, lambda, H = 12) {
  r <- numeric(H); r[1] <- r0
  for (t in 2:H) {
    eps_g <- rnorm(1, 0, sigma)
    Nj    <- rpois(1, lambda)
    jumps <- if (Nj>0) rstable(Nj, alpha=alpha, beta=0, scale=gamma) else 0
    r[t]  <- r[t-1] + drift[t] + eps_g + sum(jumps)
  }
  r
}

paths <- replicate(1000,
                   simulate_path(dt$drift, dt$ruonia_rate[1],
                                 sigma_hat, gamma_hat, alpha_hat, lambda_hat),
                   simplify = "matrix")

# 7. Оценка симуляций (CRPS, PIT, VaR/CVaR) ####
actual <- dt$ruonia_rate[2:13]
metrics <- list(
  CRPS = apply(paths, 2, ~ crps_sample(actual, .x)) %>% colMeans(),
  PIT  = apply(paths, 2, ~ pit_sample(actual, .x))  %>% colMeans(),
  VaR95= apply(paths, 2, quantile, .95),
  CVaR95= apply(paths, 2, function(x) mean(x[x>quantile(x,.95)]))
)

# 8. Визуализации ####
qs <- paths %>% apply(2, quantile, probs = c(.1,.25,.5,.75,.9)) %>% t() %>% as.data.table()
qs[, horizon := .I]

fig <- qs %>%
  plot_ly(x = ~horizon) %>%
  add_ribbons(ymin = ~`10%`, ymax = ~`90%`, fillcolor = 'rgba(0,100,200,0.2)', line = list(color = 'transparent'), name = '10–90%') %>%
  add_ribbons(ymin = ~`25%`, ymax = ~`75%`, fillcolor = 'rgba(0,100,200,0.4)', line = list(color = 'transparent'), name = '25–75%') %>%
  add_lines(y = ~`50%`, line = list(color = 'black'), name = 'Median') %>%
  layout(title = 'Fan Chart прогнозов Ruonia (Plotly)',
         xaxis = list(title = 'Горизонт (кварталы)'),
         yaxis = list(title = 'Ставка, %'))

fig
