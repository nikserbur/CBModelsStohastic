# Централизованный соурс всех разработанных функций
# Стохастическое моделирование процентных ставок ЦБ
# Последний апдейт 25.12.2024

#### Подключение библиотек ####
try(source("1.R/9.libraries.r"))

#### Подключение всех функций проекта ####
try(source("1.R/1.functions/HarVolatilityFunctions.R"))
try(source("1.R/1.functions/HmmFunctions.R"))
try(source("1.R/1.functions/BoostingFunctions.R"))
try(source("1.R/1.functions/TalebR2Functions.R"))

#### Настройки соединения ####
# Чтобы не выдавал ошибку Error in curl::curl_fetch_memory(url, handle = handle) :SSL certificate problem: certificate has expired
httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))

cat("Все библиотеки и функции успешно подключены!\n")
