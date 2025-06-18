# Централизованный соурс всех разработанных функций
# Последний апдейт 09.04.2023

try(source("1.R/9.libraries.r"))

# Чтобы не выдавал ошибку Error in curl::curl_fetch_memory(url, handle = handle) :SSL certificate problem: certificate has expired
httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))

print("Sourced")
