if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}

rstan_options(auto_write = TRUE)
options(mc.cores = 8)

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 56, 7))
N = nrow(weightF6)

narrow_weight = weightF6[1:N,] %>% mutate(ind = 1:N) %>% gather(trait, value, Weight_D0:Weight_D56) %>% mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% filter(!is.na(value))

stan_data = list(N = N,
                 M = nrow(narrow_weight),
                 ind = narrow_weight$ind,
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                 time = narrow_weight$times,
                 y = narrow_weight$value, 
                 run_estimation = 1)
# attach(stan_data)
# stan_rdump(names(stan_data), file = "logistc_fit_data.R")
partialPooledLogistic = stan(file = "fitLogistic.stan", 
                             model_name = "partial_pooled_logistic", data = stan_data, 
                             iter = 2000, chains = 8, control = list(adapt_delta = 0.99))
saveRDS(partialPooledLogistic, "./Rdatas/fit_logistics.rds")
partialPooledLogistic = readRDS("./Rdatas/fit_logistics.rds")
