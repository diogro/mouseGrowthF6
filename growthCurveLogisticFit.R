if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(grofit)){install.packages("grofit_1.1.1-1.tar.gz", repos = NULL, type="source")
  ; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}

rstan_options(auto_write = TRUE)
options(mc.cores = 2)

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
                             iter = 1000, chains = 1, control = list(adapt_delta = 0.99))
#saveRDS(partialPooledLogistic, "./Rdatas/fit_logistics.rds")
#partialPooledLogistic = readRDS("./Rdatas/fit_logistics.rds")

fake_data_matrix  <- partialPooledLogistic %>% 
  as.data.frame %>% 
  dplyr::select(dplyr::contains("y_sim"))
narrow_weight_sim = narrow_weight
narrow_weight_sim$value = t(fake_data_matrix[1,])

coefsLogistic = summary(partialPooledLogistic, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefsLogistic = summary(partialPooledLogistic, pars = c("A", "mu", "lambda"))$summary[,"mean"]

grid <- with(narrow_weight, seq(min(times), max(times), length = 100))
logistic_mean_curve <- ddply(narrow_weight, "Sex", function(df) {
  data.frame( 
    times = grid,
    curve = logistic(grid, A = 10 * (coefsLogistic["A_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["A_sex"], 0)), 
                     mu = coefsLogistic["mu_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["mu_sex"], 0), 
                     lambda = coefsLogistic["lambda_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["lambda_sex"], 0))
  )
}
)
logistic_ID_curve <- ddply(narrow_weight, .(Sex, ID), function(df) {
  id = which(weightF6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = logistic(grid, A = all_coefsLogistic[paste0("A[", id, "]")], 
                     mu = all_coefsLogistic[paste0("mu[", id, "]")], 
                     lambda = all_coefsLogistic[paste0("lambda[", id, "]")])
  )
}
)

sim_data_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + 
  geom_jitter(alpha = 0.3, size = 0.6) + 
  geom_jitter(data = narrow_weight_sim, alpha = 0.3, size = 0.6, color = "red") + 
  facet_wrap(~Sex) +
  labs(x = "Dias", y = "Peso (g)")

curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.15) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")
