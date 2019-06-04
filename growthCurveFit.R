if(!require(rstan)){install.packages("rstan"); library(rstan)}
if(!require(grofit)){install.install.packages("grofit_1.1.tar.gz", repos = NULL, type="source")
; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}

rstan_options(auto_write = TRUE)
options(mc.cores = 8)

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 56, 7))

TestRun = gcFitModel(times, weightF6[1,weight_traits])
summary(TestRun$nls)
plot(TestRun)

N = nrow(weightF6)

narrow_weight = weightF6[1:N,] %>% mutate(ind = 1:N) %>% gather(trait, value, Weight_D0:Weight_D56) %>% mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% filter(!is.na(value))

stan_data = list(N = N,
                 M = nrow(narrow_weight),
                 ind = narrow_weight$ind,
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                 time = narrow_weight$times,
                 y = narrow_weight$value)

#partialPooledFit = stan(file = "fitGompertz.stan", model_name = "partial_pooled_gompertz", data = stan_data, 
#                 iter = 2000, control = list(adapt_delta = 0.99))
#saveRDS(partialPooledFit, "./Rdatas/fit_gompertz.rds")
partialPooledFit = readRDS("./Rdatas/fit_gompertz.rds")

plot(partialPooledFit, pars = c("mu_i"))
plot(partialPooledFit, pars = c("A"))
plot(partialPooledFit, pars = c("mu"))
plot(partialPooledFit, pars = c("A_0", "A_sex", "mu_0", "mu_sex"))

coefs = summary(partialPooledFit, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefs = summary(partialPooledFit, pars = c("A", "mu", "lambda"))$summary[,"mean"]

grid <- with(narrow_weight, seq(min(times), max(times), length = 100))
gompertz_mean_curve <- ddply(narrow_weight, "Sex", function(df) {
  data.frame( 
    times = grid,
    curve = gompertz(grid, A = 10 * (coefs["A_0"] + ifelse(df$Sex[1] == "M", coefs["A_sex"], 0)), 
                       mu = coefs["mu_0"] + ifelse(df$Sex[1] == "M", coefs["mu_sex"], 0), 
                       lambda = coefs["lambda_0"] + ifelse(df$Sex[1] == "M", coefs["lambda_sex"], 0))
  )
}
)
gompertz_ID_curve <- ddply(narrow_weight, .(Sex, ID), function(df) {
  id = which(weightF6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = gompertz(grid, A = all_coefs[paste0("A[", id, "]")], 
                       mu = all_coefs[paste0("mu[", id, "]")], 
                       lambda = all_coefs[paste0("lambda[", id, "]")])
  )
}
)

ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.1) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = gompertz_ID_curve, colour = "gray", alpha = 0.1) + 
  geom_line(aes(y = curve, group = 1), data = gompertz_mean_curve, colour = "red")
  

#### Logistic

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
partialPooledLogistic = readRDS("./Rdatas/fit_logistics_short.rds")

fake_data_matrix  <- partialPooledLogistic %>% 
  as.data.frame %>% 
  dplyr::select(dplyr::contains("y_sim"))
narrow_weight_sim = narrow_weight
narrow_weight_sim$value = t(fake_data_matrix[2,])

plot(partialPooledLogistic, pars = c("mu_i"))
plot(partialPooledLogistic, pars = c("A"))
plot(partialPooledLogistic, pars = c("mu"))
plot(partialPooledLogistic, pars = c("lambda_0", "lambda_sex"))
plot(partialPooledLogistic, pars = c("A_0", "A_sex"))
plot(partialPooledLogistic, pars = c("mu_0", "mu_sex"))

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
# save_plot("~/Desktop/growth-curves.png", curves_plot,
#           base_height = 6, base_aspect_ratio = 1, ncol = 2)


ggplot(narrow_weight, aes(times, value, group = ID)) + 
  geom_text(data = filter(narrow_weight, Sheep_tagID == 1213), aes(label = Sheep_tagID)) + 
  facet_wrap(~Sex) + 
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.1) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") +
  scale_x_continuous(breaks = seq(0, 60, 10)) + 
  labs(x = "Days", "Weight (g)")
