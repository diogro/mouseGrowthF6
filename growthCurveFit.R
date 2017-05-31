library(rstan)
if(!require(grofit)){install.packages("grofit"); library(grofit)}
library(plyr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 56, 7))

TestRun = gcFitModel(times, weightF6[1,weight_traits])
summary(TestRun$nls)
plot(TestRun)

N = nrow(weightF6)

wide_weight = weightF6[1:N,] %>% mutate(ind = 1:N) %>% gather(trait, value, Weight_D0:Weight_D56) %>% mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% filter(!is.na(value))

stan_data = list(N = N,
                 M = nrow(wide_weight),
                 ind = wide_weight$ind,
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                 time = wide_weight$times,
                 y = wide_weight$value)

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

grid <- with(wide_weight, seq(min(times), max(times), length = 100))
gompertz_mean_curve <- ddply(wide_weight, "Sex", function(df) {
  data.frame( 
    times = grid,
    curve = gompertz(grid, A = 10 * (coefs["A_0"] + ifelse(df$Sex[1] == "M", coefs["A_sex"], 0)), 
                       mu = coefs["mu_0"] + ifelse(df$Sex[1] == "M", coefs["mu_sex"], 0), 
                       lambda = coefs["lambda_0"] + ifelse(df$Sex[1] == "M", coefs["lambda_sex"], 0))
  )
}
)
gompertz_ID_curve <- ddply(wide_weight, .(Sex, ID), function(df) {
  id = which(weightF6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = gompertz(grid, A = all_coefs[paste0("A[", id, "]")], 
                       mu = all_coefs[paste0("mu[", id, "]")], 
                       lambda = all_coefs[paste0("lambda[", id, "]")])
  )
}
)

ggplot(wide_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.1) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = gompertz_ID_curve, colour = "gray", alpha = 0.1) + 
  geom_line(aes(y = curve, group = 1), data = gompertz_mean_curve, colour = "red")
  

#### Logistic

N = nrow(weightF6)

wide_weight = weightF6[1:N,] %>% mutate(ind = 1:N) %>% gather(trait, value, Weight_D0:Weight_D56) %>% mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% filter(!is.na(value))

stan_data = list(N = N,
                 M = nrow(wide_weight),
                 ind = wide_weight$ind,
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                 time = wide_weight$times,
                 y = wide_weight$value)

#partialPooledLogistic = stan(file = "fitLogistic.stan", model_name = "partial_pooled_logistic", data = stan_data, 
#                        iter = 2000, control = list(adapt_delta = 0.99))
#saveRDS(partialPooledLogistic, "./Rdatas/fit_logistics.rds")
partialPooledLogistic = readRDS("./Rdatas/fit_logistics.rds")

plot(partialPooledLogistic, pars = c("mu_i"))
plot(partialPooledLogistic, pars = c("A"))
plot(partialPooledLogistic, pars = c("mu"))
plot(partialPooledLogistic, pars = c("lambda_0", "lambda_sex"))
plot(partialPooledLogistic, pars = c("A_0", "A_sex"))
plot(partialPooledLogistic, pars = c("mu_0", "mu_sex"))

coefsLogistic = summary(partialPooledLogistic, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefsLogistic = summary(partialPooledLogistic, pars = c("A", "mu", "lambda"))$summary[,"mean"]

grid <- with(wide_weight, seq(min(times), max(times), length = 100))
logistic_mean_curve <- ddply(wide_weight, "Sex", function(df) {
  data.frame( 
    times = grid,
    curve = logistic(grid, A = 10 * (coefsLogistic["A_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["A_sex"], 0)), 
                     mu = coefsLogistic["mu_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["mu_sex"], 0), 
                     lambda = coefsLogistic["lambda_0"] + ifelse(df$Sex[1] == "M", coefsLogistic["lambda_sex"], 0))
  )
}
)
logistic_ID_curve <- ddply(wide_weight, .(Sex, ID), function(df) {
  id = which(weightF6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = logistic(grid, A = all_coefsLogistic[paste0("A[", id, "]")], 
                     mu = all_coefsLogistic[paste0("mu[", id, "]")], 
                     lambda = all_coefsLogistic[paste0("lambda[", id, "]")])
  )
}
)

ggplot(wide_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.1) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.1) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") 


ggplot(wide_weight, aes(times, value, group = ID)) + geom_text(data = filter(wide_weight, ID == 4183), aes(label = ID)) + facet_wrap(~Sex) + geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.1) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") 
