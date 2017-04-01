library(rstan)
library(grofit)
library(plyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 49, 7))

mus = sapply(1:nrow(weightF6), function(i) coef(gcFitModel(times, weightF6[i,weight_traits], control = grofit.control(model.type = "gompertz"))$nls)["mu"])
hist(unlist(mus))

TestRun = gcFitModel(times, weightF6[1,weight_traits], control = grofit.control(model.type = "gompertz"))
summary(TestRun$nls)
plot(TestRun)

N = nrow(weightF6)
K = length(weight_traits)

stan_data = list(N = N,
                 K = K,
                 ind = rep(1:N, K),
                 sex = (as.numeric(factor(weightF6$Sex)) - 1)[1:N],
                 time = rep(times, each = N),
                 y = as.numeric(as.matrix((weightF6[1:N,weight_traits]))))

pooledFit = stan(file = "fitGrompertzPooled.stan", model_name = "pooled_gompertz", data = stan_data, 
                 iter = 2000, control = list(adapt_delta = 0.9))

summary(pooledFit)

partialPooledFit = stan(file = "fitGrompertz.stan", model_name = "partial_pooled_gompertz", data = stan_data, 
                 iter = 2000, control = list(adapt_delta = 0.9))

summary(partialPooledFit, pars = "sigma_mu")
summary(partialPooledFit, pars = c("A_0"))
summary(partialPooledFit, pars = c("A", "mu", "lambda"))

plot(partialPooledFit, pars = c("mu_i"))
plot(partialPooledFit, pars = c("A"))
plot(partialPooledFit, pars = c("mu_0", "mu_sex"))

coefs = summary(partialPooledFit, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefs = summary(partialPooledFit, pars = c("A", "mu", "lambda"))$summary[,"mean"]

wide_weight = weightF6[1:N,] %>% gather(trait, value, Weight_D0:Weight_D49) %>% mutate(times = rep(times, each = N))
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
  