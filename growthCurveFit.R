library(rstan)
library(grofit)
library(plyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("./read_F6_phenotypes.R")

names(weightF6)
weightF6[,weight_traits]

times = c(0, 3, seq(7, 49, 7))

mus = sapply(1:nrow(weightF6), function(i) coef(gcFitModel(times, weightF6[i,weight_traits], control = grofit.control(model.type = "gompertz"))$nls)["mu"])
hist(unlist(mus))


TestRun = gcFitModel(times, weightF6[1,weight_traits], control = grofit.control(model.type = "gompertz"))
summary(TestRun$nls)
plot(TestRun)

N = 100
K = length(weight_traits)

stan_data = list(N = N,
                 K = K,
                 ind = rep(1:N, K),
                 sex = (as.numeric(factor(weightF6$Sex)) - 1)[1:N],
                 time = rep(times, each = N),
                 y = as.numeric(as.matrix((weightF6[1:N,weight_traits]))))

pooledFit = stan(file = "fitGrompertzPooled.stan", model_name = "pooled_gompertz", data = stan_data, 
                 iter = 2000, control = list(adapt_delta = 0.99))

summary(pooledFit)

partialPooledFit = stan(file = "fitGrompertz.stan", model_name = "partial_pooled_gompertz", data = stan_data, 
                 iter = 2000, control = list(adapt_delta = 0.99))

summary(partialPooledFit, pars = "sigma_mu")
summary(partialPooledFit, pars = c("A_0", "", "lambda"))
summary(partialPooledFit, pars = c("A", "mu", "lambda"))

plot(partialPooledFit, pars = c("A", "mu", "lambda"))
plot(partialPooledFit, pars = c("A_0", "A_sex", "mu_0", "mu_sex"))


library("plyr")
coeflines <- alply(as.matrix(coefs), 1, function(coef) {
    stat_function(fun=function(x){coef[1]*x^coef[2]}, colour="grey")
  })
