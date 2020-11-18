if(!require(rstan)){install.packages("rstan", dependencies = TRUE); library(rstan)}
if(!require(cmdstanr)){install.packages("cmdstanr"); library(cmdstanr)}
if(!require(grofit)){install.packages("./grofit_1.1.1-1.tar.gz", repos = NULL, type="source")
  ; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

library(AtchleyMice)

times = c(0, 3, seq(7, 56, 7))

N = nrow(mice_weight$F6)

narrow_weight = mice_weight$F6[1:N,] %>% 
  mutate(ind = 1:N) %>% 
  gather(trait, value, Weight_D0:Weight_D56) %>% 
  mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% 
  filter(!is.na(value))

mice_weight$traits
#### Logistic no pooling

TestRun = gcFitModel(times, mice_weight$F6[10, mice_weight$traits])
summary(TestRun$nls)
plot(TestRun)

grow_models = vector("list", nrow(mice_weight$F6))
for(i in 1:nrow(mice_weight$F6)){
  grow_models[[i]] = gcFitModel(times, mice_weight$F6[i,mice_weight$traits], control = grofit.control(model.type = "logistic"))
}
saveRDS(grow_models, file = "Rdatas/no_pooling_nls_growth_fit.RDS")
mu_non_pooled = sapply(grow_models, function(x) x$parameters$mu[1])
A_non_pooled = sapply(grow_models, function(x) x$parameters$A[1])
lambda_non_pooled = sapply(grow_models, function(x) x$parameters$lambda[1])

hist(mu_non_pooled, breaks = 100)
hist(A_non_pooled, breaks = 100)
hist(lambda_non_pooled, breaks = 100)
mu_nls = data.frame(ID = mice_weight$F6$ID, Sex = mice_weight$F6$Sex, mu = mu_non_pooled, A = A_non_pooled, lambda = lambda_non_pooled)
ggplot(mu_nls, aes(mu, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)
mu_nls$mu_res[!is.na(mu_nls$mu)] = residuals(lm(mu~Sex, data = mu_nls))
mu_nls$A_res[!is.na(mu_nls$mu)] = residuals(lm(A~Sex, data = mu_nls))
mu_nls$lambda_res[!is.na(mu_nls$mu)] = residuals(lm(lambda~Sex, data = mu_nls))
ggplot(mu_nls, aes(mu_res, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)

mean(mu_nls$mu, na.rm = T)
sd(mu_nls$mu_res, na.rm = T)

weight_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_text(data = filter(narrow_weight, ID == 4167), aes(label = ID), color = "red") + 
  labs(x = "Dias", y = "Peso (g)")

##### Logistic complete pooling

stan_data_pooled = list(M = nrow(narrow_weight),
                        y = narrow_weight$value,
                        sex = as.numeric(factor(narrow_weight$Sex, levels = c("F", "M"))) - 1,
                        time = narrow_weight$times)
LogisticPooled_model = cmdstan_model(stan_file = "fitLogisticPooled.stan")
fitLogisticPooled = LogisticPooled_model$sample(data = stan_data_pooled, 
                                                chains = 4, 
                                                parallel_chains = 4,
                                                iter_warmup = 4000, 
                                                iter_sampling = 1000, 
                                                adapt_delta = 0.999, 
                                                max_treedepth = 11)

fitLogisticPooled_draws <- rstan::read_stan_csv(fitLogisticPooled$output_files())

coefsLogistic = summary(fitLogisticPooled_draws, pars = c("A", "A_s", "mu", "mu_s", "lambda", "lambda_s"))$summary[,"mean"]
all_coefsLogistic = summary(fitLogisticPooled_draws, pars = c("A", "mu", "lambda"))$summary[,"mean"]

grid <- with(narrow_weight, seq(min(times), max(times), length = 100))
logistic_mean_curve <- ddply(narrow_weight, "Sex", function(df) {
  data.frame( 
    times = grid,
    curve = logistic(grid, A = 10 * (coefsLogistic["A"] + ifelse(df$Sex[1] == "M", coefsLogistic["A_s"], 0)), 
                     mu = coefsLogistic["mu"] + ifelse(df$Sex[1] == "M", coefsLogistic["mu_s"], 0), 
                     lambda = coefsLogistic["lambda"] + ifelse(df$Sex[1] == "M", coefsLogistic["lambda_s"], 0))
  )
}
)

curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

##### Logistic partial pooling

stan_data_partial_pooled = list(N = N,
                                M = nrow(narrow_weight),
                                ind = narrow_weight$ind,
                                sex = (as.numeric(factor(mice_weight$F6$Sex[1:N])) - 1),
                                time = narrow_weight$times,
                                y = narrow_weight$value, 
                                run_estimation = 1)
Logistic_model = cmdstan_model(stan_file = "fitLogistic.stan")
fitLogistic = Logistic_model$sample(data = stan_data_partial_pooled, 
                                    chains = 4, 
                                    parallel_chains = 4,
                                    iter_warmup = 1000, 
                                    iter_sampling = 1000, 
                                    adapt_delta = 0.999, 
                                    max_treedepth = 11)
partialPooledLogistic = rstan::read_stan_csv(fitLogistic$output_files())

saveRDS(partialPooledLogistic, "./Rdatas/fit_logistics.rds")
partialPooledLogistic = readRDS("./Rdatas/fit_logistics.rds")

fake_data_matrix  <- partialPooledLogistic %>% 
  as.data.frame %>% 
  dplyr::select(dplyr::contains("y_sim"))
narrow_weight_sim = narrow_weight
narrow_weight_sim$value = t(fake_data_matrix[2,])

coefsLogistic = summary(partialPooledLogistic, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefsLogistic = summary(partialPooledLogistic, pars = c("A", "mu", "lambda"))$summary[,"mean"]

mu_i = summary(partialPooledLogistic, pars = c("mu_i"))$summary[,"mean"]
hist(mu_i, breaks = 100)

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
  id = which(mice_weight$F6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = logistic(grid, A = all_coefsLogistic[paste0("A[", id, "]")], 
                     mu = all_coefsLogistic[paste0("mu[", id, "]")], 
                     lambda = all_coefsLogistic[paste0("lambda[", id, "]")])
  )
}
)

curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.15) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

sim_data_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + 
  geom_jitter(alpha = 0.3, size = 0.6) + 
  geom_jitter(data = narrow_weight_sim, alpha = 0.3, size = 0.6, color = "red") + 
  facet_wrap(~Sex) +
  labs(x = "Dias", y = "Peso (g)")


#### Logistic partial pooled varying sigmas

stan_data = list(N = N,
                 M = nrow(narrow_weight),
                 ind = narrow_weight$ind,
                 sex = (as.numeric(factor(mice_weight$F6$Sex[1:N])) - 1),
                 time = narrow_weight$times,
                 n_times = 10,
                 time_int = as.numeric(factor(narrow_weight$times, levels = times)),
                 y = narrow_weight$value, 
                 run_estimation = 1)
LogisticVariableSigma = cmdstan_model(stan_file = "fitLogisticVariableSigma.stan")
fitLogisticVariableSigma = LogisticVariableSigma$sample(data = stan_data, 
                                                        chains = 4, 
                                                        parallel_chains = 4,
                                                        iter_warmup = 1000, 
                                                        iter_sampling = 1000, 
                                                        adapt_delta = 0.999, 
                                                        max_treedepth = 11)
partialPooledLogistiVarSigma = rstan::read_stan_csv(fitLogisticVariableSigma$output_files())
saveRDS(partialPooledLogistiVarSigma, "./Rdatas/fit_logisticsVarSigma.rds")
partialPooledLogistiVarSigma = readRDS("./Rdatas/fit_logisticsVarSigma.rds")

fake_data_matrix  <- partialPooledLogistiVarSigma %>% 
  as.data.frame %>% 
  dplyr::select(dplyr::contains("y_sim"))
narrow_weight_sim = narrow_weight
narrow_weight_sim$value = t(fake_data_matrix[1,])

coefsLogistic = summary(partialPooledLogistiVarSigma, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefsLogistic = summary(partialPooledLogistiVarSigma, pars = c("A", "mu", "lambda"))$summary[,"mean"]

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
  id = which(mice_weight$F6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = logistic(grid, A = all_coefsLogistic[paste0("A[", id, "]")], 
                     mu = all_coefsLogistic[paste0("mu[", id, "]")], 
                     lambda = all_coefsLogistic[paste0("lambda[", id, "]")])
  )
}
)

curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.15) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

sim_data_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + 
  geom_jitter(alpha = 0.3, size = 0.6) + 
  geom_jitter(data = narrow_weight_sim, alpha = 0.3, size = 0.6, color = "red") + 
  facet_wrap(~Sex) +
  labs(x = "Dias", y = "Peso (g)")





#### Logistic partial pooled varying sigmas animal model

mice_weight$F6

F6_ids = unique(narrow_weight$ID)
ginverse = inverseA(pedigree, F6_ids)
Ainv = ginverse$Ainv
#A = solve(Ainv + 1e-3*diag(nrow(Ainv)))
A = solve(Ainv)
dimnames(A) = list(rownames(Ainv), rownames(Ainv))
if(isSymmetric(A)){A = (A + t(A))/2 
} else stop("A not symmetric")
A[1,2]

stan_data_animal = list(N = N,
                 M = nrow(narrow_weight),
                 ind = narrow_weight$ind,
                 sex = (as.numeric(factor(mice_weight$F6$Sex[1:N])) - 1),
                 time = narrow_weight$times,
                 n_times = 10,
                 time_int = as.numeric(factor(narrow_weight$times, levels = times)),
                 y = narrow_weight$value, 
                 R = as.matrix(A),
                 run_estimation = 1)

partialPooledLogistiVarSigmaAnimal = stan(file = "./fitLogisticVariableSigmaAnimal.stan", 
                                    model_name = "partial_pooled_logistic", data = stan_data_animal, 
                                    iter = 2000, chains = 1, control = list(adapt_delta = 0.999))
#saveRDS(partialPooledLogistiVarSigmaAnimal, "./Rdatas/fit_logisticsVarSigmaAnimal.rds")
#partialPooledLogistiVarSigmaAnimal = readRDS("./Rdatas/fit_logisticsVarSigmaAnimal.rds")

fake_data_matrix  <- partialPooledLogistiVarSigmaAnimal %>% 
  as.data.frame %>% 
  dplyr::select(dplyr::contains("y_sim"))
narrow_weight_sim = narrow_weight
narrow_weight_sim$value = t(fake_data_matrix[2,])

coefsLogistic = summary(partialPooledLogistiVarSigmaAnimal, pars = c("A_0", "A_sex", "mu_0", "mu_sex", "lambda_0", "lambda_sex"))$summary[,"mean"]
all_coefsLogistic = summary(partialPooledLogistiVarSigmaAnimal, pars = c("A", "mu", "lambda"))$summary[,"mean"]

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
  id = which(mice_weight$F6$ID == df$ID[1])
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

mu_i = summary(partialPooledLogistiVarSigmaAnimal, pars = c("mu_i"))$summary[,"mean"]
mus = data.frame(mu_vS = mu_i, ind = 1:length(mu_i))
hist(mu_i, breaks = 100)
narrow_weight$ind

unique(left_join(mus, narrow_weight[,c("ind", "ID")], by = "ind"))
