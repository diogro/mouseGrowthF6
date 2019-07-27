if(!require(rstan)){install.packages("rstan", dependencies = TRUE); library(rstan)}
if(!require(grofit)){install.install.packages("grofit_1.1.tar.gz", repos = NULL, type="source")
  ; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 56, 7))

N = nrow(weightF6)

narrow_weight = weightF6[1:N,] %>% 
  mutate(ind = 1:N) %>% 
  gather(trait, value, Weight_D0:Weight_D56) %>% 
  mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% 
  filter(!is.na(value))

#### Logistic no pooling

TestRun = gcFitModel(times, weightF6[10,weight_traits])
summary(TestRun$nls)
plot(TestRun)

grow_models = vector("list", nrow(weightF6))
for(i in 1:nrow(weightF6)){
  grow_models[[i]] = gcFitModel(times, weightF6[i,weight_traits], control = grofit.control(model.type = "logistic"))
}
saveRDS(grow_models, file = "Rdatas/no_pooling_nls_growth_fit.RDS")
mu_non_pooled = sapply(grow_models, function(x) x$parameters$mu[1])
A_non_pooled = sapply(grow_models, function(x) x$parameters$A[1])
lambda_non_pooled = sapply(grow_models, function(x) x$parameters$lambda[1])

hist(mu_non_pooled, breaks = 100)
hist(A_non_pooled, breaks = 100)
hist(lambda_non_pooled, breaks = 100)
mu_nls = data.frame(ID = weightF6$ID, Sex = weightF6$Sex, mu = mu_non_pooled, A = A_non_pooled, lambda = lambda_non_pooled)
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
fitLogisticPooled = stan(file = "fitLogisticPooled.stan", data = stan_data_pooled, 
                         chains = 4, warmup = 4000, iter = 5000, control = list(adapt_delta = 0.999, max_treedepth = 11))

coefsLogistic = summary(fitLogisticPooled, pars = c("A", "A_s", "mu", "mu_s", "lambda", "lambda_s"))$summary[,"mean"]
all_coefsLogistic = summary(fitLogisticPooled, pars = c("A", "mu", "lambda"))$summary[,"mean"]

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
                                sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                                time = narrow_weight$times,
                                y = narrow_weight$value, 
                                run_estimation = 1)
partialPooledLogistic = stan(file = "fitLogistic.stan", 
                             model_name = "partial_pooled_logistic", data = stan_data_partial_pooled, 
                             iter = 10000, warmup = 8000, chains = 4, 
                             control = list(adapt_delta = 0.999, max_treedepth = 11))
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
  id = which(weightF6$ID == df$ID[1])
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
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
                 time = narrow_weight$times,
                 n_times = 10,
                 time_int = as.numeric(factor(narrow_weight$times, levels = times)),
                 y = narrow_weight$value, 
                 run_estimation = 1)

#partialPooledLogistiVarSigma = stan(file = "fitLogisticVariableSigma.stan", 
                                    #model_name = "partial_pooled_logistic", data = stan_data, 
                                    #iter = 2000, chains = 4, control = list(adapt_delta = 0.999))
#saveRDS(partialPooledLogistiVarSigma, "./Rdatas/fit_logisticsVarSigma.rds")
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
  id = which(weightF6$ID == df$ID[1])
  data.frame( 
    times = grid,
    curve = logistic(grid, A = all_coefsLogistic[paste0("A[", id, "]")], 
                     mu = all_coefsLogistic[paste0("mu[", id, "]")], 
                     lambda = all_coefsLogistic[paste0("lambda[", id, "]")])
  )
}
)

#### Logistic partial pooled varying sigmas animal model

weightF6

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
                 sex = (as.numeric(factor(weightF6$Sex[1:N])) - 1),
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

mu_i = summary(partialPooledLogistiVarSigmaAnimal, pars = c("mu_i"))$summary[,"mean"]
mus = data.frame(mu_vS = mu_i, ind = 1:length(mu_i))
hist(mu_i, breaks = 100)
narrow_weight$ind

unique(left_join(mus, narrow_weight[,c("ind", "ID")], by = "ind"))
