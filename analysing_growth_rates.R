if(!require(rstan)){install.packages("rstan", dependencies = TRUE); library(rstan)}
if(!require(grofit)){install.install.packages("grofit_1.1.tar.gz", repos = NULL, type="source")
  ; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())

rstan_options(auto_write = TRUE)
options(mc.cores = 8)

source("./read_F6_phenotypes.R")

times = c(0, 3, seq(7, 56, 7))

N = nrow(weightF6)

narrow_weight = weightF6[1:N,] %>% 
  mutate(ind = 1:N) %>% 
  gather(trait, value, Weight_D0:Weight_D56) %>% 
  mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% 
  filter(!is.na(value))

grow_models = readRDS(file = "Rdatas/no_pooling_nls_growth_fit.RDS")
mu_non_pooled = sapply(grow_models, function(x) x$parameters$mu[1])

mu_nls = data.frame(ID = weightF6$ID, 
                    Sex = weightF6$Sex, 
                    mu = mu_non_pooled)

ggplot(mu_nls, aes(mu, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)
mu_nls$mu_res[!is.na(mu_nls$mu)] = residuals(lm(mu~Sex, data = mu_nls))
ggplot(mu_nls, aes(mu_res, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)

weight_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_text(data = filter(narrow_weight, ID == 4167), aes(label = ID), color = "red") + 
  labs(x = "Dias", y = "Peso (g)")
