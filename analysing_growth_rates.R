if(!require(rstan)){install.packages("rstan", dependencies = TRUE); library(rstan)}
if(!require(grofit)){install.install.packages("grofit_1.1.tar.gz", repos = NULL, type="source")
  ; library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}

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

(grow_models)

mu_nls = data.frame(ID = weightF6$ID, 
                    Sex = weightF6$Sex, 
                    mu = mu_non_pooled)

mu_nls_plot = ggplot(mu_nls, aes(mu, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)
save_plot("/home/MouseScans/figures/mu_nls.png", mu_nls_plot, base_height = 5)
mu_nls$nP_mu[!is.na(mu_nls$mu)] = residuals(lm(mu~Sex, data = mu_nls))
ggplot(mu_nls, aes(nP_mu, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)

weight_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_text(data = filter(narrow_weight, ID == 4167), aes(label = ID), color = "red") + 
  labs(x = "Dias", y = "Peso (g)")

partialPooledLogistiVarSigma = readRDS("./Rdatas/fit_logisticsVarSigma.rds")

names(partialPooledLogistiVarSigma)
pPLVS_mu_summary = summary(partialPooledLogistiVarSigma, pars = c("mu_i"))$summary
pPLVS_mu = summary(partialPooledLogistiVarSigma, pars = c("mu_i"))$summary[,"mean"]
pPLVS_mu_df = dplyr::select(narrow_weight, ID, ind) %>% unique %>% mutate(pPLVS_mu = pPLVS_mu)
mu_df = left_join(mu_nls, 
                  dplyr::select(pPLVS_mu_df, ID, pPLVS_mu), 
                  by = "ID")
ggplot(mu_df, aes(nP_mu, pPLVS_mu, group = Sex, shape = Sex, color = Sex)) + geom_point(size = 3, alpha = 0.8) + theme_fivethirtyeight() + geom_abline(slope = 1, intercept = 0)
summary(lm(pPLVS_mu ~ nP_mu, data = mu_df))

ggplot(mu_df, aes(scale(pPLVS_mu))) + geom_histogram(bins = 100)

extreme_IDs = mu_df %>%
  filter(pPLVS_mu > quantile(pPLVS_mu, 0.99) | pPLVS_mu < quantile(pPLVS_mu, 0.01)) %>%
  arrange(desc(pPLVS_mu))

weight_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_text(data = filter(narrow_weight, ID %in% extreme_IDs$ID), aes(label = ID), color = "red") + 
  labs(x = "Dias", y = "Peso (g)")

weight_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_text(data = filter(narrow_weight, ID == 3199), aes(label = ID), color = "red") + 
  labs(x = "Dias", y = "Peso (g)")


curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = filter(logistic_ID_curve, ID == 3815), colour = "gray", alpha = 1) + 
  geom_text(data = filter(narrow_weight, ID == 3815), aes(label = ID), color = "red") + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

plot(gcFitModel(times, weightF6[weightF6$ID == 3815,weight_traits], control = grofit.control(model.type = "logistic")))

excluded_Ids = c("4010", "4011", "4004", "3807", "3815")
mu_df_clean = mu_df %>%
  filter(!ID %in% excluded_Ids)

ggplot(mu_df_clean, aes(scale(pPLVS_mu))) + geom_histogram(bins = 100)

mu_df_clean$mu_scaled = scale(mu_df_clean$pPLVS_mu)
mu_df_clean$growth_class = mu_df_clean$mu_scaled > 0 
ggplot(mu_df_clean, aes(1,mu_scaled, color = growth_class)) + geom_violin() + geom_jitter(width = 0.1, height = 0) 


pPLVS_mu_summary = summary(partialPooledLogistiVarSigma, pars = c("mu_i"), probs = c(0.20, 0.80))$summary
pPLVS_mu_df = dplyr::select(narrow_weight, ID, ind, Sheep_tagID) %>% 
  unique %>% 
  arrange(ind) %>%
  mutate(mean = pPLVS_mu_summary[,'mean'],
         upper =  pPLVS_mu_summary[,'80%'],
         lower = pPLVS_mu_summary[,'20%'])

classifyMu <- function(x) {
  if(x['upper'] > 0 && x['lower'] > 0)
    return("fast")
  else if(x['upper'] < 0 && x['lower'] < 0)
    return("slow")
  else 
    return("0")
}

tps_folder = "/home/MouseScans/Dicom/tps"
measured = gsub(".tps", "", list.files(tps_folder, pattern = ".\\.tps", full.names = FALSE))
measured = gsub('-.*', "", measured, perl = TRUE)
measured_twice = names(table(measured)[table(measured) >= 2])
measured = unique(measured)

ddply(pPLVS_mu_df, .(ID), classifyMu)

write_csv(left_join(pPLVS_mu_df, ddply(pPLVS_mu_df, .(ID), classifyMu), by = "ID"), 
          path = "/home/MouseScans/Mouse_phenotypes/F6_growth_rates.csv")
pPLVS_mu_df$measured = pPLVS_mu_df$Sheep_tagID %in% measured

mu_df = left_join(pPLVS_mu_df, ddply(pPLVS_mu_df, .(ID), classifyMu))
mus_plot = ggplot(mu_df, aes(1, mean, color = V1)) + geom_jitter(width = 0.1, height = 0) + 
  scale_color_manual(values = c("gray","#F8766D","#00BFC4"), label = c("0", "Fast", "Slow")) +
  labs(x = "Indivíduo", y = "Desvio da taxa de crescimento", color = "Classe") + theme(axis.text.x = element_blank())
table(mu_df$V1)

logistic_ID_curve_class = inner_join(mu_df, logistic_ID_curve, by = "ID")

curves_plot_w = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  #geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.3) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

curves_plot = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = logistic_ID_curve, colour = "gray", alpha = 0.3) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "red") + 
  labs(x = "Dias", y = "Peso (g)")

curves_plot_classified = ggplot(narrow_weight, aes(times, value, group = ID)) + geom_jitter(alpha = 0.3, size = 0.6) + facet_wrap(~Sex) +
  geom_line(aes(y = curve), data = filter(logistic_ID_curve_class, V1 == "fast"), colour = "#F8766D", alpha = 0.3) + 
  geom_line(aes(y = curve), data = filter(logistic_ID_curve_class, V1 == "slow"), colour = "#00BFC4", alpha = 0.3) + 
  #geom_line(aes(y = curve), data = filter(logistic_ID_curve_class, V1 == "0"), colour = "gray", alpha = 0.15) + 
  geom_line(aes(y = curve, group = 1), data = logistic_mean_curve, colour = "black") + 
  labs(x = "Dias", y = "Peso (g)")

save_plot("/home/MouseScans/figures/taxas_mu.png", mus_plot, base_height = 6)

  save_plot("/home/MouseScans/figures/growth_curves.png", curves_plot, base_height = 6)
  save_plot("/home/MouseScans/figures/growth_curves_weight.png", curves_plot_w, base_height = 6)
  
  save_plot("/home/MouseScans/figures/growth_curves_class.png", curves_plot_classified, base_height = 6)
  
write_csv(mu_df, "/home/MouseScans/Mouse_phenotypes/F6_growth_rates.csv")
