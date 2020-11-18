if(!require(rstan)){install.packages("rstan", dependencies = TRUE); library(rstan)}
if(!require(grofit)){install.install.packages("grofit_1.1.tar.gz", repos = NULL, type="source")
  ;library(grofit)}
if(!require(shinystan)){install.packages("shinystan"); library(shinystan)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
theme_set(theme_cowplot())
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}

if(!require(qvalue)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
  library(qvalue)
}

rstan_options(auto_write = TRUE)
options(mc.cores = 8)

times = c(0, 3, seq(7, 56, 7))

N = nrow(mice_weight$F6)

narrow_weight = mice_weight$F6[1:N,] %>% 
  mutate(ind = 1:N) %>% 
  gather(trait, value, Weight_D0:Weight_D56) %>% 
  mutate(times = as.numeric(unlist(gsub("[^0-9]", "", unlist(trait))))) %>% 
  filter(!is.na(value))

grow_models = readRDS(file = "Rdatas/no_pooling_nls_growth_fit.RDS")
mu_non_pooled = sapply(grow_models, function(x) x$parameters$mu[1])

(grow_models)

mu_nls = data.frame(ID = mice_weight$F6$ID, 
                    Sex = mice_weight$F6$Sex, 
                    mu = mu_non_pooled)

ggplot(mu_nls, aes(mu, group = Sex)) + geom_histogram(bins = 100) + facet_wrap(~Sex, ncol = 1)
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
ggplot(mu_df, aes(nP_mu, pPLVS_mu, group = Sex, shape = Sex, color = Sex)) + 
  geom_point(size = 3, alpha = 0.8) + theme_minimal() +
  labs(x = "No pooling", y = "Pooled") + geom_abline(slope = 1, intercept = 0)
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

plot(gcFitModel(times, mice_weight$F6[mice_weight$F6$ID == 3815,mice_weight$traits], 
                control = grofit.control(model.type = "logistic")))

excluded_Ids = c("4010", "4011", "4004", "3807", "3815")
mu_df_clean = mu_df %>%
  filter(!ID %in% excluded_Ids)

ggplot(mu_df_clean, aes(scale(pPLVS_mu))) + geom_histogram(bins = 100)

mu_df_clean$mu_scaled = scale(mu_df_clean$pPLVS_mu)
mu_df_clean$growth_class = mu_df_clean$mu_scaled > 0 
ggplot(mu_df_clean, aes(1,mu_scaled, color = growth_class)) + geom_violin() + geom_jitter(width = 0.1, height = 0) 

pPLVS_mu_summary = summary(partialPooledLogistiVarSigma, pars = c("mu_i"), probs = c(0.25, 0.75))$summary
pPLVS_mu_df = dplyr::select(narrow_weight, ID, ind, TPS_ID) %>% 
  unique %>% 
  arrange(ind) %>%
  mutate(mean = pPLVS_mu_summary[,'mean'],
         upper =  pPLVS_mu_summary[,'75%'],
         lower = pPLVS_mu_summary[,'25%'])

classifyMu <- function(x) {
  if(x['upper'] > 0 && x['lower'] > 0)
    return("fast")
  else if(x['upper'] < 0 && x['lower'] < 0)
    return("slow")
  else 
    return("0")
}

pPLVS_mu_df$measured = pPLVS_mu_df$TPS_ID %in% mice_crania$ed$TPS_ID

mu_df = left_join(pPLVS_mu_df, ddply(pPLVS_mu_df, .(ID), classifyMu)) %>%
  filter(measured)
ggplot(mu_df, aes(V1, mean, color = V1, group = V1)) + geom_violin() + geom_jitter(width = 0.1, height = 0) 
table(mu_df$V1)
write_csv(mu_df, "/home/MouseScans/Mouse_phenotypes/F6_growth_rates.csv")

mu_gwas = univariateGEMMA(pPLVS_mu ~ 1, data = mu_df)

if(!require(ggman)){install.packages("ggman"); library(ggman)}
qv = qvalue(mu_gwas$p_wald, fdr.level = 0.01)
table(qv$significant)
mu_gwas = mu_gwas %>%
  mutate(qvalue = qvalue(p_wald)$qvalue)

mu.gwas.sig <- mu_gwas[qv$significant,]

p1 <- ggman(mu_gwas, snp = "rs", pvalue = "qvalue", bp = "ps", pointSize = 1, chrom = "chr", sigLine = 0)
ggmanLabel(p1, labelDfm = as.data.frame(mu.gwas.sig), snp = "rs", type = "label", label = "rs")
