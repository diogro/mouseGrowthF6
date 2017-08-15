source("./read_F6_phenotypes.R")

if(!require(ggjoy)){install.packages("ggjoy"); library(ggjoy)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}

necropsy_plot = necropsyF6 %>%
  dplyr::select(-Kidney_L, -Kidney_R) %>%
  rename(Rins = Kidney_total,
         "Fígado" = Liver,
         Baço = Spleen,
         Coração = Heart,
         Gordura = Fat,
         Sexo = Sex) %>%
  gather(trait, value, Fígado:Gordura) %>%
  ggplot(aes(value, trait, fill = Sexo)) + geom_joy(alpha = 0.7) + 
  scale_fill_manual(values = viridis(2)) + 
  labs(x = "Peso (g)", y = "") + scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5))
save_plot("~/Dropbox/labbio/relatorios/fapesp/2017-08-10-Doutorado-Parcial-2/necropsy_dist.png", necropsy_plot,
          base_width = 6, base_aspect_ratio = 1)
