if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(MasterBayes)){install.packages("MasterBayes"); library(MasterBayes)}
if(!require(MCMCglmm)){install.packages("MCMCglmm"); library(MCMCglmm)}
if(!require(pedantics)){install.packages("pedantics"); library(pedantics)}

full_data_F6 = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex, 
         Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
         Birth_litter_size_weaning, Foster_litter_size_weaning,
         Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "F6")

full_data_Strain = read_csv("./data/Mouse phenotypes.csv") %>%
  dplyr::select(Litter_ID_new:Sex, 
                Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                Birth_litter_size_weaning, Foster_litter_size_weaning,
                Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "Strain")

full_data_F6$ID[full_data_F6$ID == 3202] = "3302"

pedigree = as.data.frame(read.csv("./data/Intercross_pedigree.csv")) %>% 
  rename(id = animal) %>% orderPed

full_data_F6 <- mutate(full_data_F6, 
                       growth_D0D14  = Weight_D14 - Weight_D0,
                       growth_D14D28 = Weight_D28 - Weight_D14,
                       growth_D28D42 = Weight_D42 - Weight_D28,
                       growth_D42D56 = Weight_D56 - Weight_D42)

                       #growth_D0D3   = Weight_D3 - Weight_D0,
                       #growth_D3D7   = Weight_D7 - Weight_D3,
                       #growth_D7D14  = Weight_D14 - Weight_D7,
                       #growth_D14D21 = Weight_D21 - Weight_D14,
                       #growth_D21D28 = Weight_D28 - Weight_D21,
                       #growth_D28D35 = Weight_D35 - Weight_D28,
                       #growth_D35D42 = Weight_D42 - Weight_D35,
                       #growth_D42D49 = Weight_D49 - Weight_D42, 
                       #growth_D49D56 = Weight_D56 - Weight_D49)

full_data_Strain <- mutate(full_data_Strain, 
                       growth_D0D14  = Weight_D14 - Weight_D0,
                       growth_D14D28 = Weight_D28 - Weight_D14,
                       growth_D28D42 = Weight_D42 - Weight_D28,
                       growth_D42D56 = Weight_D56 - Weight_D42)

                       #growth_D0D3   = Weight_D3 - Weight_D0,
                       #growth_D3D7   = Weight_D7 - Weight_D3,
                       #growth_D7D14  = Weight_D14 - Weight_D7,
                       #growth_D14D21 = Weight_D21 - Weight_D14,
                       #growth_D21D28 = Weight_D28 - Weight_D21,
                       #growth_D28D35 = Weight_D35 - Weight_D28,
                       #growth_D35D42 = Weight_D42 - Weight_D35,
                       #growth_D42D49 = Weight_D49 - Weight_D42, 
                       #growth_D49D56 = Weight_D56 - Weight_D49)

weight_traits = c("Weight_D0", "Weight_D3", "Weight_D7", "Weight_D14", "Weight_D21",
                  "Weight_D28", "Weight_D35", "Weight_D42", "Weight_D49", "Weight_D56")

#growth_traits = c("growth_D0D3", "growth_D3D7", "growth_D7D14", "growth_D14D21", "growth_D21D28",
                  #"growth_D28D35", "growth_D35D42", "growth_D42D49", "growth_D49D56")

growth_traits = c("growth_D0D14", "growth_D14D28", "growth_D28D42", "growth_D42D56") 
growthF6 = full_data_F6 %>% dplyr::select(Litter_ID_new:Sex, 
                        Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                        Birth_litter_size_weaning, Foster_litter_size_weaning, 
                        growth_D0D14:growth_D42D56, Final_weight) %>% na.omit

growthStrain = full_data_Strain %>% dplyr::select(Litter_ID_new:Sex, 
                                          Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                                          Birth_litter_size_weaning, Foster_litter_size_weaning, 
                                          growth_D0D14:growth_D42D56, Final_weight) %>% na.omit


weightF6 = full_data_F6 %>% dplyr::select(Litter_ID_new:Sex, 
                                   Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                                   Birth_litter_size_weaning, Foster_litter_size_weaning, 
                                   Weight_D0:Weight_D56)

eVec = eigen(cov(growthF6[,growth_traits]))$vectors
growthF6$fast = as.matrix(growthF6[, growth_traits]) %*% eVec[,1]
growthF6$fast[growthF6$Sex == "M"] = scale(growthF6$fast[growthF6$Sex == "M"])
growthF6$fast[growthF6$Sex == "F"] = scale(growthF6$fast[growthF6$Sex == "F"])
m_full_F6 = gather(growthF6, period, growth, growth_D0D14:growth_D42D56)

m_full_F6$period = factor(m_full_F6$period, levels = growth_traits)
ggplot(m_full_F6, aes(period, growth, group = ID)) + geom_line(alpha = 0.1)

necropsy_traits = c("Liver", "Spleen", "Kidney_total", "Heart", "Fat")

necropsyF6 = full_data_F6 %>% dplyr::select(Litter_ID_new:Sex, 
                                            Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                                            Birth_litter_size_weaning, Foster_litter_size_weaning, 
                                            Liver:Fat)
# 
# if(!require(superheat)){install.packages("superheat"); library(superheat)}
# if(!require(corrplot)){install.packages("corrplot"); library(corrplot)}
# 
# ggpairs(weightF6, mapping = aes(color = Sex), columns = weight_traits, upper = "blank")
# corrP = cor(residuals(lm(as.matrix(growthF6[,growth_traits])~growthF6$Sex)))
# gcorr_M = cor(growthF6[growthF6$Sex == "M",growth_traits])
# gcorr_F = cor(growthF6[growthF6$Sex == "F",growth_traits])
# plot_grid(superheat(gcorr_M), superheat(gcorr_F))
# colnames(gcorr) = gsub("growth_", "", growth_traits)
# png(filename = "~/Dropbox/labbio/relatorios/fapesp/2017-08-10-Doutorado-Parcial-2/corr_crescimento.png", width = 600, height = 600)
# corrplot = corrplot.mixed(gcorr, upper = "ellipse")
# dev.off()
