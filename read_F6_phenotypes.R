library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


full_data_F6 = read_csv("./data/Mouse phenotypes.csv") %>%
  select(Litter_ID_new:Sex, 
         Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
         Birth_litter_size_weaning, Foster_litter_size_weaning,
         Weight_D0:Weight_D70, Final_weight, Liver:Fat) %>%
  filter(Gen == "F6")

full_data_F6[full_data_F6$ID == 4033,"Weight_D42"] = 26.94

full_data_F6 = full_data_F6 %>%
  mutate(
    growth_D3D0   = Weight_D3 - Weight_D0,
    growth_D7D3   = Weight_D7 - Weight_D3,
    growth_D14D7  = Weight_D14 - Weight_D7,
    growth_D21D14 = Weight_D21 - Weight_D14,
    growth_D28D21 = Weight_D28 - Weight_D21,
    growth_D35D28 = Weight_D35 - Weight_D28,
    growth_D42D35 = Weight_D42 - Weight_D35,
    growth_D49D42 = Weight_D49 - Weight_D42)

growth_traits = c("growth_D3D0", "growth_D7D3", "growth_D14D7", "growth_D21D14", "growth_D28D21",
                  "growth_D35D28", "growth_D42D35", "growth_D49D42")

filter(full_data_F6, growth_D42D35 > 100)


growthF6 = full_data_F6 %>% select(Litter_ID_new:Sex, 
                        Gen, Pat_ID, Mat_ID, Nurse_ID, Litter_size_birth, 
                        Birth_litter_size_weaning, Foster_litter_size_weaning, growth_D3D0:growth_D49D42, Final_weight) %>% na.omit

cor(growthF6[,growth_traits])
eVec = eigen(cov(growthF6[,growth_traits]))$vectors
growthF6$fast = as.matrix(growthF6[, growth_traits]) %*% eVec[,1]
m_full_F6 = gather(growthF6, variable, value, growth_D3D0:growth_D49D42)

cor(growthF6$fast, growthF6$Final_weight)

m_full_F6$variable = factor(m_full_F6$variable, levels = growth_traits)
m_full_F6$growth = c(3, 7, 14, 21, 28, 35, 42, 49)[as.numeric(m_full_F6$variable)]

ggplot(m_full_F6, aes(variable, value, group = variable)) + geom_boxplot()
ggplot(m_full_F6, aes(growth, value, group = ID, color = fast)) + geom_line(alpha = 0.1) + scale_color_viridis()

