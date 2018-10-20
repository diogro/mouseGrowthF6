
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

mean_growth_strain = ddply(growthStrain, .(Strain), numcolwise(mean))
mean_growth_F6 = ddply(growthF6, .(Litter_ID_new), numcolwise(mean))

g_means = colMeans(mean_growth_strain[,growth_traits])
m_growth_strain = gather(mean_growth_strain, period, growth, growth_D0D14:growth_D42D56)
m_growth_f6 = gather(growthF6, period, growth, growth_D0D14:growth_D42D56)
m_growth_strain$growth_std = m_growth_strain$growth/g_means[m_growth_strain$period]

m_growth_strain$period = factor(m_growth_strain$period, levels = growth_traits)
stain_growth = ggplot(m_growth_strain, aes(period, growth, group = Strain, color = Strain)) + geom_line(size = 1) + 
  scale_x_discrete(labels = gsub("growth_", "", growth_traits)) + labs(y = "Growth per interval (g)", x = "Intervals") + 
  geom_abline(slope = 0) +
  annotate("segment", x = 1, xend = 3, y= 4, yend = 4, size = 2) + annotate("text", x = 2, y = 4.8, label = "Early \n 0 to 10 days", size = 7) + 
  annotate("segment", x = 6, xend = 9, y= 7, yend = 7, size = 2) + annotate("text", x = 7.5, y = 7.8, label = "Late \n 28 to 56 days", size = 7) +
  theme(legend.position="none")
stain_growth_std = ggplot(m_growth_strain, aes(period, growth_std, group = Strain, color = Strain)) + geom_line(size = 1) + 
  scale_x_discrete(labels = gsub("growth_", "", growth_traits)) + labs(y = "Mean scaled growth (g)", x = "Intervals") + 
  geom_abline(slope = 0, intercept = 1)  + scale_color_discrete(labels = c("A13 - E+", "A22 - E-", "A23 - E-", "A31 - L+", "A41 - L-", "A42 - L-")) + 
  theme(legend.position=c(0.85, 0.85))
#save_plot("~/Desktop/strain_growth-curves.png", plot_grid(stain_growth, stain_growth_std),
#          base_height = 7, base_aspect_ratio = 1, ncol = 2)
library(kinship2)
dev.off()
save_plot("~/Desktop/pedigree.png", drawPedigree(pedigree), base_height = 8)

(f6_growth = ggplot() + geom_line(data = m_growth_f6, 
                                  aes(period, growth, group = ID), 
                                  size = 0.2, alpha = 0.3, color = 'grey'))

(stain_growth = f6_growth + geom_line(data = m_growth_strain, aes(period, growth, group = Strain, color = Strain), size = 1.5) + 
    labs(y = "Growth per interval (g)", x = "") + 
    scale_color_discrete(labels = c("A13 - E+", "A22 - E-", "A23 - E-", "A31 - L+", "A41 - L-", "A42 - L-")) + 
    scale_x_discrete(labels = paste0(c("0 to 14", "14 to 28", "28 to 42", "42 to 56"), "\ndays")) + 
  theme_cowplot() + theme(legend.position=c(0.75, 0.80), text = element_text(size=20), axis.text = element_text(size=20)))


save_plot("~/Dropbox/labbio/posters/2018 - 07 - 19 - Evolution2018/strain_growth-curves_2w.png", stain_growth,
          base_height = 10, base_aspect_ratio = 1)

