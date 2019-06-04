source("./read_F6_phenotypes.R")
library(evolqg)
library(corrplot)

growthF5F6 = growthF5F6[complete.cases(growthF5F6[,growth_traits]),]
growthF5F6$animal = growthF5F6$ID

ginverse = inverseA(pedigree)
f6 = rownames(ginverse$Ainv) %in% growthF5F6$ID
Ainv = ginverse$Ainv[f6, f6]

length(levels(factor(growthF5F6$ID)))
length(rownames(Ainv))

levels(factor(growthF5F6$ID))[!levels(factor(growthF5F6$ID)) %in% rownames(Ainv)]

n_traits = length(growth_traits)
g_formula = paste0("cbind(", paste(growth_traits, collapse = ", "), ") ~ trait + trait:Sex - 1")

growthF5F6_std = growthF5F6
growthF5F6_std[,growth_traits] = scale(growthF5F6_std[,growth_traits])
growthF5F6_sd = apply(growthF5F6[,growth_traits], 2, sd)
prior_bi <- list(G = list(G1 = list(V = diag(n_traits), n = 1.002)),
                 R = list(V = diag(n_traits), n = 1.002))
model_growth <- MCMCglmm(as.formula(g_formula),
                         random = ~us(trait):animal,
                         rcov = ~us(trait):units, family = rep("gaussian", n_traits),
                         pedigree = pedigree, data = as.data.frame(growthF5F6_std), prior = prior_bi,
                         nitt = 200000 + 1500, thin = 50, burnin = 1500, verbose = TRUE)
summary(model_growth)
save(model_growth, growthF5F6_sd, file = "./data/growthGfitF5F6_2w.Rdata")
load("./data/growthGfit_2w.Rdata")
Gs = array(model_growth$VCV[,grep("animal", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), n_traits, n_traits))
Gs = aaply(Gs, 1, `*`, outer(growthF5F6_sd, growthF5F6_sd))
Rs = array(model_growth$VCV[,grep("units", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), n_traits, n_traits))
Rs = aaply(Rs, 1, `*`, outer(growthF5F6_sd, growthF5F6_sd))
G = apply(Gs, 2:3, mean)
R = apply(Rs, 2:3, mean)
corrGs = aaply(Gs, 1, cov2cor)
corrG = apply(corrGs, 2:3, mean)
P = CalculateMatrix(lm(as.matrix(growthF5F6[,growth_traits]) ~ growthF5F6$Sex))
corrP = cor(residuals(lm(as.matrix(growthF5F6[,growth_traits])~growthF5F6$Sex)))

colnames(corrG) = c("0 to 14\ndays", "14 to 28\ndays", "28 to 42\ndays", "42 to 56\ndays")
diag(corrG) = 0
png("/home/diogro/Dropbox/labbio/posters/2018 - 07 - 19 - Evolution2018/growth_Gmatrix_2w_4t.png", width = 2200, height = 2200)
par(mfrow=c(1, 1), cex = 4)
corrplot.mixed(corrG, upper = "ellipse", mar = c(0, 0, 0, 0), cl.lim = c(-0.8, 0.8), addgrid.col = "black", is.corr = FALSE)
dev.off()

par(mfrow = c(1, 2))
corrplot.mixed(corrG, upper = "ellipse", main = "G")
corrplot.mixed(corrP, upper = "ellipse", main = "P")

lt = function(x) x[lower.tri(x, diag = TRUE)]
data.frame(P = lt(P), G = lt(G), R = lt(R), T = lt(G + R)) %>%
  ggplot(aes(P, T)) + geom_point() + geom_abline(slope = 1, intercept = 0)

corridor = as.matrix(read_delim('~/Dropbox/labbio/data/evomod_julia_data/corridor.txt',delim = "\t", col_names = FALSE))
divergent = as.matrix(read_delim('~/Dropbox/labbio/data/evomod_julia_data/divergent.txt',delim = "\t", col_names = FALSE))
stabilizing = as.matrix(read_delim('~/Dropbox/labbio/data/evomod_julia_data/stabilizing.txt',delim = "\t", col_names = FALSE))
colnames(corridor) = colnames(divergent) =  colnames(stabilizing) = paste0(1:4)

png("/home/diogro/Dropbox/labbio/posters/2018 - 07 - 19 - Evolution2018/corridor.png", width = 800, height = 800)
par(mfrow = c(1, 1), cex = 3)
corrplot.mixed(corridor, upper = "ellipse")
dev.off()
png("/home/diogro/Dropbox/labbio/posters/2018 - 07 - 19 - Evolution2018/divergent.png", width = 800, height = 800)
par(mfrow = c(1, 1), cex = 3)
corrplot.mixed(divergent, upper = "ellipse")
dev.off()
png("/home/diogro/Dropbox/labbio/posters/2018 - 07 - 19 - Evolution2018/stabilizing.png", width = 800, height = 800)
par(mfrow = c(1, 1), cex = 3)
corrplot.mixed(stabilizing, upper = "ellipse")
dev.off()
