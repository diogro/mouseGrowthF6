source("./read_F6_phenotypes.R")
library(evolqg)
library(corrplot)

growthF6 = growthF6[complete.cases(growthF6[,growth_traits]),]
growthF6$animal = growthF6$ID

ginverse = inverseA(pedigree)
f6 = rownames(ginverse$Ainv) %in% growthF6$ID
Ainv = ginverse$Ainv[f6, f6]

length(levels(factor(growthF6$ID)))
length(rownames(Ainv))

n_traits = length(growth_traits)
g_formula = paste0("cbind(", paste(growth_traits, collapse = ", "), ") ~ trait + trait:Sex - 1")

growthF6_std = growthF6
growthF6_std[,growth_traits] = scale(growthF6_std[,growth_traits])
growthF6_sd = apply(growthF6[,growth_traits], 2, sd)
prior_bi <- list(G = list(G1 = list(V = diag(n_traits), n = 1.002)),
                 R = list(V = diag(n_traits), n = 1.002))
model_growth <- MCMCglmm(as.formula(g_formula),
                         random = ~us(trait):animal,
                         rcov = ~us(trait):units, family = rep("gaussian", n_traits),
                         pedigree = pedigree, data = as.data.frame(growthF6_std), prior = prior_bi,
                         nitt = 200000 + 1500, thin = 50, burnin = 1500, verbose = TRUE)
summary(model_growth)
save(model_growth, growthF6_sd, file = "./data/growthGfit_2w.Rdata")
load("./data/growthGfit_2w.Rdata")
Gs = array(model_growth$VCV[,grep("animal", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), n_traits, n_traits))
Gs = aaply(Gs, 1, `*`, outer(growthF6_sd, growthF6_sd))
Rs = array(model_growth$VCV[,grep("units", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), n_traits, n_traits))
Rs = aaply(Rs, 1, `*`, outer(growthF6_sd, growthF6_sd))
G = apply(Gs, 2:3, mean)
R = apply(Rs, 2:3, mean)
corrGs = aaply(Gs, 1, cov2cor)
corrG = apply(corrGs, 2:3, mean)
P = CalculateMatrix(lm(as.matrix(growthF6[,growth_traits]) ~ growthF6$Sex))
corrP = cor(residuals(lm(as.matrix(growthF6[,growth_traits])~growthF6$Sex)))

colnames(corrG) = c("0 to 14\ndays", "14 to 28\ndays", "28 to 42\ndays", "42 to 56\ndays")
png("./data/figures/growth_Gmatrix_2w_4t.png", width = 800, height = 800)
par(mfrow=c(1, 1), cex = 2)
corrplot.mixed(corrG, upper = "ellipse")
dev.off()

par(mfrow = c(1, 2))
corrplot.mixed(corrG, upper = "ellipse", main = "G")
corrplot.mixed(corrP, upper = "ellipse", main = "P")

lt = function(x) x[lower.tri(x, diag = TRUE)]
data.frame(P = lt(P), G = lt(G), R = lt(R), T = lt(G + R)) %>%
  ggplot(aes(P, T)) + geom_point() + geom_abline(slope = 1, intercept = 0)

