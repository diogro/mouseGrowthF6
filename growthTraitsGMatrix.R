source("./read_F6_phenotypes.R")

drawPedigree(pedigree)

growthF6$animal = growthF6$ID
ginverse = inverseA(pedigree)

which(pedigree$dam %in% pedigree$sire)

which(!levels(factor(growthF6$ID)) %in% rownames(ginverse$Ainv))

f6 = rownames(ginverse$Ainv) %in% growthF6$ID
Ainv = ginverse$Ainv[f6, f6]

length(levels(factor(growthF6$ID)))
length(rownames(Ainv))

n_traits = length(growth_traits)
g_formula = paste0("cbind(", paste(growth_traits, collapse = ", "), ") ~ trait + trait:Sex")

prior_bi <- list(G = list(G1 = list(V = diag(n_traits), n = 1.002)),
                 R = list(V = diag(n_traits), n = 1.002))
#model_bi <- MCMCglmm(as.formula(g_formula),
#                     random = ~us(trait):animal,
#                     rcov = ~us(trait):units, family = rep("gaussian", n_traits),
#                     pedigree = pedigree, data = as.data.frame(growthF6), prior = prior_bi,
#                     nitt = 650000, thin = 50, burnin = 15000, verbose = TRUE)
summary(model_bi)
#save(model_bi, file = "./data/GrowthGfit.Rdata")
load("./data/GrowthGfit.Rdata")
Gs = array(model_bi$VCV[,grep("animal", colnames(model_bi$VCV))], c(nrow(model_bi$VCV), n_traits, n_traits))
G = apply(Gs, 2:3, mean)
corrGs = aaply(Gs, 1, cov2cor)
corrG = apply(corrGs, 2:3, mean)
corrP = cor(residuals(lm(as.matrix(growthF6[,growth_traits])~growthF6$Sex)))
corrplot.mixed(corrG, upper = "ellipse")
library(evolqg)
MatrixCompare(gcorr, corrG)
par(mfrow = c(1, 2))
corrplot.mixed(corrG, upper = "ellipse", main = "G")
corrplot.mixed(corrP, upper = "ellipse", main = "P")
