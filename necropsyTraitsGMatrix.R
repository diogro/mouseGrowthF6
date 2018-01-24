source("./read_F6_phenotypes.R")

necropsyF6 = necropsyF6[complete.cases(necropsyF6[,necropsy_traits]),]
necropsyF6$animal = necropsyF6$ID

drawPedigree(pedigree)

f6 = rownames(ginverse$Ainv) %in% necropsyF6$ID
Ainv = ginverse$Ainv[f6, f6]

length(levels(factor(necropsyF6$ID)))
length(rownames(Ainv))

n_traits = length(necropsy_traits)
g_formula = paste0("cbind(", paste(necropsy_traits, collapse = ", "), ") ~ trait + trait:Sex - 1")

necropsyF6_std = necropsyF6
necropsyF6_std[,necropsy_traits] = scale(necropsyF6_std[,necropsy_traits])
necropsyF6_sd = apply(necropsyF6[,necropsy_traits], 2, sd)
prior_bi <- list(G = list(G1 = list(V = diag(n_traits), n = 1.002)),
                 R = list(V = diag(n_traits), n = 1.002))
model_necropsy <- MCMCglmm(as.formula(g_formula),
                     random = ~us(trait):animal,
                     rcov = ~us(trait):units, family = rep("gaussian", n_traits),
                     pedigree = pedigree, data = as.data.frame(necropsyF6_std), prior = prior_bi,
                     nitt = 200000 + 1500, thin = 50, burnin = 1500, verbose = TRUE)
summary(model_necropsy)
save(model_necropsy, necropsyF6_sd, file = "./data/necropsyGfit.Rdata")
load("./data/necropsyGfit.Rdata")
Gs = array(model_necropsy$VCV[,grep("animal", colnames(model_necropsy$VCV))], 
           c(nrow(model_necropsy$VCV), n_traits, n_traits))
Gs = aaply(Gs, 1, `*`, outer(necropsyF6_sd, necropsyF6_sd))
Rs = array(model_necropsy$VCV[,grep("units", colnames(model_necropsy$VCV))], 
           c(nrow(model_necropsy$VCV), n_traits, n_traits))
Rs = aaply(Rs, 1, `*`, outer(necropsyF6_sd, necropsyF6_sd))
G = apply(Gs, 2:3, mean)
R = apply(Rs, 2:3, mean)
corrGs = aaply(Gs, 1, cov2cor)
corrG = apply(corrGs, 2:3, mean)
P = CalculateMatrix(lm(as.matrix(necropsyF6[,necropsy_traits]) ~ necropsyF6$Sex))
corrP = cor(residuals(lm(as.matrix(necropsyF6[,necropsy_traits])~necropsyF6$Sex)))
corrplot.mixed(corrG, upper = "ellipse")
library(evolqg)
par(mfrow = c(1, 2))
corrplot.mixed(corrG[order, order], upper = "ellipse", main = "G")
corrplot.mixed(corrP[order, order], upper = "ellipse", main = "P")

lt = function(x) x[lower.tri(x, diag = TRUE)]
data.frame(P = lt(P), G = lt(G), R = lt(R), T = lt(G + R)) %>%
  ggplot(aes(P, T)) + geom_point() + geom_abline(slope = 1, intercept = 0)
