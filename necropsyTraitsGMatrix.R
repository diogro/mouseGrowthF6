source("./read_F6_phenotypes.R")

necropsyF6 = necropsyF6[complete.cases(necropsyF6[,necropsy_traits]),]
necropsyF6$animal = necropsyF6$ID

f6 = rownames(ginverse$Ainv) %in% necropsyF6$ID
Ainv = ginverse$Ainv[f6, f6]

length(levels(factor(necropsyF6$ID)))
length(rownames(Ainv))

n_traits = length(necropsy_traits)
g_formula = paste0("cbind(", paste(necropsy_traits, collapse = ", "), ") ~ trait + trait:Sex")

prior_bi <- list(G = list(G1 = list(V = diag(n_traits), n = 1.002)),
                 R = list(V = diag(n_traits), n = 1.002))
#model_necropsy <- MCMCglmm(as.formula(g_formula),
#                     random = ~us(trait):animal,
#                     rcov = ~us(trait):units, family = rep("gaussian", n_traits),
#                     pedigree = pedigree, data = as.data.frame(necropsyF6), prior = prior_bi,
#                     nitt = 65000, thin = 50, burnin = 1500, verbose = TRUE)
#summary(model_necropsy)
#save(model_necropsy, file = "./data/necropsyGfit.Rdata")
load("./data/necropsyGfit.Rdata")
Gs = array(model_necropsy$VCV[,grep("animal", colnames(model_necropsy$VCV))], 
           c(nrow(model_necropsy$VCV), n_traits, n_traits))
G = apply(Gs, 2:3, mean)
corrGs = aaply(Gs, 1, cov2cor)
corrG = apply(corrGs, 2:3, mean)
corrP = cor(residuals(lm(as.matrix(necropsyF6[,necropsy_traits])~necropsyF6$Sex)))
corrplot.mixed(corrG, upper = "ellipse")
library(evolqg)
MatrixCompare(gcorr, corrG)
par(mfrow = c(1, 2))
necropsy_traits
order = c(1, 3, 5, 2, 4)
corrplot.mixed(corrG[order, order], upper = "ellipse", main = "G")
corrplot.mixed(corrP[order, order], upper = "ellipse", main = "P")
