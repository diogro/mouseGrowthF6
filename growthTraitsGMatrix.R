library(MCMCglmm)
library(pedantics)

source("./read_F6_phenotypes.R")

#drawPedigree(pedigree)

growthF6$animal = growthF6$ID
ginverse = inverseA(pedigree)

which(!growthF6$animal %in% rownames(ginverse$Ainv))

prior_bi <- list(G = list(G1 = list(V = diag(2), n = 1.002)),
                 R = list(V = diag(2), n = 1.002))
model_bi <- MCMCglmm(cbind(growth_D3D7, Final_weight) ~ trait + trait:Sex,
                     random = ~us(trait):animal,
                     rcov = ~us(trait):units, family = c("gaussian", "gaussian"),
                     pedigree = pedigree, data = growthF6, prior = prior_bi,
                     nitt = 65000, thin = 50, burnin = 15000, verbose = TRUE)
