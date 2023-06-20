source("./read_F6_phenotypes.R")

library(AtchleyMice)
library(asreml)
library(snpReady)

load("./data/growthGfitF5F6_2w.Rdata")
Gs = array(model_growth$VCV[,grep("animal", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), 4, 4))
G_mcmc = apply(Gs, 2:3, mean)

Rs = array(model_growth$VCV[,grep("units", colnames(model_growth$VCV))], 
           c(nrow(model_growth$VCV), 4, 4))
R_mcmc = apply(Rs, 2:3, mean)

growthF5F6 = mice_growth$F5F6
growth_traits = mice_growth$traits
growthF5F6 = growthF5F6[complete.cases(growthF5F6[,growth_traits]),]
growthF5F6$animal = growthF5F6$ID

F6_ids = growthF5F6$ID[growthF5F6$Gen == "F6"]
F5_ids = unique(c(growthF5F6$Pat_ID[growthF5F6$ID %in% F6_ids], growthF5F6$Mat_ID[growthF5F6$ID %in% F6_ids]))
F5F6_ids = c(F6_ids, F5_ids)
growthF5F6 = growthF5F6[growthF5F6$ID %in% F5F6_ids,]

pedigree = mice_pedigree
pedigree$id = factor(pedigree$id)
pedigree$dam = factor(pedigree$dam)
pedigree$sire = factor(pedigree$sire)
names(pedigree)[1] = "ID"

Ainv = ainverse(pedigree)

n_traits = length(growth_traits)
g_formula = paste0("cbind(", paste(growth_traits, collapse = ", "), ") ~ trait + trait:Sex")

growthF5F6_std = growthF5F6
growthF5F6_std[,growth_traits] = scale(growthF5F6_std[,growth_traits])
growthF5F6_sd = apply(growthF5F6[,growth_traits], 2, sd)

data_growth = as.data.frame(growthF5F6_std)
names(data_growth)

traits = c(1, 2, 3, 4)
current_traits = growth_traits[traits]
n_traits = length(current_traits)
g_formula = paste0("cbind(", paste(current_traits, collapse = ", "), ") ~ trait + trait:Sex")
growth_animal_asr_sv = asreml(as.formula(g_formula),
                           random = ~ us(trait):ped(ID), 
                           rcov = ~ units:us(trait),
                           ginverse = list(ID = Ainv), 
                           data = data_growth,
                           start.values  = TRUE)
sv = growth_animal_asr_sv$gammas.table 

G_mcmc = apply(Gs, 2:3, mean)
G_mcmc = G_mcmc[traits, traits]
sv$Value[1:((n_traits*n_traits-n_traits)/2+n_traits)] = G_mcmc[upper.tri(G_mcmc, diag = TRUE)]
sv$Value[21] = 0.7  
sv$Constraint[21] = "F"
sv$Value[17] = 0.7
sv$Constraint[17] = "F"
data_growth$animal = as.factor(data_growth$animal)
growth_animal_asr = asreml(as.formula(g_formula),
                           random = ~ us(trait):vm(animal, Ainv), 
                           rcov = ~ units:us(trait),
                           data = data_growth)
G = matrix(NA, n_traits, n_traits)
G[upper.tri(G, diag = TRUE)] = growth_animal_asr$G.param$`trait:ped(ID)`$trait$initial
G[lower.tri(G)] = t(G)[lower.tri(G)]
G_init = G
G = G * outer(growthF5F6_sd, growthF5F6_sd)

R = matrix(NA, n_traits, n_traits)
R[upper.tri(R, diag = TRUE)] = growth_animal_asr$R.param$R$trait$initial
R[lower.tri(R)] = t(R)[lower.tri(R)]
R = R * outer(growthF5F6_sd, growthF5F6_sd)

corrG = cov2cor(G)
colnames(corrG) = c("0 to 14\ndays", "14 to 28\ndays", "28 to 42\ndays", "42 to 56\ndays")
diag(corrG) = 0
corrplot.mixed(corrG, upper = "ellipse", mar = c(0, 0, 0, 0), cl.lim = c(-0.8, 0.8), addgrid.col = "black", is.corr = FALSE)

calcG_asreml = function(traits){
    current_traits = growth_traits[traits]
    n_traits = length(current_traits)
    g_formula = paste0("cbind(", paste(current_traits, collapse = ", "), ") ~ trait + trait:Sex")
    growth_animal_asr_sv = asreml(as.formula(g_formula),
                                  random = ~ us(trait):ped(ID), 
                                  rcov = ~ units:us(trait),
                                  ginverse = list(ID = Ainv), 
                                  data = data_growth,
                                  start.values  = TRUE)
    sv = growth_animal_asr_sv$gammas.table 
    G_mcmc = G_init[traits, traits]
    sv$Value[1:((n_traits*n_traits-n_traits)/2+n_traits)] = G_mcmc[upper.tri(G_mcmc, diag = TRUE)]
    growth_animal_asr = asreml(as.formula(g_formula),
                               random = ~ us(trait):ped(ID), 
                               rcov = ~ units:us(trait),
                               ginverse = list(ID = Ainv), 
                               data = data_growth,
                               G.param = sv)
    G = matrix(NA, n_traits, n_traits)
    G[upper.tri(G, diag = TRUE)] = growth_animal_asr$G.param$`trait:ped(ID)`$trait$initial
    G[lower.tri(G)] = t(G)[lower.tri(G)]
    R = matrix(NA, n_traits, n_traits)
    R[upper.tri(R, diag = TRUE)] = growth_animal_asr$R.param$R$trait$initial
    R[lower.tri(R)] = t(R)[lower.tri(R)]
    corrG = cov2cor(G)
    dimnames(G) = dimnames(R) = dimnames(corrG) = list(current_traits, current_traits)
    list(G = G, R = R, corrG = corrG)
}
G_init
calcG_asreml(c(2,3,4))
calcG_asreml(c(1,2,3))
calcG_asreml(c(1,2,4))
calcG_asreml(c(1,2))
calcG_asreml(c(1,3))
calcG_asreml(c(1,4))
calcG_asreml(c(2,3))
calcG_asreml(c(2,4))
calcG_asreml(c(1,2,3,4))

par(mfrow = c(1, 2))
corrplot.mixed(corrG, upper = "ellipse", main = "G")
corrplot.mixed(corrP, upper = "ellipse", main = "P")

P = CalculateMatrix(lm(as.matrix(growthF5F6[,growth_traits]) ~ growthF5F6$Sex))
corrP = cor(residuals(lm(as.matrix(growthF5F6[,growth_traits])~growthF5F6$Sex)))

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
