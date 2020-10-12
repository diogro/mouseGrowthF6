source("./read_F6_phenotypes.R")
growthF5F6 = mice_growth$F5F6

library(stanAnimal)

growthF5F6 = growthF5F6[complete.cases(growthF5F6[,mice_growth$traits]),]
growthF5F6$animal = growthF5F6$ID

F6_ids = growthF5F6$ID[growthF5F6$Gen == "F6"]
F5_ids = unique(c(growthF5F6$Pat_ID[growthF5F6$ID %in% F6_ids], growthF5F6$Mat_ID[growthF5F6$ID %in% F6_ids]))
F5F6_ids = c(F6_ids, F5_ids)
A = makeA(mice_pedigree)

length(levels(factor(growthF5F6$ID)))
length(rownames(A))

levels(factor(growthF5F6$ID))[!levels(factor(growthF5F6$ID)) %in% rownames(A)]

g_formula = paste0("cbind(", paste(mice_growth$traits, collapse = ", "), ") ~ Sex")
model_input = genAnimalModelInput(g_formula, growthF5F6, A)
stan_animal_model  = lmm_animal(model_input$Y, model_input$X, model_input$A, chains = 4, cores = 4)
#write_rds(list(fit = stan_animal_model, data = model_input), file = "./Rdatas/growth_9trait_stanAnimalModel.rds")
summary(stan_animal_model, pars = "corrG")

colMeans(stan_animal_model@G)
growth_animal = rstan::extract(stan_animal_model)
colMeans(growth_animal$G)
corrplot.mixed(colMeans(growth_animal$corrG), upper = "ellipse")


n_traits = length(mice_growth$traits)
g_formula = paste0("cbind(", paste(mice_growth$traits, collapse = ", "), ") ~ trait + trait:Sex - 1")

growthF5F6_std = growthF5F6
growthF5F6_std[,mice_growth$traits] = scale(growthF5F6_std[,mice_growth$traits])
growthF5F6_sd = apply(growthF5F6[,mice_growth$traits], 2, sd)
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

par(mfrow = c(1, 2))
corrplot.mixed(colMeans(growth_animal$corrG), upper = "ellipse")
corrplot.mixed(cov2cor(G), upper = "ellipse")


P = CalculateMatrix(lm(as.matrix(growthF5F6[,mice_growth$traits]) ~ growthF5F6$Sex))
corrP = cor(residuals(lm(as.matrix(growthF5F6[,mice_growth$traits]) ~ growthF5F6$Sex)))

colnames(corrG) = c("0 to 14\ndays", "14 to 28\ndays", "28 to 42\ndays", "42 to 56\ndays")
diag(corrG) = 0
png("/home/diogro/Dropbox/labbio/articles/TeseDoutorado/chapter_atchley/media/growth_Gmatrix_2w_4t.png", width = 2200, height = 2200)
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
