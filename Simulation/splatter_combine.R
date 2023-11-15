search()

library(splatter)
library(scater)
library(fitdistrplus)

# estimate parameters from dataset
# counts <- t(apply(Y, 1, function(x)((x-min(x))/(max(x)-min(x))*100)))

# colmean <- colMeans(counts)
# rowmean <- rowMeans(counts)

for (i in 1:10){
  name <- paste("counts_ID", i, sep = "")
  assign(name, Y_nl[i,])
  
  shape <- paste("shape", i, sep = "")
  rate <- paste("rate", i, sep = "")
  
  fit.gamma <- fitdist(Y_nl[i,], distr = "gamma", method = 'mle')
  summary(fit.gamma)
  assign(shape, fit.gamma$estimate['shape'])
  assign(rate, fit.gamma$estimate['rate'])
  # plot(fit.gamma, demp=TRUE)
  
}

#Setting the parameters
params <- newSplatParams()
params <- setParams(params, nGenes = 1500, batchCells = c(4000),
                    group.prob = c(0.8,0.2),
                    de.downProb = c(0.5, 0.1))
                    # de.facLoc = c(0.01,0.5,0.01,0.5),
                    # de.facScale = c(0.01,0.5,0.01,0.5))

# params1 <- setParams(params, update = list(mean.shape=shape1, mean.rate=rate1))
# params2 <- setParams(params, update = list(mean.shape=shape2, mean.rate=rate2))

#Simulating the experiments
sim <- splatSimulate(params, method = "groups")
sim_norm <- normalizeSCE(sim)
sim_count_table <- counts(sim)

#PCA plot:
plotPCA(sim_norm, colour_by = "Group")

#t-SNE plot:
plotTSNE(sim_norm, colour_by = "Group")
