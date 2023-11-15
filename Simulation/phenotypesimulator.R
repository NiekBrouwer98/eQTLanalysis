# devtools::install_github("HannahVMeyer/PhenotypeSimulator")

library(PhenotypeSimulator)

# simulate simple bi-allelic genotypes and estimate kinship
# genotypes <- simulateGenotypes(N = 50 , NrSNP = 1E5,
#                                frequencies = c(0.05,0.1,0.3,0.4),
#                                verbose = FALSE, sampleID = "ID_")

genotypes <- plink_genotypes

# Set parameters
genVar <- 0.6
noiseVar <- 1- genVar
h2s <- 1 #genetic variant effect
phi <- 0.5 #observational noise
rho <- 0.0 #correlation effect
delta <- 0.5 #non-genetic covariant effect
shared <- 0.8
independent <- 1- shared

# kinship estimate based on standardised SNPs
# genotypes_sd <- standardiseGenotypes(genotypes$genotypes)
# kinship <- getKinship(N=50, X=genotypes_sd, verbose = FALSE)

# simulate genectic variant effect (from non-standardised SNP genotypes)
causalSNPs <- getCausalSNPs(N=50, NrCausalSNPs = 1500, genotypes = genotypes$genotypes)
genFixed <- geneticFixedEffects(X_causal = causalSNPs, N=50, P=1500)

# simulate infinitesimal genetic effect
# genBg <- geneticBgEffects(N = 50, kinship = kinship, P = 1500)

#simulate confounder effects
noiseFixed <- noiseFixedEffects(N = 50, P = 1500, NrFixedEffects = 1,
                                NrConfounders = 1, distConfounders = "cat_unif",
                                catConfounders = 2)

# simulate correlated effects with max correlation of 0.8
# correlatedBg <- correlatedBgEffects(N = 10, P = 1500, pcorr = 0.6)

# simulate observational noise effects
noiseBg <- noiseBgEffects(N = 50, P = 1500)

# rescale phenotype components
genFixed_shared_scaled <- rescaleVariance(genFixed$shared, shared * h2s *genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent,
                                               independent * h2s *genVar)
# genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
# genBg_independent_scaled <- rescaleVariance(genBg$independent,
#                                             independent * (1-h2s) * genVar)
noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * phi* noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent,
                                              independent * phi* noiseVar)
# correlatedBg_scaled <- rescaleVariance(correlatedBg$correlatedBg,
                                       # shared * rho * noiseVar)
noiseFixed_shared_scaled <- rescaleVariance(noiseFixed$shared, shared * delta *
                                              noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent,
                                                 independent * delta * noiseVar)

# Total variance proportion have to add up to 1
total <- shared * h2s *genVar + independent * h2s * genVar +
  shared * phi* noiseVar + independent * phi* noiseVar +
  shared * delta * noiseVar + independent * delta * noiseVar

total == 1
#> [1] TRUE

# combine components into final phenotype
# Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component +
#              genBg_independent_scaled$component + noiseBg_independent_scaled$component +
#              genFixed_shared_scaled$component + noiseFixed_shared_scaled$component +
#              genFixed_independent_scaled$component + noiseFixed_independent_scaled$component +
#              correlatedBg_scaled$component)

Y <- scale(noiseFixed_shared_scaled$component +
                       genFixed_shared_scaled$component + 
                       genFixed_independent_scaled$component +
             noiseBg_independent_scaled$component +
             noiseBg_shared_scaled$component)

#Transforming phenotypes
f_custom <- function(x){(x-min(x))/(max(x)-min(x))}

Y_nl <- transformNonlinear(Y, alpha=1, method="custom", f=f_custom)