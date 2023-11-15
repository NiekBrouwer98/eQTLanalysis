devtools::install_github("HannahVMeyer/PhenotypeSimulator")

library(PhenotypeSimulator)

# read genotypes from external file
# use one of the sample genotype file provided in the
# extdata/genotypes/subfolders (e.g.extdata/genotypes/hapgen )
genotypefile <- readStandardGenotypes(N=100, "file",
                                      format= 'delim', delimiter = '\t')
  
  
# simulate phenotype with the same five phenotype components
# and settings as above; display progress via verbose=TRUE
phenotype <- runSimulation(N = 100, P = 15, genotypefile = genotypefile,
                           format = "delim", genoDelimiter = '\t', cNrSNP = 30, genVar = genVar, h2s = h2s,
                           phi = 0.6, delta = 0.3, distBetaGenetic = "unif", mBetaGenetic = 0.5,
                           sdBetaGenetic = 1, NrFixedEffects = 4, NrConfounders = c(1,
                                                                                    2, 1, 2), pIndependentConfounders = c(0, 1, 1, 0.5),
                           distConfounders = c("bin", "cat_norm", "cat_unif", "norm"),
                           probConfounders = 0.2, catConfounders = c(3, 4), pcorr = 0.8,
                           verbose = TRUE)

#save
out <- savePheno(phenotype, directory="/tmp",
                 outstring="test_simulation",
                 format=c("csv", "plink"), verbose=FALSE)