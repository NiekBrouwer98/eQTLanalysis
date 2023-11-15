# Load package
library(splatter)

# Create mock data
library(scater)
set.seed(1)
sce <- mockSCE()

params <- splatEstimate(sce)

counts(sce)[1:5, 1:5]
#information about genes
head(rowData(sce))
#information about cells
head(colData(sce))

sim <- splatSimulate(params)
