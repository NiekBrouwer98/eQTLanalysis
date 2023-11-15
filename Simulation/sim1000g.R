library("sim1000G")

#Initialisation
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir,"region.vcf.gz")
#or provide own vcf file:
# vcf_file = "C:/Users/niekb/Documents/BEP datasets/hapmap_for_sim1000g/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

vcf = readVCF( vcf_file , maxNumberOfVariants = 2000 , min_maf = 0.1 , max_maf = 1 ) 

  # readGeneticMap( chromosome = 4 )

startSimulation( vcf )

ids = generateUnrelatedIndividuals(10)

genotype_sim100g = retrieveGenotypes(ids)

#Simulating genotypes for 10 unrelated individuals
# # The following function will erase all the indivudals that have been simulated and start over.
# # SIM$reset()
# 
# id = c()
# for(i in 1:10) id[i] = SIM$addUnrelatedIndividual()
# 
# # Show haplotype 1  of first 5 individuals
# #print(SIM$gt1[1:5,1:6])
# 
# # Show haplotype 2
# #print(SIM$gt1[1:5,1:6])
# 
# 
# genotypes = SIM$gt1[1:20,] + SIM$gt2[1:20,]
# 
# print(dim(genotypes))
# 
# str(genotypes)
# 
# library(gplots)
# 
# heatmap.2(cor(genotypes)^2, col=rev(heat.colors(100)) , trace="none",Rowv=F,Colv=F)

#Simulating families
time10families = function() {
  
  
  fam = lapply(1:10, function(x) newFamilyWithOffspring(x,5) )
  fam = do.call(rbind, fam)
  fam
  
}

fam <- time10families() 

# Uncomment the following line to write a plink compatible ped/map file

# writePED(vcf, fam, "C:/Users/niekb/OneDrive/Bachelor Eindproject/Sim_datasets_results/plink_genotype_inputs/test_plink_file")
