# additive and dominance GWAS
This repository includes scripts for performing joint association analyses that test both additive and dominance genetic effects, with example simulated genotype and phenotype data. This association tests are used in [Dominance genetic architecture of complex traits in the Japanese population.](link).

## Steps
1. Generation of genotypes\
   The genotype simulation is based on a Gaussian-copula haplotype model. The LD matirx is an AR(1) correlation with parameter p. Each haplotype allele is obtained by thresholding latent variables. Two independent haplotypes are combined per individual to form diploid genotypes. The resulting additive and domiancne genotypes are standardized and mutually uncorrelated.
2. Generation of phenotypes assuming varying degree of inbreeding\
   
3. Association test\
   
