# additive and dominance GWAS
This repository includes scripts for joint association testing for additive and dominance genetic effects. To check the the simulated genotype and quantitative phenotype data are also provided. The association testings in [Dominance genetic architecture of complex traits in the Japanese population.](link) are .

## Steps
1. Generation of genotypes
The genotype simulation is based on a Gaussian-copula haplotype model. The LD matirx is an AR(1) correlation with parameter p. Each haplotype allele is obtained by thresholding latent variables. Two independent haplotypes are combined per individual to form diploid genotypes. The resulting additive and domiancne genotypes are mutually uncorrelated.
2. Generation of phenotypes
3. Association test
