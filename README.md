# additive and dominance GWAS
This repository provides scripts for joint association testing for additive and dominance genetic effects. The provided simulated genotype data will be used as an input.

## Steps
1. Generation of genotypes
The genotype simulation is based on a Gaussian-copula haplotype model. The LD matirx is an AR(1) correlation with parameter p. Each haplotype allele is obtained by thresholding latent variables. Two independent
haplotypes are combined per individual to form diploid genotypes.
The resulting additive and domiancne genotypes are mutually uncorrelated.
3. Association test
