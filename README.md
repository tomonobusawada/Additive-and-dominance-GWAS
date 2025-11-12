## additive and dominance GWAS
This repository includes scripts for performing joint association analyses that test both additive and dominance genetic effects, with example simulated genotype and phenotype data. This association tests are used in [Dominance genetic architecture of complex traits in the Japanese population.](link).

### Steps
1. Generation of genotypes\
   The genotype simulation is based on a Gaussian-copula haplotype model which generates haplotype from correlated latent normal variables whose covariance matrix defines LD structure (AR(1) correlation is used here). Two independent haplotypes are combined per individual to form diploid genotypes. The resulting additive and domiancne genotypes are standardized and mutually uncorrelated.
3. Generation of phenotypes assuming varying degree of inbreeding\
   In this simulation, additive and dominance causal variants are selected independently, one for each effect type. The additive heritability and domiannce heritability are set to 15% and 1%, respectively. All other simulation settings follow those described in [the reference](link).  
4. Association test\
   The simulated phenotypes are jointly regressed on the additive and dominance genotypes. We computed marginal association statistics using OLS: $\frac{1}{n}X^{T}_{A, j}\mathbf{y}$ and $\frac{1}{n}X^{T}_{D, j}\mathbf{y}$ for each SNP $j$.
   
