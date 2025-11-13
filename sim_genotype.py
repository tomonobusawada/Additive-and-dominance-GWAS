import numpy as np
import pandas as pd
from statsmodels.distributions.copula.api import GaussianCopula

#1. AR(1) correlation matrix for LD
def ar1_corr(m, rho):
    ind = np.arange(m)
    return rho ** np.abs(ind[:, None] - ind[None, :])

#2. Main simulation function
def sim_genotypes(n, ps, rho= 0.5, seed=None):
    rng = np.random.default_rng(seed)
    p = np.asarray(ps, dtype=float)          
    m = len(p)

    # AR(1) LD structure
    R = ar1_corr(m, rho)

    # Build Gaussian copula
    try:
        gc = GaussianCopula(m)
        gc.corr = R        
    except ValueError:
        gc = GaussianCopula(R)   # statsmodels >= 0.14

    # Draw uniforms with Gaussian dependence 
    n_haps = 2 * n      # Two haplotypes per individual
    try:
        U = gc.random(n_haps, random_state=rng)
    except AttributeError:    
        U = gc.rvs(n_haps, random_state=rng)

    # Convert to haplotype allele counts [0 or 1] and combine 
    hap = (U < p[None, :]).astype(int)
    hap1, hap2 = hap[:n], hap[n:]
    X = hap1 + hap2                                  
    XA = (X - 2*p) / np.sqrt(2*p*(1-p))              # standardized additive matrix
    XD = (-X**2 + (1+2*p)*X - 2*p**2) / (2*p*(1-p))  # orthogonal dominance matrix       

    return XA, XD

n, m = 5000, 500
rng = np.random.default_rng(42)
ps = rng.uniform(0.05, 0.5, size=m)      # allele frequency rangeing from 5% to 50%
XA, XD = sim_genotypes(n, ps, rho= 0.5, seed=1234)

np.savetxt("sim_data/XA_n{}_m{}.tsv".format(n, m), XA)
np.savetxt("sim_data/XD_n{}_m{}.tsv".format(n, m), XD)
