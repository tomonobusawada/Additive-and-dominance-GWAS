import numpy as np
import pandas as pd

def sim_phenotypes(b, n=5000, M=100, h2add=0.10, h2dom=0.01):

    sigma2e = 1 - h2add - h2dom

    # Generation of random effects
    # Inbreeding depression affects the mean of dominance effects
    beta = np.random.normal(0,h2add**0.5,1)
    delta = np.random.normal(b,h2dom**0.5,1)
    epsilon = np.random.normal(0,sigma2e**0.5, n) 

    betas = np.zeros(M)
    deltas = np.zeros(M)

    causal_index_add = np.random.choice(list(range(M)), 1)
    causal_index_dom = np.random.choice(list(range(M)), 1)
    betas[causal_index_add] = beta
    deltas[causal_index_dom] = delta

    phen = betas@A + deltas@D + epsilon
    np.savetxt('sim_phenotype/phen_ID{}.tsv'.format(b), phen, fmt='%.8e')


np.random.seed(1234)
A = np.loadtxt("sim_genotype/XA_n5000_m100.tsv").T
D = np.loadtxt("sim_genotype/XD_n5000_m100.tsv").T
bs = [1, 3, 5, 7]
bs.sim_phenotypes(b)
