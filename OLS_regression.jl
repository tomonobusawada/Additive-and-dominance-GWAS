using PyCall, LinearAlgebra
using CSV, DataFrames, Distributed, Printf
using CodecZlib, Mmap,  DelimitedFiles, TranscodingStreams
using Folds

py"""
import numpy as np
A = np.loadtxt('sim_genotype/XA_n5000_m500.tsv').T
D = np.loadtxt('sim_genotype/XD_n5000_m500.tsv').T
"""

A = PyArray(py"A"o)
D = PyArray(py"D"o)
m = size(D, 1)
n = size(D, 2)

trait = "simulation"
ID = parse(Int, ENV["b"])
y = readdlm(@sprintf("sim_phenotype/phen_ID%d.tsv", ID))

function joint_regression(j)
    Q = hcat(A[j,:], D[j, :])
    beta = (transpose(Q)*Q) \ (transpose(Q)*y)

    res = y - Q*beta
    sigma2_hat = dot(res, res) / (n - 2)
    var = sigma2_hat * ((transpose(Q)*Q) \ I)

    return beta[1],sqrt(var[1,1]),(beta[1])^2/var[1,1],beta[2],sqrt(var[2,2]),(beta[2])^2/var[2,2] 
end 

df = DataFrame(Folds.map(joint_regression, 1:m))
df = DataFrame(df, [:Beta, :se_Beta, :Chisq_ADD, :Delta, :se_Delta, :Chisq_DOM])

CSV.write(@sprintf("sumstat/sumstat.%s.ID%d.tsv", trait, ID), df, delim="\t")
