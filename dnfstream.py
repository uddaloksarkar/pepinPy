import gmpy2 as gp
import numpy as np
import math


def isSAT(dnfclause, sol):         # to fix clause and sol
    val = True
    # for lit in sol:
    #     if lit > 0:
            
    #     if lit in dnfclause:
    return val


def ComputeNumSamples(t, p):

    return


def GenerateSamples(N, Formula):
    sampSet = []

    k = 0
    for j in range(1, N)

    return sampSet


def constructNewFile(currClause, )

# file handling
inputFile = "test.dnf"
f = open(inputFile, "r")
lines = f.readlines()
f.close()

initLine = lines[0].strip().split()


if initLine[0] == "p":
    nVars = initLine[2]
    nClause = initLine[3]

print(nVars, nClause)


# parameters / initialization
eps = 0.8
delta = 0.36
m = nClause
thresh = max(12 * math.log(24/delta) / eps^2, 6*(math.log(6/delta) + math.log(m)))
p = 1
solset = []



for i in range(1, nClause):
    
    currClause = lines[i].strip().split()
    clauseWidth = len(currClause)
    t = 2**(nVars-clauseWidth)
    
    for s in solset:
        if isSAT(currClause, s):
            solset = solset.remove(s)
    
    while p >= thresh / t:
        for sol in solset:
            if np.random.uniform(0,1) > 0.5 :
                solset.remove(sol)
        p = p / 2

    N_i = ComputeNumSamples(t, p)


    while N_i + len(solset) > thresh:
        for sol in solset:
            if np.random.uniform(0,1) > 0.5 :
                solset.remove(sol)
        N_i = N_i / 2 # to change to Binomial (N_i, 1/2)
        p = p / 2

    sol = GenerateSamples(N_i, currClause, delta, m)
    solset.append(sol)

print(len(solset)/p)
