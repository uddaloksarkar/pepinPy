import gmpy2 as gp
from gmpy2 import mpz, mpq, mpfr
import numpy as np
import math
import os, sys
import random
import time
import argparse
import cProfile
from mprandom import mpbinomial

gp.get_context().precision=20000


def isSAT(Box, sol):
    
    b = True
    for dim in range(len(sol)):
        interval = sol[dim]
        ptDim = np.random.uniform(interval[0], interval[1])
        if ptDim > Box[dim][1] or ptDim < Box[dim][0]:
            b *= False

    return b
    


def ComputeNumSamples(t, p, thresh, m, delta, method):

    print(f"sample n : {t}, p : {p}")

    if method == 1:
        # vanilla version
        try:
            N = np.random.binomial(t, p)
        except OverflowError:
            print("SAMPLING FAILURE!")
            exit("SAMPLING FAILURE!")

    elif method == 2:
        # improved version
        thresh1 = 12 * thresh**2 * m / delta
        thresh2 = (delta / (6 * m))**0.5

        thresh1, thresh2 = mpfr(thresh1), mpfr(thresh2)

        if t * p >= thresh2 :
            if t <= thresh1:
                print("binomial")
                N = np.random.binomial(int(t), float(p))
            else:
                print("poisson")
                N = np.random.poisson(float(t * p))
        else:
            print("small binomial")
            N = np.random.binomial(1, float(t*p))

    elif method == 3:
        # mp version
        N = mpbinomial(int(t), p, err=delta / (6 * m))

    print(f"ni : {N}")
    
    return N



def getSolutionFromVanillaSampler(dnfClause, nVars):
    sol = []
    tmpRand = np.random.uniform(0,1,nVars)
    for i in range(1, nVars+1):
        if i in dnfClause:
            sol.append(i)
        else:
            if tmpRand[i-1] > 0.5:
                sol.append(-i)
            else:
                sol.append(i)
    return sol



def constructLazySample(dnfClause):
    # sol = []
    # for lit in dnfClause:
    #     sol.append(lit)
    return dnfClause


def GenerateSamples(N, dnfClause, delta, m, nVars, thresh):
    sampSet = []

    # tmpFile = open("tmpClause.cnf", 'w')
    # tmpFile.write('p cnf ' + str(nVars) + ' ' + str(len(dnfClause)) + '\n')
    # varstr = 'c ind '
    # for i in range(1, nVars+1):
    #     varstr += str(i) + ' '
    #     if i % 10 == 0 and i < nVars:
    #         varstr += '0\nc ind '
    # tmpFile.write(varstr)
    # tmpFile.write('0\n')
    # # clauseStr = ''
    # for lit in dnfClause:
    #     tmpFile.write(str(lit) + ' 0\n')
    #     # clauseStr += str(lit) + ' 0\n'
    # # tmpFile.write(clauseStr)
    # tmpFile.close()

    if False: #nVars - len(dnfClause) - 2 * math.log2(1+thresh) <= math.log2(6*m/delta):

        if N > 0 :
            k = 0
            lmt = int(N * (math.log(N) + math.log(6/delta) + math.log(m)))
            for j in range(int(lmt)):
                s = getSolutionFromVanillaSampler(dnfClause, nVars)
                if s not in sampSet:
                    sampSet.append(s)
                    k += 1
                if k == N: break

            # s = getSolutionFromSTS("tmpClause.cnf", lmt)
            
            # s = getSolutionFromQuickSampler("tmpClause.cnf", lmt)
            
    else:
        
        for j in range(N):
            sampSet.append(constructLazySample(dnfClause))

    return sampSet


def parseInput(inputFile):
    f = open(inputFile, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.startswith('c dim'):
            ndim = line.strip().split()[-1]
        if line.startswith('c nbox'):
            nbox = line.strip().split()[-1]
    return ndim, nbox


def calculateVolume(box):
    volume = mpfr('1')
    for interval in box:
        volume = gp.mul(volume, (interval[1] - interval[0]))
    return volume




def dnfstream():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--eps", type=float, help="default = 0.2", default=0.2, dest="eps"
    )
    parser.add_argument(
        "--delta", type=float, help="default = 0.1", default=0.1, dest="delta"
    )
    parser.add_argument("--seed", type=int, dest="seed", default=10)
    parser.add_argument("--samp", type=int, dest="samp", default=1)
    parser.add_argument("input", help="input file")

    args = parser.parse_args()

    # file handling
    inputFile = args.input #"test4.dnf" # args.input
    f = open(inputFile, "r")
    lines = f.readlines()
    f.close()

    seed = args.seed
    sampMethod = args.samp

    np.random.seed(seed)


    ndim, nbox = parseInput(inputFile)

    print(ndim, nbox)


    # parameters / initialization
    eps = args.eps
    delta = args.delta
    m = int(nbox)
    ndim = int(ndim)
    thresh = max(12 * math.log(24/delta) / eps**2, 6*(math.log(6/delta) + math.log(m)))
    # thresh = 4*math.log2(m+1)/(eps**2)*math.log2(1.0/delta) 
    p = 1
    solset = []

    # multi-precision conversion
    thresh, p = mpfr(thresh), mpfr(p)

    line = 0  # line 0 corresponds to p dnf
    cl = 0
    # for i in range(1, m+1):
    while True:
        if lines[line].startswith("c") or lines[line].startswith("p") or lines[line].startswith("w"):
            line += 1
            continue 
        currBox = lines[line].strip().split()
        line += 1
        currBox = list(map(float, currBox))

        tmpBox = []
        count = 0

        assert(len(currBox)%2 == 0)

        while count < len(currBox):
            tmpBox.append([currBox[count], currBox[count + 1]])
            count += 2
        currBox = tmpBox
        print(f"adding box {currBox}")

        t = calculateVolume(currBox)

        t = mpz(t)  # for binomial sampling
    
        for s in solset:
            if isSAT(currBox, s):
                solset.remove(s)
    
        if cl == 1 and p >= thresh / t:
            # pow = gp.ceil(gp.log2(p * t / thresh))
            pow = gp.ceil(gp.log2(gp.div(gp.mul(p, t), thresh)))
            # p = p / 2**pow
            p = gp.div(p, 2**pow)

        while p * t >= thresh:
            for sol in solset:
                if np.random.uniform(0,1) > p : # this was 0.5 before
                    solset.remove(sol)
            # p = p / 2
            p = gp.div(p,2)
        print(f"p: {p} | thresh : {int(thresh)} | bucket : {len(solset)}")

        N_i = ComputeNumSamples(t, p, thresh, m, delta, sampMethod)

        Npast = N_i
        while N_i + len(solset) > thresh:
            for sol in solset:
                if np.random.uniform(0,1) > p : # this was 0.5 before
                    solset.remove(sol)
            N_i = np.random.binomial(N_i , 1/2)
            p = p / 2
            print(f"bucket reduced to : {len(solset)}")

        print(f"old ni : {Npast}, new ni: {N_i}")

        sol = GenerateSamples(N_i, currBox, delta, m, ndim, thresh) # ndim just a placeholder :: not in work
        solset += sol
        cl += 1
        if cl == m : break

        seed += 1
        
    print(1/p)

    modelCount = int(len(solset)/p)
    
    return modelCount




if __name__ == "__main__":
    
    start_time = time.time()

    modelCount = dnfstream()

    end_time = time.time()

    print("time used by counter (seconds) :", end_time - start_time)
    print("Approx-count : ", int(modelCount))
    print("Approx-count (log) : 2^", (math.log2(int(modelCount))))
