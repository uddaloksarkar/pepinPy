import gmpy2 as gp
from gmpy2 import mpq, mpfr
import numpy as np
import math
import os, sys
import random
import time

gp.get_context().precision=20000


def isSAT(dnfclause, sol):         
    for lit in dnfclause:
        if sol[lit-1] == '0':
            return False
    return True
    



def ComputeNumSamples(t, p, thresh, m, delta):

    # # vanilla version
    # N = np.random.binomial(t, p)

    # improved version
    thresh1 = 12 * thresh**2 * m / delta
    thresh2 = (delta / (6 * m))**0.5

    thresh1, thresh2 = mpfr(thresh1), mpfr(thresh2)

    if t * p >= thresh2 :
        if t <= thresh1:
            N = np.random.binomial(int(t), float(p))
        else:
            N = np.random.poisson(float(t * p))
    else:
        N = np.random.binomial(1, float(t*p))
    return N




def getSolutionFromVanillaSampler(dnfClause, nVars):
    sol = ''
    tmpRand = np.random.uniform(0,1,nVars)
    for i in range(1, nVars+1):
        if i in dnfClause:
            sol += '1'
        else:
            if tmpRand[i-1] > 0.5:
                sol += '0'
            else:
                sol += '1'
    return sol


def getSolutionFromSTS(dnfClauseFile, numSolutions):
    kValue = 50
    samplingRounds = int(numSolutions/kValue) + 1
    outputFile = "tmpSTSoutput.out"
    cmd = './samplers/STS -k='+str(kValue)+' -nsamples='+str(samplingRounds)+' '+str(dnfClauseFile)
    cmd += ' > '+str(outputFile)
    os.system(cmd)

    with open(outputFile, 'r') as f:
        lines = f.readlines()

    solList = []
    shouldStart = False
    for j in range(len(lines)):
        if(lines[j].strip() == 'Outputting samples:' or lines[j].strip() == 'start'):
            shouldStart = True
            continue
        if (lines[j].strip().startswith('Log') or lines[j].strip() == 'end'):
            shouldStart = False
        if (shouldStart):
            i = 0
            sol = []
            # valutions are 0 and 1 and in the same order as c ind.
            for x in list(lines[j].strip()):
                if (x == '0'):
                    sol.append(-1*indVarList[i])
                else:
                    sol.append(indVarList[i])
                i += 1
            solList.append(sol)

    solreturnList = solList
    if len(solList) > numSolutions:
        solreturnList = random.sample(solList, numSolutions)
    elif len(solList) < numSolutions:
        print(len(solList))
        print("STS Did not find required number of solutions")
        sys.exit(1)

    os.unlink(outputFile)
    return solreturnList



def getSolutionFromQuickSampler(dnfClauseFile, numSolutions):
    cmd = (
        "./samplers/quicksampler -n "
        + str(numSolutions * 5)
        + " "
        + str(dnfClauseFile)
        )
    print(cmd)
    os.system(cmd)
    cmd = "./samplers/z3 " + str(dnfClauseFile) #+ " > /dev/null 2>&1"
    print(cmd)
    os.system(cmd)
    i = 0
    if numSolutions > 1:
        i = 0

    f = open(dnfClauseFile + ".samples", "r")
    lines = f.readlines()
    f.close()
    f = open(dnfClauseFile + ".samples.valid", "r")
    validLines = f.readlines()
    f.close()
    solList = []
    for j in range(len(lines)):
        # if validLines[j].strip() == "0":
        #     continue
        fields = lines[j].strip().split(":")
        solList.append(fields[1])
        # sol = []
        # i = 0
        # for x in list(fields[1].strip()):
        #     if x == "0":
        #         sol.append(-1*indVarList[i])
        #     else:
        #         sol.append(indVarList[i])
        #     i += 1
        # solList.append(sol)

    solreturnList = solList
    if len(solList) > numSolutions:
        solreturnList = random.sample(solList, numSolutions)
    elif len(solreturnList) < numSolutions:
        print("Did not find required number of solutions")
        exit(1)

    os.unlink(dnfClauseFile+'.samples')
    os.unlink(dnfClauseFile+'.samples.valid')

    return solreturnList


def GenerateSamples(N, dnfClause, delta, m, nVars):
    sampSet = []

    tmpFile = open("tmpClause.cnf", 'w')
    tmpFile.write('p cnf ' + str(nVars) + ' ' + str(len(dnfClause)) + '\n')
    varstr = 'c ind '
    for i in range(1, nVars+1):
        varstr += str(i) + ' '
        if i % 10 == 0 and i < nVars:
            varstr += '0\nc ind '
    tmpFile.write(varstr)
    tmpFile.write('0\n')
    # clauseStr = ''
    for lit in dnfClause:
        tmpFile.write(str(lit) + ' 0\n')
        # clauseStr += str(lit) + ' 0\n'
    # tmpFile.write(clauseStr)
    tmpFile.close()

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
            

    return sampSet





def dnfstream():

    # file handling
    inputFile = "test1.dnf"
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
    m = int(nClause)
    n = int(nVars)
    thresh = max(12 * math.log(24/delta) / eps**2, 6*(math.log(6/delta) + math.log(m)))
    p = 1
    solset = []

    # multi-precision conversion
    thresh, p = mpfr(thresh), mpfr(p)


    for i in range(1, m):

        # print(i)
    
        currClause = lines[i].strip().split()[:-1]
        currClause = list(map(int, currClause))
        clauseWidth = len(currClause)
        t = mpfr(2**(n-clauseWidth))
    
        for s in solset:
            if isSAT(currClause, s):
                solset.remove(s)
    
        while p >= thresh / t:
            for sol in solset:
                if np.random.uniform(0,1) > 0.5 :
                    solset.remove(sol)
            p = p / 2

        N_i = ComputeNumSamples(t, p, thresh, m, delta)


        while N_i + len(solset) > thresh:
            for sol in solset:
                if np.random.uniform(0,1) > 0.5 :
                    solset.remove(sol)
            N_i = N_i / 2 # to change to Binomial (N_i, 1/2)
            p = p / 2

        sol = GenerateSamples(N_i, currClause, delta, m, n)
        solset += sol

    modelCount = int(len(solset)/p)
    
    return modelCount




if __name__ == "__main__":
    
    start_time = time.time()

    modelCount = dnfstream()

    end_time = time.time()

    print("time used by counter (seconds) :", end_time - start_time)
    print("Approx-count : 2^", int(math.log2(modelCount)))