import sys
import os
import tempfile
from math import log2 

_parallel = False

tempDir = tempfile.gettempdir()
processrankStr = os.environ.get("OMPI_COMM_WORLD_RANK")
if (processrankStr == None):
    processrankStr = 0

processrank = int(processrankStr)
indVarLen = 200
all_file_name = str(tempDir)+'/accurcay_exp_'+str(processrank)+'.log'
cmd = 'find accuracy_exp/5/ -name \*randomDNF_' + str(indVarLen) + '_*.dnf > '+all_file_name
os.system(cmd)

f = open(all_file_name,'r')
allfiles = f.readlines()
f.close()

if _parallel :
    filepos = allfiles[processrank%len(allfiles)].strip()
    print(filepos)
    fileSuffix = filepos.split('/')[-1][:-4]
    if not os.path.isdir("outDir"):
        cmd ='mkdir outDir'
        os.system(cmd)

    outFile = open("allResults" + str(indVarLen), "a")
    toWrite = ""

    # # make equivalent CNF DIMACS 
    # fileposcnf = filepos[:-4] + '.cnf'
    # cmd = 'cp ' + filepos + ' ' + fileposcnf
    # os.system(cmd)

    # with open(fileposcnf, 'r') as file :
    #     filedata = file.read()
    # filedata = filedata.replace('dnf', 'cnf')
    # with open(fileposcnf, 'w') as file:
    #     file.write(filedata)

    # cmd = './ganak ' + fileposcnf + ' > outDir/'+fileSuffix+'.out'
    # os.system(cmd)

    # os.unlink(fileposcnf)

    cmd = './ganak_a ' + filepos + ' > outDir/'+fileSuffix+'.out'
    os.system(cmd)

    f = open("outDir/" + fileSuffix + ".out")
    lines = f.readlines()
    f.close()

    indVarLen = 100
    for line in lines:
        if line.startswith("c Sampling set size:"):
            indVarLen = int(line.strip().split(" ")[-1])
        if line.startswith("s mc"):
            mc = 2**indVarLen - int(line.strip().split(" ")[-1])
            mclog = log2(mc)
            toWrite += str(mc) + "\t" + str(mclog)
            toWrite += "\t"

    cmd = 'python3 dnfstream.py ' + filepos+' > outDir/'+fileSuffix+'.out'
    os.system(cmd)

    f = open("outDir/" + fileSuffix + ".out")
    lines = f.readlines()
    f.close()

    for line in lines:
        if line.startswith("Approx-count : "):
            appmc = int(line.strip().split(":")[1])
            toWrite += str(appmc)
            toWrite += "\t"
        if line.startswith("Approx-count (log) : "):
            toWrite += str(float(line.strip().split("2^")[1]))
            toWrite += "\t"

    if appmc < 1.2 * mc and appmc > 0.8 * mc:
        toWrite += "Success"
    else:
        toWrite += "Failure"
    toWrite += "\t"

    cmd = 'python3 dnfstream2.py ' + filepos+' > outDir/'+fileSuffix+'.out'
    os.system(cmd)

    f = open("outDir/" + fileSuffix + ".out")
    lines = f.readlines()
    f.close()

    for line in lines:
        if line.startswith("Approx-count : "):
            appmc = int(line.strip().split(":")[1])
            toWrite += str(appmc)
            toWrite += "\t"
        if line.startswith("Approx-count (log) : "):
            toWrite += str(float(line.strip().split("2^")[1]))
            toWrite += "\t"

    if appmc < 1.2 * mc and appmc > 0.8 * mc:
        toWrite += "Success"
    else:
        toWrite += "Failure"

    toWrite += "\n"

    outFile.write(toWrite)
    outFile.close()

else:
    for filepos in allfiles:
        filepos = filepos.strip()
        print(filepos)
        fileSuffix = filepos.split('/')[-1][:-4]
        if not os.path.isdir("outDir"):
            cmd ='mkdir outDir'
            os.system(cmd)

        outFile = open("allResults" + str(indVarLen), "a")
        toWrite = filepos + '\t'

        ## make equivalent CNF DIMACS 
        #fileposcnf = filepos[:-4] + '.cnf'
        #cmd = 'cp ' + filepos + ' ' + fileposcnf
        #os.system(cmd)

        #with open(fileposcnf, 'r') as file :
        #    filedata = file.read()
        #filedata = filedata.replace('dnf', 'cnf')
        #with open(fileposcnf, 'w') as file:
        #    file.write(filedata)

        #cmd = './ganak ' + fileposcnf + ' > outDir/'+fileSuffix+'.out'
        #os.system(cmd)

        #os.unlink(fileposcnf)

        cmd = './ganak_a ' + filepos + ' > outDir/' +fileSuffix+ '.out'
        os.system(cmd)

        f = open("outDir/" + fileSuffix + ".out")
        lines = f.readlines()
        f.close()

        flag = False
        for line in lines:
            if line.startswith("c Sampling set size:"):
                indVarLen = int(line.strip().split(" ")[-1])
            if line.startswith("s mc"):
                mc = 2**indVarLen - int(line.strip().split(" ")[-1])
                mclog = log2(mc)
                toWrite += str(mc) + "\t" + str(mclog)
                toWrite += "\t"
                flag = True
        if not flag:
            toWrite += 'X\tX\t'
        
        
        cmd = './pepin -e 0.2 -d 0.1 ' + filepos + ' > outDir/' +fileSuffix+ '.out'
        os.system(cmd)

        f = open("outDir/" + fileSuffix + ".out")
        lines = f.readlines()
        f.close()


        for line in lines:
            if line.startswith("c [dnfs] Low-precision approx num points:"):
                appmc = int(line.strip().split(":")[1])
                toWrite += str(appmc)
                toWrite += "\t"
            if line.startswith("c [dnfs] bucket_size/sampl_prob "):
                toWrite += str(float(line.strip().split("2**")[1]))
                toWrite += "\t"

        if appmc < 1.2 * mc and appmc > 0.8 * mc and flag:
            toWrite += "Success"
        elif flag:
            toWrite += "Failure"

        toWrite += "\t"

        cmd = 'python3 dnfstream.py --samp 3 ' + filepos+' > outDir/'+fileSuffix+'.out'
        os.system(cmd)

        f = open("outDir/" + fileSuffix + ".out")
        lines = f.readlines()
        f.close()

        for line in lines:
            if line.startswith("Approx-count : "):
                appmc = int(line.strip().split(":")[1])
                toWrite += str(appmc)
                toWrite += "\t"
            if line.startswith("Approx-count (log) : "):
                toWrite += str(float(line.strip().split("2^")[1]))
                toWrite += "\t"

        if appmc < 1.2 * mc and appmc > 0.8 * mc and flag:
            toWrite += "Success"
        elif flag:
            toWrite += "Failure"

        toWrite += "\t"

        cmd = 'python3 dnfstream.py --samp 1 ' + filepos+' > outDir/'+fileSuffix+'.out'
        os.system(cmd)

        f = open("outDir/" + fileSuffix + ".out")
        lines = f.readlines()
        f.close()

        flag2 = True
        for line in lines:
            if line.startswith("SAMPLING FAILURE"):
                toWrite += "X\tX\t"
                flag2 = False
                break
            if line.startswith("Approx-count : "):
                appmc = int(line.strip().split(":")[1])
                toWrite += str(appmc)
                toWrite += "\t"
            if line.startswith("Approx-count (log) : "):
                toWrite += str(float(line.strip().split("2^")[1]))
                toWrite += "\t"

        if flag2:
            if appmc < 1.2 * mc and appmc > 0.8 * mc and flag:
                toWrite += "Success"
            elif flag:
                toWrite += "Failure"

        else:
            toWrite += "Failure"

        toWrite += "\n"

        outFile.write(toWrite)
        outFile.close()

