import sys
import os
import tempfile

tempDir = tempfile.gettempdir()
processrankStr = os.environ.get("OMPI_COMM_WORLD_RANK")
if (processrankStr == None):
    processrankStr = 0

processrank = int(processrankStr)
all_file_name = str(tempDir)+'/benchmarks_'+str(processrank)+'.log'
cmd = 'find benchmark -name \*.dnf > '+all_file_name

os.system(cmd)
f = open(all_file_name,'r')
lines = f.readlines()
f.close()
#samplerType = 4-(processrank/len(lines))
#samplerType = 1
#if (samplerType < 1):
#    exit(-1)
filepos = lines[processrank%len(lines)].strip()
fileSuffix = filepos.split('/')[-1][:-4]
if not os.path.isdir("outDir"):
    cmd ='mkdir outDir'
    os.system(cmd)
#cmd = 'ulimit -v 4000000; python Verifier.py --sampler '+str(samplerType)+' --inverted 1 --reverse 0 --seed 0 --exp 1000 '+filepos+'  outDir/sampler_'+str(samplerType)+'_'+fileSuffix+'.out'
#print(cmd)
cmd = 'python3 dnfstream.py ' + filepos+' > '+fileSuffix+'.out'
os.system(cmd)
#print(cmd)
