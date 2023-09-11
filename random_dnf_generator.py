#!/usr/bin/python2

import sys
import os
import random
import math

random.seed(1)
nLow = int(sys.argv[1])
nHigh =  int(sys.argv[2])
nStep = int(sys.argv[3])
mDensityLow = float(sys.argv[4])
mDensityHigh = float(sys.argv[5])
mDensityStep = float(sys.argv[6])
mSizeType = int(sys.argv[7])
mSizeLow = float(sys.argv[8])
mSizeHigh = float(sys.argv[9])
mSizeStep = float(sys.argv[10])

instances = int(sys.argv[11])

outputDir = sys.argv[12]
monotone = int(sys.argv[13])

for n in range(nLow,nHigh,nStep):
	m = mDensityLow
	while m < mDensityHigh:
		mSizeLow1 = mSizeLow
		mSizeHigh1 = mSizeHigh
		mSizeStep1 = mSizeStep
		if mSizeType==0:
			mSizeLow1 = n*mSizeLow
			mSizeHigh1 = n*mSizeHigh
			mSizeStep1 = n*mSizeStep
			print n,mSizeLow,mSizeHigh,mSizeStep
		for k in range(int(mSizeLow1),int(mSizeHigh1),int(mSizeStep1)):
			for i in range(instances):
				opStr = "p cnf "+str(n)+" "+str(int(n*m))+'\n'
				for j in range(int(n*m)):

					sampVars = random.sample(xrange(1,n+1),k)
					if (j%200 == 0):
						sampVars = random.sample(xrange(1,n+1),min(n, 300*k))
					for l in sampVars:
						if monotone == 1:
							opStr += str(l)+" "
						else:
							if random.random()>0.5:
								opStr += str(l)+" "
							else:
								opStr += "-"+str(l)+" "
					opStr += "0\n"
				f = open(outputDir+"/randomDNF_"+str(n)+"_"+str(int(n*m))+"_"+str(k)+"_"+str(i)+".dnf",'w')
				f.write(opStr)
				f.close()

		m += mDensityStep

