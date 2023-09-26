from mprandom import vanbinomialvariate as vbn, mpbinomial as mpb
from decimal import Decimal
import time 

N = Decimal(2** 7000)

p = Decimal(1)
for i in range(5000):
    p /= 2

f = open("timeFile", "w")

toWrite = ""

xtimes = 1 

for b in range(1000, 7001, 200):
    toWrite += str(b) +"\t"
    for pb in range(500, b+1, 500):
        totimevbn, totimempb = 0, 0
        p = Decimal(1)
        for i in range(pb):
            p /= 2
        N = Decimal(2**b)
        for _ in range(xtimes):
            start_time = time.time()
            # vbn(N, p)
            print("vbn", vbn(N, p))
            totimevbn += (time.time() - start_time)
            start_time = time.time()
            # mpb(N, p)
            print("mpb",mpb(N, p))
            totimempb += (time.time() - start_time)

        toWrite += str(totimevbn / xtimes) + "\t" + str(totimempb / xtimes) + "\t"
    toWrite += "\n"

f.write(toWrite)