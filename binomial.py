from random import *
import decimal
from decimal import Decimal, Context, getcontext, setcontext, MAX_PREC, MAX_EMAX, MIN_EMIN
from operator import index as _index
from math import log2 as _log2, log10 as _log10, fabs as _fabs, lgamma as _lgamma, log as _log, floor as _floor
from math import sqrt as _sqrt, exp as _exp
import gmpy2 as gp
import time
import matplotlib.pyplot as plt
log = lambda t: gp.log10(gp.mpfr(str(t)))
    

#----------------- Arbitrary Precision Random number generator -----------


def mpuniform(precision=30):
    """ multiprecision uniform number
    """
    precision2 = int(_log2(10 ** precision))
    univariate = getrandbits(precision2)
    getcontext().prec = precision
    univariate = univariate / Decimal(10**precision)
    return univariate

    """multi precision Normal distribution.

        uses Charles F. F. Karney method
        ref: https://arxiv.org/pdf/1303.6257.pdf
    """
    getcontext().prec = precision
    while True:
        k = 0
        # N1
        while not _half_exp_bernoulli(): pass
        while _half_exp_bernoulli():
            k += 1
        # N2
        check = True
        for i in range(k*(k-1)): check *= _half_exp_bernoulli()
        if not check: continue
        #N3
        x = mpuniform(precision = precision)
        #N4
        check = True
        for i in range(k+1): check *= _general_exp_variate(x, (2*k+x)/(2*k+2), precision)
        if check: break
    #N5
    x
    y = Decimal(k) + Decimal(x)
    #N6
    if random() > 0.5 : sgn = 1
    else: sgn = -1
    return Decimal(sgn) * y 

# def _logfact(k):
#     """logfactorial
#     """
#     a = [8.333333333333333e-02, -2.777777777777778e-03,
#         7.936507936507937e-04, -5.952380952380952e-04,
#         8.417508417508418e-04, -1.917526917526918e-03,
#         6.410256410256410e-03, -2.955065359477124e-02,
#         1.796443723688307e-01, -1.39243221690590e+00]
    
#     if k == 1 or k == 2 : 
#         return 0
#     elif k < 7 :
#         n = 7 -k
#     else:
#         n = 0
#     k0 = k + n
#     k2 = (1 / k0) ** 2
#     lg2pi = 1.8378770664093453e+00
#     gl0 = a[9]
#     for i in range(9, 0, -1):
#         gl0 *= k2
#         gl0 += a[i]
#     gl = gl0 / k0 + 0.5 * lg2pi + (k0 - 0.5) * _log(k0) - k0
#     if k < 7:
#         gl -= _log(k0 - 1)
#         k0  -= 1
#     return gl

def lanczpoisson(lambd = 10):
    """PTRS
    """
    start = time.time()
    lambd = gp.mpfr(str(lambd))
    if lambd < 10 :
        exlam = gp.exp(-lambd)
        k = 0
        prod = gp.mpfr(1)
        while True:
            U = gp.mpfr(str(mpuniform()))
            prod *= U
            if prod > exlam:
                k += 1
            else:
                print("lanczpois time taken : ", time.time() - start)
                return k
    elif lambd >= 10 :
        lnlam = gp.log(lambd)
        b = 0.931 + 2.53 * gp.sqrt(lambd)
        a = - 0.059 + 0.02483 * b
        vr = 0.9277 - 3.6224 / (b - 2)
        invalpha = 1.1239 + 1.1328 / (b - 3.4)
        
        while True:
            U, V = gp.mpfr(str(mpuniform())) - 0.5, gp.mpfr(str(mpuniform()))
            us = 0.5 - _fabs(U)
            k = _floor(( 2 * a / us + b) * U + lambd + 0.43)
            # if (us >= 0.07) and (V <= vr):
            #     print("lanczpois time taken : ", time.time() - start)
            #     return k
            if (k <= 0) or (us < 0.013 and V > us):
                continue
            if (gp.log(V) + gp.log(invalpha) - gp.log(a / us**2 + b)) <= (-lambd + k * lnlam - gp.lgamma(k)[0]):
                print("lanczpois time taken : ", time.time() - start)
                return k 
    
            

def lanczbinom(n=1, p=0.5):
        """Binomial random variable.

        Gives the number of successes for *n* independent trials
        with the probability of success in each trial being *p*:

            sum(random() < p for i in range(n))

        Returns an integer in the range:   0 <= X <= n

        """
        start = time.time()
        gp.set_context(gp.context(precision = 3000))
        p, n = gp.mpfr(str(p)), gp.mpfr(str(n))
        
        # Error check inputs and handle edge cases
        if n < 0:
            raise ValueError("n must be non-negative")
        if p <= 0.0 or p >= 1.0:
            if p == 0.0:
                return 0
            if p == 1.0:
                return n
            raise ValueError("p must be in the range 0.0 <= p <= 1.0")

        # Exploit symmetry to establish:  p <= 0.5
        if p > 0.5:
            return n - lanczbinom(n, 1.0 - p)

        # BTRS: Transformed rejection with squeeze method by Wolfgang Hörmann
        # https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.47.8407&rep=rep1&type=pdf
        setup_complete = False

        spq = gp.sqrt(gp.mul(gp.mul(n, p), (1 - p)))  # Standard deviation of the distribution
        b = gp.add(gp.mpfr('1.15'),  gp.mul(gp.mpfr('2.53'), spq))
        a = gp.add(gp.add(gp.mpfr('-0.0873'), gp.mul(gp.mpfr('0.0248'), b)), gp.mul(gp.mpfr('0.01'), p))
        c = gp.add(gp.mul(n, p), gp.mpfr('0.5'))
        vr = gp.sub(gp.mpfr('0.92'), gp.div(gp.mpfr('4.2'), b))


        while True:

            u = random()
            u -= 0.5
            us = 0.5 - _fabs(u)
            k = int(gp.mul(gp.add(gp.div(gp.mul(2, a), gp.mpfr(str(us))), b), gp.mpfr(str(u))) + c)
            #    print(f"k:{k}")
            if k < 0 or k > n:
                continue

            # The early-out "squeeze" test substantially reduces
            # the number of acceptance condition evaluations.
            v = gp.mpfr(random())
            # if us >= 0.07 and v <= vr:
            #     print("LanczBinom time taken : ", time.time() - start)
            #     return k
            
            #    print(f"u : {u}, us : {us}")

            # Acceptance-rejection test.
            # Note, the original paper errorneously omits the call to log(v)
            # when comparing to the log of the rescaled binomial distribution.
            if not setup_complete:
                alpha = gp.mul(gp.add(gp.mpfr('2.83'), gp.div(gp.mpfr('5.1'), b)), gp.mpfr(spq))
                ratio = p / (1 - p)
                lpq = gp.log2(gp.mpfr(str(ratio))) * gp.mpfr(_log(2)) 
                m = _floor((n + 1) * p)         # Mode of the distribution
                h1 = gp.mpfr(str(m + 1)); h2 = gp.mpfr(str(n - m + 1))
                h = gp.lgamma(h1)[0] + gp.lgamma(h2)[0]
                setup_complete = True           # Only needs to be done once
            v *= alpha / (a / (us * us) + b)
            logv = gp.log2(gp.mpfr(str(v))) * gp.mpfr(_log(2))       # mpfr log
            #    a1 = v.log10() * Decimal(_log(10))     # decimal log
            h_log = gp.sub(h, gp.add(gp.lgamma(gp.mpfr(str(k + 1)))[0], gp.lgamma(gp.mpfr(str(n - k + 1)))[0]))
            #    if _log(v) <= h - _lgamma(k + 1) - _lgamma(n - k + 1) + (k - m) * lpq:
            if logv <= h_log + (k - m) * lpq:
                print("LanczBinom time taken : ", time.time() - start)
                return k


def stirlbinom(n=1, p=0.5):
        """Binomial random variable.

        Gives the number of successes for *n* independent trials
        with the probability of success in each trial being *p*:

            sum(random() < p for i in range(n))

        Returns an integer in the range:   0 <= X <= n

        """
        start = time.time()
        gp.set_context(gp.context(precision = 3000))
        p, n = gp.mpfr(str(p)), gp.mpfr(str(n))
        
        # Error check inputs and handle edge cases
        if n < 0:
            raise ValueError("n must be non-negative")
        if p <= 0.0 or p >= 1.0:
            if p == 0.0:
                return 0
            if p == 1.0:
                return n
            raise ValueError("p must be in the range 0.0 <= p <= 1.0")

        # Exploit symmetry to establish:  p <= 0.5
        if p > 0.5:
            return n - stirlbinom(n, 1.0 - p)
        
        # BTRS: Transformed rejection with squeeze method by Wolfgang Hörmann
        # https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.47.8407&rep=rep1&type=pdf
        #assert n*p >= 10.0 and p <= 0.5
        setup_complete = False

        spq = gp.sqrt(gp.mul(gp.mul(n, p), (1 - p)))  # Standard deviation of the distribution
        b = gp.add(gp.mpfr('1.15'),  gp.mul(gp.mpfr('2.53'), spq))
        a = gp.add(gp.add(gp.mpfr('-0.0873'), gp.mul(gp.mpfr('0.0248'), b)), gp.mul(gp.mpfr('0.01'), p))
        c = gp.add(gp.mul(n, p), gp.mpfr('0.5'))
        vr = gp.sub(gp.mpfr('0.92'), gp.div(gp.mpfr('4.2'), b))


        while True:

            u = random()
            u -= 0.5
            us = 0.5 - _fabs(u)
            k = int(gp.mul(gp.add(gp.div(gp.mul(2, a), gp.mpfr(str(us))), b), gp.mpfr(str(u))) + c)
            #    print(f"k:{k}")
            if k < 0 or k > n:
                continue

            # The early-out "squeeze" test substantially reduces
            # the number of acceptance condition evaluations.
            v = gp.mpfr(random())
            # if us >= 0.07 and v <= vr:
            #     print("StirlBinom time taken : ", time.time() - start)
            #     return k
            
            #    print(f"u : {u}, us : {us}")

            # Acceptance-rejection test.
            # Note, the original paper errorneously omits the call to log(v)
            # when comparing to the log of the rescaled binomial distribution.
            if not setup_complete:
                alpha = gp.mul(gp.add(gp.mpfr('2.83'), gp.div(gp.mpfr('5.1'), b)), gp.mpfr(spq))
                ratio = p / (1 - p)
                lpq = gp.log2(gp.mpfr(str(ratio))) * gp.mpfr(_log(2)) 
                m = _floor((n + 1) * p)         # Mode of the distribution
                h1 = gp.mpfr(str(m + 1)); h2 = gp.mpfr(str(n - m + 1))
                h = gp.log(gp.sqrt(2 * gp.const_pi())) + gp.mul((h1 + 1/2),gp.log(h1)) - gp.mpfr(h1) + gp.log(gp.sqrt(2 * gp.const_pi())) + gp.mul((h2 + 1/2),gp.log(h2)) - gp.mpfr(h2)
                #h = gp.lgamma(h1)[0] + gp.lgamma(h2)[0]
                setup_complete = True           # Only needs to be done once
            v *= alpha / (a / (us * us) + b)
            logv = gp.log2(gp.mpfr(str(v))) * gp.mpfr(_log(2))       # mpfr log
            strilk = gp.log(gp.sqrt(2 * gp.const_pi())) + gp.mul((gp.mpfr(str(k)) + 1/2),gp.log(gp.mpfr(str(k)))) - gp.mpfr(gp.mpfr(str(k)))
            strilnk1 = gp.log(gp.sqrt(2 * gp.const_pi())) + gp.mul((gp.mpfr(str(n)) - gp.mpfr(str(k)) + 3/2),gp.log(gp.mpfr(str(n)) - gp.mpfr(str(k)) + 1)) - gp.mpfr(gp.mpfr(str(n)) - gp.mpfr(str(k)) + 1)
            h_log = gp.sub(h, gp.add(strilk, strilnk1))
            # h_log = gp.sub(h, gp.add(gp.lgamma(gp.mpfr(str(k + 1)))[0], gp.lgamma(gp.mpfr(str(n - k + 1)))[0]))
            if logv <= h_log + (k - m) * lpq:
                print("StirlBinom time taken : ", time.time() - start)
                return k


def bernbinom(n=1, p=0.5):
        """Binomial random variable.

        Gives the number of successes for *n* independent trials
        with the probability of success in each trial being *p*:

            sum(random() < p for i in range(n))

        Returns an integer in the range:   0 <= X <= n

        """
        start = time.time()
        u = gp.mpfr(str(mpuniform()))
        if u <= gp.mul(gp.mpfr(str(n)), gp.mpfr(str(p))):
            print("BernBinom time taken : ", time.time() - start)
            return 1
        else:
            print("BernBinom time taken : ", time.time() - start)
            return 0


def binomialvariate(n=1, p=0.5):

    zeta = 2 * 10**(-10)
    
    sampler = ''

    if n*p < _sqrt(zeta):
        start = time.time()
        k = bernbinom(n,p)
        tottime = time.time()-start
        sampler = 'BernBinom'
    elif n*p*(1-p) > 4 / zeta:
        if p < zeta**2/16 :
            start = time.time()
            k = lanczpoisson(n*p)
            tottime = time.time()-start
            sampler = 'LanczPois'
        else:
            start = time.time()
            k = stirlbinom(n, p)
            tottime = time.time()-start
            sampler = 'StirlBinom'
    else:
        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        sampler = 'LanczBinom'

    print(sampler, ' ', tottime, ' ', k)

    return sampler + ' ' + str(tottime) + ' \n'# + str(k)


def testbinom2(n=1, p=0.5):

    zeta = 2 * 10**(-10)
    
    sampler = ''
    returnstr = ''
    minerr = 1
    if n*p < 1:

        start = time.time()
        k = bernbinom(n,p)
        tottime = time.time()-start
        err = n**2 * p**2
        assert err <=1
        if err < minerr: 
            sampler = 'bernbinom'
            minerr = err
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        
        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4*(n+2)**2 * p / ((n+1)*(1-p))
        if err < minerr: 
            sampler = 'stirlbinom'
            minerr = err
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        
        start = time.time()
        k = lanczpoisson(n*p)
        tottime = time.time()-start
        err = 2 * n * p**2 + 2 * n * p
        assert err <=1
        if err < minerr: 
            sampler = 'poisbinom'
            minerr = err
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        
        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        err = 30 * zeta
        assert err <=1
        if err < minerr: 
            sampler = 'lanczbinom'
            minerr = err
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        returnstr += ' ' + sampler + '\n'
    
    else:
        
        returnstr += '0 0 '
        
        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4/ (n * p *(1-p))
        assert err <=1
        if err < minerr: 
            sampler = 'stirlbinom'
            minerr = err
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        
        if n*p**2 > 1:
            returnstr += '0 0 '
        else:
            start = time.time()
            k = lanczpoisson(n*p)
            tottime = time.time()-start
            err = 2 * n * p**2 + 2 /( n * p)
            print(n*p**2, err)
            assert err <=1
            if err < minerr: 
                sampler = 'poisbinom'
                minerr = err
            returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '

        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        assert err <=1
        if err < minerr: 
            sampler = 'lanczbinom'
            minerr = err
        err = 30 * zeta
        returnstr += str(tottime) + ' ' + str(gp.log10(gp.mpfr(str(err)))) + ' '
        returnstr += ' ' + sampler + '\n'

    return returnstr

def testbinom(n=1, p=0.5, reg = 'r1'):

    zeta = 2 * 10**(-10)

    times = []
    errs = []

    if reg == 'r1':

        start = time.time()
        k = bernbinom(n,p)
        tottime = time.time()-start
        err = n**2 * p**2
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4*(n+2)**2 * p / ((n+1)*(1-p))
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczpoisson(n*p)
        tottime = time.time()-start
        err = 2 * n * p**2 + 2 * n * p
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        err = 30 * zeta
        times.append(tottime)
        errs.append(err)

        return list(map(log,errs)), times


    elif reg == 'r2':

        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4 / (n*p*(1-p))
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        err = 30 * zeta
        times.append(tottime)
        errs.append(err)

        return list(map(log,errs)), times
    
    elif reg == 'r3':

        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4 / (n * p *(1-p))
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczpoisson(n*p)
        tottime = time.time()-start
        err = 2 * n * p**2 + 2 / (n * p)
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        err = 30 * zeta
        times.append(tottime)
        errs.append(err)

        return list(map(log,errs)), times 
            
    elif reg == 'r4':

        start = time.time()
        k = stirlbinom(n,p)
        tottime = time.time()-start
        err = 4 / (n*p*(1-p))
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczpoisson(n*p)
        tottime = time.time()-start
        err = 2 * n * p**2 + 2 / (n * p)
        times.append(tottime)
        errs.append(err)

        start = time.time()
        k = lanczbinom(n,p)
        tottime = time.time()-start
        err = 30 * zeta
        times.append(tottime)
        errs.append(err)

        return list(map(log,errs)), times

   

if __name__ == '__main__':
    context = Context(prec=3000)
    setcontext(context)
    zeta = 2 * 10**(-10)    

    ncalls = 100
    f = open("results", "w")
    r1_error = {} #bern
    r2_error = {} #stirl
    r3_error = {} #pois
    r4_error = {} #lancz
    
    r1_time = {} #bern
    r2_time = {} #stirl
    r3_time = {} #pois
    r4_time = {} #lancz
    
    for i in range(100, 5000, 100):
        N = 2**i
        p = Decimal(0.5)
        for j in range(0, 10000, 500):
            p /= 2**j
            
            if N*p < _sqrt(zeta):
                r1_error[log(N*p)], r1_time[log(N*p)] = testbinom(N,p,'r1')
            elif N*p*(1-p) > 4 / zeta:
                if p > zeta**2/16 :
                    r2_error[log(N*p)], r2_time[log(N*p)] = testbinom(N,p,'r2')
                else:
                    r3_error[log(N*p)], r3_time[log(N*p)] = testbinom(N,p,'r3')
            else:
                r4_error[log(N*p)], r4_time[log(N*p)] = testbinom(N,p,'r4')


    fig, axs = plt.subplots(4, 1)

    x, y0, y1, y2, y3 = [], [], [], [], []
    for item, values in r1_error.items():
        x.append(item)
        y0.append(values[0])
        y1.append(values[1])
        y2.append(values[2])
        y3.append(values[3])        
    axs[0].plot(x,y0, label = 'bernbinom')
    axs[0].plot(x,y1, label = 'stirlbinom')
    axs[0].plot(x,y2, label = 'poisbinom')
    axs[0].plot(x,y3, label = 'lanczbinom')
    axs[0].set_title('Error in Region - 1')
    axs[0].set(xlabel='np values in $\log_{10}$', ylabel ='Error in $\log_{10}$')
    axs[0].legend()
    axs[0].grid(True)

    x, y0, y1, y2, y3 = [], [], [], [], []
    for item, values in r2_error.items():
        x.append(item)
        y0.append(values[0])
        y1.append(values[1])
    axs[1].plot(x,y0, label = 'stirlbinom')
    axs[1].plot(x,y1, label = 'lanczbinom')
    axs[1].set_title('Error in Region - 1')
    axs[1].set(xlabel='np values in $\log_{10}$', ylabel ='Error in $\log_{10}$')
    axs[1].legend()
    axs[1].grid(True)
    
    x, y0, y1, y2, y3 = [], [], [], [], []
    for item, values in r3_error.items():
        x.append(item)
        y0.append(values[0])
        y1.append(values[1])
        y2.append(values[2])
    axs[2].plot(x,y0, label = 'stirlbinom')
    axs[2].plot(x,y1, label = 'poisbinom')
    axs[2].plot(x,y2, label = 'lanczbinom')
    axs[2].set_title('Error in Region - 3')
    axs[2].set(xlabel='np values in $\log_{10}$', ylabel ='Error in $\log_{10}$')
    axs[2].legend()
    axs[2].grid(True)
    
    x, y0, y1, y2, y3 = [], [], [], [], []
    for item, values in r4_error.items():
        x.append(item)
        y0.append(values[0])
        y1.append(values[1])
        y2.append(values[2])
    axs[3].plot(x,y0, label = 'stirlbinom')
    axs[3].plot(x,y1, label = 'poisbinom')
    axs[3].plot(x,y2, label = 'lanczbinom')
    axs[3].set_title('Error in Region - 4')
    axs[3].set(xlabel='np values in $\log_{10}$', ylabel ='Error in $\log_{10}$')
    axs[3].legend()
    axs[3].grid(True)
    

    # x, y0, y1, y2, y3 = [], [], [], [], []
    # for item, values in r1_time.items():
    #     print(values)
    #     x.append(item)
    #     y0.append(values[0])
    #     y1.append(values[1])
    #     y2.append(values[2])
    #     y3.append(values[3])        
    # axs[0,1].scatter(x,y0, label = 'bernbinom')
    # axs[0,1].scatter(x,y1, label = 'stirlbinom')
    # axs[0,1].scatter(x,y2, label = 'poisbinom')
    # axs[0,1].scatter(x,y3, label = 'lanczbinom')
    # axs[0,1].set(xlabel='np values in $\log_{10}$', ylabel ='Time in $\log_{10}$')
    # axs[0,1].set_title('Time in Region - 1')
    # axs[0,1].legend()
    # axs[0,1].grid(True)


    plt.show()