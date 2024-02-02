from math import pi, sqrt, e, comb, factorial

def stirling(n):
    return sqrt(2 * pi * n) * (n / e) ** n


def bindist(n, p):
    dist1 = {}
    for k in range(n+1):
        dist1[k] = comb(n, k) * p**k * (1-p)**(n-k)
    
    dist2 = {}
    dist2[0] = (1-p)**n
    distsum = dist2[0]
    for k in range(1,n):
        try:
            dist2[k] = stirling(n) / (stirling(k) * stirling(n-k)) * p**k * (1-p)**(n-k)
            distsum += dist2[k]
        except ZeroDivisionError:
            print("ERROR:", stirling(k), stirling(n-k))
            exit
    dist2[n] = p**n
    distsum += dist2[n]

    fk_sum = distsum

    # for k in range(1+n):
    #     dist2[k] = dist2[k] / distsum

    for k in range(1, n):
        print(k, "->",  dist2[k], " LB:", dist1[k] * e**( 1/(12*(n-k)+1) + 1/(12*k+1) - 1/(12*n)), "UB:", dist1[k] * e**(1/(12*(n-k)) + 1/(12*k) - 1/(12*n+1)))
        assert dist2[k] >= dist1[k] * (1 + 1/(12*(n-k)+1) + 1/(12*k+1) - 1/(12*n)) and dist2[k] <= dist1[k] * (1 + 2/(12*(n-k)) + 2/(12*k) - 2/(12*n+1))

    dist3 = {}
    dist3[0] = (1-p)**n
    distsum = dist3[0]
    for k in range(1,n):
        try:
            dist3[k] = dist1[k] * e**(1/(12*(n-k)) + 1/(12*k) - 1/(12*n))
            distsum += dist3[k]
        except ZeroDivisionError:
            print("ERROR:", stirling(k), stirling(n-k))
            exit
    dist3[n] = p**n
    distsum += dist3[n]

    for k in range(1+n):
        dist3[k] = dist3[k] / distsum
    print("distsum", distsum)
    
    return dist1, dist2, fk_sum

if __name__=='__main__':
    # for n in range(1, 1000)
    n = 150
    p = 0.01
    q = 1-p

    for k in range(1,n):
        # print("UB:", stirling(k) * e**(1/(12*k)))
        # print("LB:", stirling(k) * e**(1/(12*k+1)))
        # print("k! : ", factorial(k))
        # print(k, "->", comb(n,k), " LB:", (stirling(n) / (stirling(n-k) * stirling(k))) * e**( -1/(12*(n-k)) - 1/(12*k) + 1/(12*n+1)), "UB:", (stirling(n) / (stirling(n-k) * stirling(k))) * e**(-1/(12*(n-k)+1) - 1/(12*k+1) + 1/(12*n)))
        # print(k, "->", (stirling(n) / (stirling(n-k) * stirling(k))), " LB:", comb(n,k) * e**( 1/(12*(n-k)+1) + 1/(12*k+1) - 1/(12*n)), "UB:", comb(n,k) * e**(1/(12*(n-k)) + 1/(12*k) - 1/(12*n+1)))
        print(k, "->", (stirling(n) / (stirling(n-k) * stirling(k))) * p**k * (1-p)**(n-k), " LB:", comb(n,k) * e**( 1/(12*(n-k)+1) + 1/(12*k+1) - 1/(12*n)) * p**k * (1-p)**(n-k), "UB:", comb(n,k) * e**(1/(12*(n-k)) + 1/(12*k) - 1/(12*n+1)) * p**k * (1-p)**(n-k))
        

    print("\n\n\n")

    dist1, dist2, dist2_sum = bindist(n,p)
    l1 = 0
    for i in range(n+1): #list(dist1.keys()):
        l1 += abs(dist1[i]-dist2[i])
        print("e->", e**( 1/(12*(n-i)+1) + 1/(12*i+1) - 1/(12*n)))
    
    print('*'*100 , '\n', dist1, '\n', '*'*100 , '\n', '*'*100 , '\n', dist2, '\n', '*'*100)
    
    lmn = (1-p**(n+2)-(1-p)**(n+2)) / ((n+1)*p*(1-p)) - p**(n+1) - (1-p)**(n+1) - p**(n+1)/(n+1) - q**(n+1)/(n+1)
    # print("esdfs", 1/(1-p) + 1/((n+1)*(1-p)))
    M = lmn/12 - 1/(12* n) + 1
    N = lmn - 1/(12* n) + 1

    # print("1/npq", lmn/ (1-1/(12*n)+lmn), "M:", M, "N:", N)
    # print("The l1 distance:", l1, "  UB:", lmn * 2/M, "  LB:", lmn * (2/12)/N)
    # print("1/nq", 1 / ((n+1) * (1-p)))

    fk, fnk = 0, 0
    fk_ =0 
    for k in range(n+1):
        fk += dist1[k] / (k + 1)
        fnk += dist1[k] / (n - k + 1)
        fk_ += dist2[k]

    # print(M, "<", dist2_sum, "<", N)

    # print(fk/12 + fnk/12 - 1/(12*n) + 1, (1 - p**(n+1))/((n+1)*(1-p)))

    # print("lmn", lmn)

    l1_ = 0; l1 = 0; l1_2 = 0
    for i in range(1, n+1):
        l1 += abs(dist1[i] - dist2[i])
        l1_ += dist1[i]/M * abs(lmn - 1/(12 * n) + 1 - (1+ 1/(12*(n-i)+1) + 1/(12*i+1) - 1/(12*n)))
        # l1_2 += 
    #     print("lsdf", dist1[i]/dist2_sum * abs(dist2_sum - e**( 1/(12*(n-i)+1) + 1/(12*i+1) - 1/(12*n))))
    #     print("sd", abs(dist1[i] - dist2[i]))
    #     print("l1_", l1_)
    # print("l1", l1, "\nl1_", l1_, dist2_sum)

    print("l1", l1)

    # print("fk", fk, "1/npq", (1-p**n-(1-p)**n)/(n*p*(1-p)), 1-p, (1-p)**n)

    # print("np/q", n*p/q)