import numpy
import operator  

def math_C(n,k):  
        return reduce(operator.mul, range(n - k + 1, n + 1)) /reduce(operator.mul, range(1, k +1))  

def formOccu(ne,no,ns,Z):
	# Build the complete dictionary for index to occupation
    #-------------
    # Explanation: shavitt map
    # 1 - 1 - 1 - 1 ne+1
    # |   |   |   |
    # 1 - 2 - 3 - 4
    # |   |   |   |
    # 1 - 3 - 6 - 10*
    # no-ne+1
    #--------------
    occu = [[0]*ns,{}]
    for l in xrange(ns):
        pi, pj = no-ne, ne
        cu = 0
        ll = l
        #print ll
        for i in xrange(no):
            if pi == 0: # at top line
                cu = cu << 1
                cu += 1
                pj -= 1

            else:
                if ll >= Z[pi-1,pj]:
                    ll -= Z[pi-1,pj]
                    cu = cu << 1
                    cu += 1
                    pj -= 1 # move left
                else:
                    cu = cu << 1
                    pi -= 1 # move up
            #print pi, ',', pj, 'll=', ll
        res = bin(cu)[:1:-1]
        for i in xrange(no-len(res)): res+='0'
        occu[0][l] = int(res,2) # index to occu, e.g. 19 = 0b10011
        #print int(res,2)
        occu[1][str(occu[0][l])] = l # occu to index
        #print res
        assert sum([int(i) for i in res if i == '1']) == ne ##debug##
    #print occu
    return occu

def constructZ(ne,no):
    Z = numpy.ones([no-ne+1,ne+1])
    #Full CI
    for i in xrange(1,no-ne+1):
        for j in xrange(1,ne+1):
            Z[i,j] = Z[i-1,j] + Z[i,j-1]
    return Z

def sign(res,p,q):
    # even +1; odd -1
    kp = 1 - (int(bin(sum([int(i) for i in res[:p] if i=='1']))[-1]) << 1)
    kq = 1 - (int(bin(sum([int(i) for i in res[:q] if i=='1']))[-1]) << 1)
    if p > q:
        return kp * kq
    elif p < q:
        return kp * kq * (-1)
    else:
        print res
        print 'p,q=',p,q
        print 'error'
        assert 0
        return 1

