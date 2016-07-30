import numpy
import operator  

def math_C(n,k):  
        return reduce(operator.mul, range(n - k + 1, n + 1)) /reduce(operator.mul, range(1, k +1))  

def getOccu(ne,no,l,Z):
    #-------------
    # Explanation: shavitt map
    # 1 - 1 - 1 - 1 ne+1
    # |   |   |   |
    # 1 - 2 - 3 - 4
    # |   |   |   |
    # 1 - 3 - 6 - 10*
    # no-ne+1
    #--------------
    pi, pj = no-ne, ne
    occu = numpy.zeros(no)
    for i in xrange(no):
        if pi == 0: # at top line
            pj -= 1
            occu[no-1-i] = 1

        else:
            if l >= Z[pi-1,pj]:
                l -= Z[pi-1,pj]
                occu[no-1-i] = 1
                pj -= 1 # move left
            else:
                pi -= 1 # move up
        #print pi, ',', pj, 'l=', l
    return occu

def getOrder(occu,Z):
    c = 0
    pi, pj = 0, 0
    for i in occu:
        if i == 1:
            pj += 1 # move right
            if pi > 0:
                c += Z[pi-1,pj]
        else:
            pi += 1 # move down
        #print pi, ',', pj, 'c=', c
    return c


def opr_E(ne,no,r,s,k,l,Z):
    occu = getOccu(ne,no,l,Z) 
    #Ers = ar+ . as
    if occu[s] == 0:
        return 0
    else:
        occu[s] = 0    
    if occu[r] == 1:
        return 0
    else:
        occu[r] = 1
    #print occu

    c = getOrder(occu,Z)
    if c == k:
        return 1
    else:
        return 0

def getStr(ne,no,l,Z):
    spstr = numpy.zeros(ne)
    pi, pj = no-ne, ne
    ei = ne - 1
    for i in xrange(no):
        if pi == 0: # at top line
            pj -= 1
            spstr[ei] = no - 1 - i
            ei -= 1
        else:
            if l >= Z[pi-1,pj]:
                l -= Z[pi-1,pj]
                occu[no-1-i] = 1
                spstr[ei] = no - 1 - i
                ei -= 1
                pj -= 1 # move left
            else:
                pi -= 1 # move up
    return spstr

def opr_EE(ne,no,p,q,r,s,k,l,Z):
    occu = getOccu(ne,no,l,Z) 
    #Ers = ar+ . as
    if occu[s] == 0:
        return 0
    else:
        occu[s] = 0    
    if occu[r] == 1:
        return 0
    else:
        occu[r] = 1

    if occu[q] == 0:
        return 0
    else:
        occu[q] = 0    
    if occu[p] == 1:
        return 0
    else:
        occu[p] = 1
    #print occu

    c = getOrder(occu,Z)
    if c == k:
        return 1
    else:
        return 0

def constructZ(ne,no):
    Z = numpy.ones([no-ne+1,ne+1])
    #Full CI
    for i in xrange(1,no-ne+1):
        for j in xrange(1,ne+1):
            Z[i,j] = Z[i-1,j] + Z[i,j-1]
    return Z

def constructE(ne,no,ns):
    Z = constructZ(ne,no)
    E = numpy.zeros([ns,ns])
    for i in xrange(ns):
        for j in xrange(ns):
            E(i,j)

if __name__ =='__main__':
    Z = constructZ(1,2)
    print Z
    print opr_E(1,2,0,0,0,0,Z)
    
