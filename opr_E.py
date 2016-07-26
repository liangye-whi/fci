import numpy
import operator  

def math_C(n,k):  
        return  reduce(operator.mul, range(n - k + 1, n + 1)) /reduce(operator.mul, range(1, k +1))  

def getOrder(ne,no,l,Z):
    pi, pj = no-ne, ne
    order = numpy.zeros(no)
    for i in range(no):
        if pi == 0:
            pj -= 1
            order[no-1-i] = 1

        else:
            if l >= Z[pi-1,pj]:
                l -= Z[pi-1,pj]
                order[no-1-i] = 1
                pj -= 1
            else:
                pi -= 1
        #print pi, ',', pj, 'l=', l
    return order


def opr_E(ne,no,r,s,k,l,Z):
    pi, pj = no-ne, ne
    order = numpy.zeros(no)
    for i in range(no):
        if pi == 0:
            pj -= 1
            order[no-1-i] = 1

        else:
            if l >= Z[pi-1,pj]:
                l -= Z[pi-1,pj]
                order[no-1-i] = 1
                pj -= 1
            else:
                pi -= 1
        #print pi, ',', pj, 'l=', l
    #print order
    
    #Ers = ar+ . as
    if order[s] == 0:
        return 0
    elif order[l] == 1 and s != l:
        return 0
    else:
        order[s] = 0
        order[r] = 1
    #print order

    c = 0
    pi, pj = 0, 0
    for i in order:
        if i == 1:
            pj += 1
            if pi > 0:
                c += Z[pi-1,pj]
        else:
            pi += 1
        #print pi, ',', pj, 'c=', c
    #print c

    if c == k:
        return 1
    else:
        return 0


def constructZ(ne,no):
    Z = numpy.ones([no-ne+1,ne+1])
    #Full CI
    for i in range(1,no-ne+1):
        for j in range(1,ne+1):
            Z[i,j] = Z[i-1,j] + Z[i,j-1]
    return Z

if __name__ =='__main__':
    Z = constructZ(1,2)
    print Z
    print opr_E(1,2,0,0,0,0,Z)
    
