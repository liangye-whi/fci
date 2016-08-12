import numpy as np
import itertools as it
import operator

def math_C(n,k):  
        return reduce(operator.mul, range(n - k + 1, n + 1)) / reduce(operator.mul, range(1, k +1))  

def sign(string0, string1):
    ss = string1 - string0
    def count_bit1(n):
        # see Hamming weight problem and K&R C program
        return bin(n).count('1')
    if ss > 0:
        # string1&ss gives the number of 1s between two strings
        return (-1) ** (count_bit1(string1&ss))
    elif ss == 0:
        return 1
    else:
        return (-1) ** (count_bit1(string0&(-ss)))

def construct_string_data(N):
    ne,no,ns = N
    # unpack numbers

    occ = tuple(it.combinations(range(no),ne))
    vir = tuple(tuple(set(range(no)).difference(i)) for i in occ)
    # occ and vir for i_th string

    spstr = []
    # i_th spin string
    sps2i = {}
    # spin string to its order i
    for i,iocc in enumerate(occ):
        sps_tmp = 0
        for j in iocc:
            sps_tmp |= 1 << j
        spstr.append(sps_tmp)
        sps2i[sps_tmp] = i

    aclist = [tuple(it.product(occ[i],vir[i])) for i in xrange(ns)]
    # (annihilation, creation)
    Jlist = []
    for i,iac in enumerate(aclist):
        J0 = []
        str0 = spstr[i]
        for a,c in iac:
            str1 = str0 ^ (1<<a|1<<c)
            J0.append(sps2i[str1])
        Jlist.append(tuple(J0))

    for i, ac in enumerate(aclist):
        str0 = spstr[i]
        for p,q in ac:
            str1 = str0 ^ (1 << p | 1 << q)
            
    signlist = []
    for i,jlist in enumerate(Jlist):
        s_list = []
        for j in jlist:
            s_list.append(sign(spstr[i],spstr[j]))
        signlist.append(s_list)

    return occ, vir, spstr, sps2i, aclist, Jlist, signlist

