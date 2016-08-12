import numpy as np
import operator
import itertools as it

def math_C(n,k):  
        return reduce(operator.mul, range(n - k + 1, n + 1)) / reduce(operator.mul, range(1, k +1))  

def construct_string_data(N):
    ne,no,ns = N
    # unpack numbers

    occ = tuple(it.combinations(range(no),ne))
    vir = tuple(tuple(set(range(no)).difference(i)) for i in occ)
    # occ and vir for i_th string

    spstr = []
    # i_th spin string
    dict_sps2i = {}
    # spin string to its order i
    for i in occ:
        sps_tmp = 0
        for j in i:
            sps_tmp |= 1 << j
        spstr.append(sps_tmp)
        dict_sps2i[sps_tmp] = i

    aclist = [tuple(it.product(occ[i],vir[i])) for i in xrange(ns)]
    # (annihilation, creation)
    dict_
    for i, ac in enumerate(aclist):
        str0 = spstr[i]
        for p,q in ac:
            str1 = str0 ^ (1 << p | 1 << q)
            
    return occ, vir, sps

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

