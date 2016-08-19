'''

This is a program calculating the ground state energy of 
a small system with Full CI) method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy as np
from pdb import set_trace as stop

def HC(Cx, k_mtx, g_mtx, K, G, E_pq, G_pq, N, string_data):
    ne,no,ns = N
    occ,vir,aclist,Jlist,signlist = string_data
    # unpacking
    C0 = Cx.reshape(ns,ns)
    ###########################################################
    sig = np.dot(K,C0) + np.dot(C0,K.T)
    del K
    ###########################################################
    sig += 0.5 * (np.dot(G,C0) + np.dot(C0,G.T))   
    del G
    #==================
    D_pq = np.einsum('pqij,jk->pqik',E_pq,C0)
    del E_pq
    sig22pq = np.einsum('pqij,pqkj->ik',D_pq,G_pq)
    sig += sig22pq
    result = sig.reshape(-1,1)
    #-----------------------------------------------------------------------
    return result 

