import numpy as np
from quimb import *

# tested: YES
def SU2_Sp(n,sparse=False):
    # dimension
    D  = int(2*n+1)
    Sp = qu(np.zeros(shape=(D,D),dtype=float),qtype='dop',sparse=sparse)
    mj = [-n+j for j in range(D)]

    # filling elements
    for j , m in enumerate(mj):
        if j+1 < D:
            Sp[j,j+1] = np.sqrt(n * (n+1) - m * (m+1))
    return Sp

def SU2_Sm(n,sparse=False):
    return np.transpose(np.conj(SU2_Sp(n,sparse)))

def SU2_Sx(n,sparse=False):
    Sp = SU2_Sp(n,sparse)
    Sm = SU2_Sm(n,sparse)
    return (Sp+Sm)/2

def SU2_Sy(n,sparse=False):
    Sp = SU2_Sp(n,sparse)
    Sm = SU2_Sm(n,sparse)
    return (Sp-Sm)/(2j)

# tested: YES
def SU2_Sz(n,sparse=False):
    # dimension
    D  = int(2*n+1)
    Sz = qu(np.zeros(shape=(D,D),dtype=float),qtype='dop',sparse=sparse)
    mj = [n-j for j in range(D)]

    # filling elements
    for j , m in enumerate(mj):
        Sz[j,j] = m
    return Sz