import numpy as np
from quimb import *

# bosonic operators

def Ab(maxocc,sparse=False):
    return qu(np.diag([np.sqrt(j) for j in range(1,maxocc+1)],k=1),qtype='dop',sparse=sparse).real.astype(float)


def Abd(maxocc,sparse=False):
    # Ad  = qu(np.diag([np.sqrt(j) for j in range(1,MaxOcc+1)],k=-1),qtype='dop',sparse=sparse).real.astype(float)
    return np.transpose(np.conj(Ab(maxocc,sparse)))

def Nb(maxocc,sparse=False):
    A  = Ab(maxocc,sparse)
    Ad = Abd(maxocc,sparse)
    Nb = Ad @ A

    Nb_ = qu(np.diag([j for j in range(maxocc+1)]),qtype='dop',sparse=sparse).real.astype(float)

    if np.amax(np.abs(Nb-Nb_)) > 1E-10:
        print("Mismatch between definitions of number operator")
        exit(0)
    
    return Nb