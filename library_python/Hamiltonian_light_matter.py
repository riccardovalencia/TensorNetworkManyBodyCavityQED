import numpy as np
from quimb import *
import su2_spins_operators as su2
import bosonic_operators as bo

# dicke Hamiltonian

def H_dicke(args,sparse=False):
 
    N = args['N']
    MaxOcc = args['MaxOcc']
    omega_ph = args['omega_ph']
    omega_z  = args['omega_z']
    g        = 2 * args['g'] / np.sqrt(N)

    Sx = su2.SU2_Sx(N/2,sparse)
    Sz = su2.SU2_Sz(N/2,sparse)

    Id_matter = qu(np.eye(int(N+1)),qtype='dop',sparse=sparse)
    Nph = bo.Nb(MaxOcc,sparse)
    Xph = bo.Ab(MaxOcc,sparse) + bo.Abd(MaxOcc,sparse)
    Iph = qu(np.eye(MaxOcc+1),qtype='dop',sparse=sparse)

    H  = omega_ph * Nph & Id_matter
    H += omega_z  * Iph & Sz
    H += g        * Xph & Sx

    return H
