import sys
sys.path.append("../../library_python")
import Hamiltonian_light_matter as Hlm
import bosonic_operators as bo
import su2_spins_operators as su2
from quimb.tensor import *
from quimb import *
from quimb.evo import *
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math

# Dynamics dissipative (photon losses) matter interacting via a single-mode cavity (bosonic)
# It is possible to select either Dicke or Tavis-Cummings coupling via the string variable
# photon_matter_coupling concerning the unitary part.
# photon_matter_coupling = "dicke"
#       H = omega0 a^dag a + h \sum_j \sigma_j^z + g/sqrt{N} (a+a^\dag)\sum_j \sigma_j^x
#         = omega0 a^dag a + h S^z + g/sqrt{N} (a+a^\dag) S^x

# The whole dynamics is obtained by integrating the Lindbland master equation
# \dot{rho} = -i[H,\rho] + k (a \rho a^dag - a^\dag a \rho - \rho a^\dag a),

# As the system respects permutation symmetry, we can work with collective operators.
# Specifically, S^alpha are SU(2) spin-N/2 operators.


folder = "/home/ricval/Documenti/projects_in_progress/dicke_model_correlated/unitary_dicke_model/dicke_model_test_ED_vs_TN"

fig, axs = plt.subplots(1,2,figsize=(12,6))


def mine_lindblad_eq(t,y,args):

    ham     = args[0]
    ls      = args[1]
    gamma   = args[2]

    d = ham.shape[0]
    ham_sparse = issparse(ham) or sparse
    idt = eye(d, sparse=ham_sparse)
    evo_superop = -1.0j * ((ham & idt) - (idt & ham.T))

    def gen_lb_terms():
        for l in ls:
            lb_sparse = issparse(l) or sparse
            idt = eye(d, sparse=lb_sparse)
            yield ((l & l.conj()) - 0.5 * ((idt & dot(dag(l), l).T) + (dot(dag(l), l) & idt)))
                                           
    evo_superop += gamma * sum(gen_lb_terms())
    return dot(evo_superop, y)

# parameters
sparse = False
N = 3
MaxOcc = 2
om_ph  = 1
om_z   = 1
kappa  = 1
gc     = 0.5 * np.sqrt( om_z * (om_ph**2 + kappa**2) / om_ph )

args = {}
args['N'] = N
args['MaxOcc'] = MaxOcc
args['omega_ph'] = om_ph
args['omega_z']  = om_z

g_list = [gc * j / 10 for j in range(20)]
g_list = [1.5]
g_list = [3]
# dynamics parameters
dt  = 0.01
T   = 15
ts = np.linspace(dt,T,num=int(T/dt))
dims = [MaxOcc+1 , N+1] 

# observables
Nb = bo.Nb(MaxOcc,sparse)
Sx = su2.SU2_Sx(N/2,sparse)

# initiail state
psi_ph = qu([1] + [0 for j in range(MaxOcc)], qtype='ket', sparse=sparse)
_ , eigenvectors = eigh(-Sx)
psi_sp = qu(eigenvectors[:,0], qtype='ket', sparse=sparse)
psi_t0 = psi_ph & psi_sp


Sx = ikron(Sx,dims=dims,inds=[1])
Nb = ikron(Nb,dims=dims,inds=[0])

np_infty = []

for g in g_list:

    args['g'] = g
    H1 = Hlm.H_dicke(args,sparse)

    # list where to save results
    sx_t   = []
    np_t   = []
    data_t = []
    
    if kappa <= 1E-10:
        evo = Evolution(psi_t0, H1, progbar=True,method='integrate')
        for idx, psi_t in enumerate(evo.at_times(ts)):
            np_t += [np.real(expec(psi_t,Nb))/N]
            sx_t += [np.real(expec(psi_t,Sx))/N]
            data_t += [[ts[idx],sx_t[-1],np_t[-1]]]
            
    else:
        # build up the density matrix
        rho_t0 = outer(psi_t0,psi_t0.H)
        rho_t0 /= trace(rho_t0)
        D = prod(dims)

        l = bo.Ab(MaxOcc,sparse)
        Lj = [ikron(l, dims, inds = [0])]

        # setting the integrator
        integrator_solver = ode(mine_lindblad_eq).set_integrator('zvode', method='bdf',rtol=1E-10)
        integrator_solver.set_f_params([H1,Lj,kappa])
        integrator_solver.set_initial_value(rho_t0.reshape(D**2))

        # integration
        tau = dt
        n_measure = 10

        for idx_t, t in enumerate(ts):
            print(t)    
            rhot = integrator_solver.integrate(t)
            rhot = rhot.reshape(D,D)

            rhot /= trace(rhot)

            np_t += [np.real(expec(rhot,Nb))/N]
            sx_t += [np.real(expec(rhot,Sx))/N]
            data_t += [[t,sx_t[-1],np_t[-1]]]

    np.savetxt(f"{folder}/obs_dicke_ED_N{N}_maxocc{MaxOcc}_omega{om_ph:.2f}_h{om_z:.2f}_g{g:.2f}_kappa{kappa:.2f}.txt",data_t)

#     axs[0].plot(ts,np_t,linestyle='-',label=f'{g:.2f}')
#     np_infty += [np_t[-1]]
# axs[0].legend()
# axs[1].plot(g_list,np_infty,marker='o',linestyle='--',color='black')
# axs[1].plot(gc,0,marker='x',color='tab:red')
# plt.show()

# folder = "test_ising/"
# folder = "/home/ricval/Documenti/projects_in_progress/edge_modes_tfic/data/non_integrable/ED_data/"
# np.savetxt(f'{folder}ED_Ising_Mixed_N{L}_Jxx{Jxx:.3f}_hx{hx:.3f}_gamma{gamma:.3f}_exact{exact}.txt',data_t)


