"""
    @Author: Yannick Dengler
    @Date:   2024-Jul-18
    
    This file contains scripts to calculate the phase shift from energy levels (obtained from "energy_levels"). We work with unitless quantities in lattice units, so everything that is named "_prime" is scaled with the pion mass. 
    Input: E_tot, dvec, N_L, m_inf
    Execute with: "python3 scattering/scattering.py"
 """
 
import numpy as np
from gz_leskovec_module.generalized_zeta import generalized_zeta as gz   # input as (l,m,d1,d2,d3,m1,m2,q2,N_L,size=10)

def q_gamma(E_tot, dvec, N_L, m_inf):
    print(dvec[0],dvec[1],dvec[2])
    P_tot2 = (2*np.pi/N_L)**2*(dvec[0]**2+dvec[1]**2+dvec[2]**2)
    E_CM = np.sqrt(E_tot**2-P_tot2)
    Pstar = np.sqrt(E_CM**2/4 - m_inf**2)
    q = Pstar*N_L/(2*np.pi)
    gamma = E_tot/E_CM
    # print(E_CM, Pstar)
    # print(q, q**2,gamma)
    # print(2*np.sqrt(m_inf**2+(2*np.pi*q/N_L)**2))
    return q, gamma

def tan_PS_T1m(q):
    return np.pi**(3/2)*q/