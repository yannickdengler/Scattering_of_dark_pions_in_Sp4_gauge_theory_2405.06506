"""
    @Author: Yannick Dengler
    @Date:   2024-Mayy-21
    
    This file contains scripts to calculate the phase shift from energy levels (obtained from "energy_levels"). We work with unitless quantities in lattice units, so everything that is named "_prime" is scaled with the pion mass
    Execute with: "python3 scattering/scattering.py"
 """

import h5py
import warnings
import numpy as np
import generalizedzeta as gz
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from tqdm import tqdm



def calculations(N_Ls, E_pis, E_pipis):
    """
    This function calculates the phase shift and related quantities from a list of lattice extends, pion- and two pion energies for one ensembles and returns a dictionary containing the results.
    """
    result = {}
    E_pi_prime = []
    E_pipi_prime = []
    E_pipi1_prime = []
    N_L_prime = []
    P_pipi_lat_prime = []
    P_pipi_prime = []
    P2_pipi_prime = []
    P_cot_PS_pipi_prime = []
    q = []
    q2 = []
    s_pipi_prime = []
    sigma_pipi_prime = []
    num_E = len(N_Ls)
    m_pi_inf, A_R, m_pi_inf_fit = inf_vol_and_mass(N_Ls, E_pis)
    for i in range(num_E):
        N_L_prime.append(N_Ls[i]*m_pi_inf)
        E_pi_prime.append(E_pis[i]/m_pi_inf)
        E_pipi_prime.append(E_pipis[i]/m_pi_inf)
    for i in range(num_E):
        P_pipi_lat_prime.append(P_of_E_pipi_lat(E_pipis[i], m_pi_inf)/m_pi_inf)
    for i in range(num_E):
        tmp = phase_shift(E_pipi_prime[i], N_L_prime[i])
        P_pipi_prime.append(tmp[0])
        P2_pipi_prime.append(tmp[1])
        P_cot_PS_pipi_prime.append(tmp[2])
        q.append(tmp[3])
        q2.append(tmp[4])
    for i in range(num_E):
        s_pipi_prime.append(E_pipi_prime[i]**2)
        sigma_pipi_prime.append(sigma_prime(P_pipi_prime[i], P_cot_PS_pipi_prime[i]))
    a0 = float(fit_phase_shift_0(P2_pipi_prime, P_cot_PS_pipi_prime)[0])
    a2, b2 = fit_phase_shift(P2_pipi_prime, P_cot_PS_pipi_prime) 
    try:
        a_Ad_free, c_Ad_free = fit_phase_shift_Adler_c_free(P2_pipi_prime, P_cot_PS_pipi_prime)
    except RuntimeError:
        a_Ad_free, c_Ad_free = [np.NaN,np.NaN]
    try:
        a_Ad_fixed, c_Ad_fixed = fit_phase_shift_Adler_c_fixed(P2_pipi_prime, P_cot_PS_pipi_prime) 
    except RuntimeError:
        a_Ad_fixed, c_Ad_fixed = [np.NaN,np.NaN]
    Adler_fixed_inter_P_cot_PS = get_interpolation_points_Adler_fixed(a_Ad_fixed,c_Ad_fixed)
    Adler_free_inter_P_cot_PS = get_interpolation_points_Adler_free(a_Ad_free,c_Ad_free)
    UTE_inter_P_cot_PS = get_interpolation_points_UTE(a2,b2)
    sigma_inter_sigma = get_interpolation_points_sigma(a2,b2)
    result["N_Ls"] = N_Ls
    result["m_pi_inf"] = [m_pi_inf,]
    result["A_R"] = [A_R,]
    result["m_pi_inf_fit"] = [m_pi_inf_fit,]
    result["P_pipi_lat_prime"] = P_pipi_lat_prime
    result["E_pi_prime"] = E_pi_prime
    result["E_pipi_prime"] = E_pipi_prime
    result["P_cot_PS_pipi_prime"] = P_cot_PS_pipi_prime
    result["q"] = q
    result["q2"] = q2
    result["P_pipi_prime"] = P_pipi_prime
    result["P2_pipi_prime"] = P2_pipi_prime
    result["s_pipi_prime"] = s_pipi_prime
    result["sigma_pipi_prime"] = sigma_pipi_prime
    result["a0"] = [a0,]
    result["a2"] = [a2,]
    result["b2"] = [b2,]
    result["a_Ad_free"] = [a_Ad_free,]
    result["c_Ad_free"] = [c_Ad_free,]
    result["a_Ad_fixed"] = [a_Ad_fixed,]
    result["c_Ad_fixed"] = [c_Ad_fixed,]
    result["UTE_inter_P_cot_PS"] = UTE_inter_P_cot_PS
    result["sigma_inter_sigma"] = sigma_inter_sigma
    result["Adler_free_inter_P_cot_PS"] = Adler_free_inter_P_cot_PS
    result["Adler_fixed_inter_P_cot_PS"] = Adler_fixed_inter_P_cot_PS
    return result

def inf_mass_fit_Goldstone(N_L, m_inf, A):
    return m_inf*(1+A*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))

def infinite_volume_pi(N_Ls, E_pis):
    popt, pcov = curve_fit(inf_mass_fit_Goldstone,N_Ls,E_pis)
    return popt

def inf_vol_and_mass(N_Ls, E_pis):
    """
    Calculates the infinite volume from extrapolation. Further calculations use min(m_pi_inf_fit, min(E_pis))
    """
    popt = infinite_volume_pi(N_Ls, E_pis)
    m_pi_inf_fit, A_R = popt[0],popt[1]
    return min(m_pi_inf_fit, min(E_pis)), A_R, m_pi_inf_fit

def P_of_E_pipi_lat(E_pipi, mass_meson):
    """
    Lattice dispersion relation. Only used for consistency checks
    """
    arg = 0.5*(np.cosh(E_pipi/2)-np.cosh(mass_meson))
    if arg > 0:
        return 2*np.arcsin(np.sqrt(arg))
    else:
        return np.NaN

def P_of_E_pipi_prime(E_pipi_prime):
    """
    Unitless continuum dispersion relation
    """
    return np.sqrt(E_pipi_prime*E_pipi_prime/4-1)

def q_of_P(P, N_L):
    """
    generalized momentum
    """
    return P*N_L/(2*np.pi)

def phase_shift(E_pipi_prime, N_L):              # all primed
    """
    Calculates the phase shift and related quantities from pipi energy levels and the lattice extend
    """
    if E_pipi_prime < 2 or E_pipi_prime > 4:
        return 1e-40,1e-40,1e-40,1e-40,1e-40
    P = P_of_E_pipi_prime(E_pipi_prime)
    q = q_of_P(P,N_L)
    tan_PS = np.pi**(3/2)*q/gz.Zeta(q*q)
    return P, P**2, P/tan_PS,q, q**2

def sigma_prime(P_prime, P_cot_PS_prime):
    """
    Calculates the cross-section for a single data point.
    """
    return 4*np.pi/(P_cot_PS_prime**2+P_prime**2)

def Adler_A0(c,P2):
    """
    Additional term for Adler-zeros
    """
    return c*np.sqrt(P2+1)/(2*P2+c)

def UTE_A0_c_fixed(P2,a,c):
    """
    Adler-zero fit with c=1
    """
    return Adler_A0(1,P2)*(-1/a+c*P2)

def UTE_A0_c_free(P2,a,c):
    """
    Adler-zero fit with c free
    """
    return -Adler_A0(c,P2)/a

def UTE(P2, a, b):
    """
    Second order univesal threshold expansion
    """
    return a+P2*b

def UTE_0(P2, a):
    """
    First order universal threshold expansion
    """
    return a+P2*0

def fit_phase_shift(P2s, P_cot_PSs):
    """
    fits the phase shift to second order in UTE
    """
    P2_tmp = []
    P_cot_PSs_tmp = []
    for i in range(len(P2s)):
        if P2s[i] != 1e-40:
            P2_tmp.append(P2s[i])
            P_cot_PSs_tmp.append(P_cot_PSs[i])
    if len(P2_tmp) < 2:
        return [1e-40,1e-40]
    warnings.simplefilter("ignore", OptimizeWarning)
    popt, pcov = curve_fit(UTE, P2_tmp, P_cot_PSs_tmp)
    return popt

def fit_phase_shift_0(P2s, P_cot_PSs):
    """
    fits the phase shift to first order in UTE
    """
    popt, pcov = curve_fit(UTE_0, P2s, P_cot_PSs)
    return popt

def fit_phase_shift_Adler_c_fixed(P2s, P_cot_PSs):
    """
    fits the phase shift to Adler-zero with c fixed
    """
    popt, pcov = curve_fit(UTE_A0_c_fixed, P2s, P_cot_PSs)
    return popt

def fit_phase_shift_Adler_c_free(P2s, P_cot_PSs):
    """
    fits the phase shift to Adler-zero with c free
    """
    popt, pcov = curve_fit(UTE_A0_c_free, P2s, P_cot_PSs)
    return popt

# def get_interpolation_points_UTE(a, b, start = 1e-4, stop = 3, num_inter = 2000):
def get_interpolation_points_UTE(a, b, start = -3, stop = 3, num_inter = 2000):
    """
    Interpolates the universal threshold expansion for an error estimate at each value of P
    """
    P2_arr = np.logspace(np.log10(start), np.log10(stop), num_inter)
    P_cot_PS_arr = []
    for P2 in P2_arr:
        P_cot_PS_arr.append(UTE(P2, a, b))
    return np.asarray(P_cot_PS_arr)

# def get_interpolation_points_Adler_fixed(a, c, start = 1e-4, stop = 3, num_inter = 2000):
def get_interpolation_points_Adler_fixed(a, c, start = -3, stop = 3, num_inter = 2000):
    """
    Interpolates the fixed Adler expansion for an error estimate at each value of P
    """
    P2_arr = np.logspace(np.log10(start), np.log10(stop), num_inter)
    P_cot_PS_arr = []
    for P2 in P2_arr:
        P_cot_PS_arr.append(UTE_A0_c_fixed(P2, a, c))
    return np.asarray(P_cot_PS_arr)

# def get_interpolation_points_Adler_free(a, c, start = 1e-4, stop = 3, num_inter = 2000):
def get_interpolation_points_Adler_free(a, c, start = -3, stop = 3, num_inter = 2000):
    """
    Interpolates the fixed Adler expansion for an error estimate at each value of P
    """
    P2_arr = np.logspace(np.log10(start), np.log10(stop), num_inter)
    P_cot_PS_arr = []
    for P2 in P2_arr:
        P_cot_PS_arr.append(UTE_A0_c_free(P2, a, c))
    return np.asarray(P_cot_PS_arr)

def get_interpolation_points_sigma(a, b, start = 4, stop = 15, num_inter = 2000):
    """
    Interpolates the universal threshold expansion for an error estimate at each value of s
    """
    s_arr = np.logspace(np.log10(start), np.log10(stop), num_inter)
    sigma_arr = []
    for s in s_arr:
        P = P_of_E_pipi_prime(np.sqrt(s))
        sigma_arr.append(sigma_prime(P, UTE(P**2,a,b)))
    return sigma_arr

def get_filelist_list(filename, num_lines = 12):
    filelist_list = []
    with open("input/"+filename) as files:
        for i in range(num_lines):
            filelist_list.append(files.readline().split())
    return filelist_list

def save_to_hdf(res,res_sample, filename):
    with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
        for key, val in res.items():
            hfile.create_dataset("orig_"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset("sample_"+key, data = val)
        
def read_from_hdf(filename):
    res, res_tmp = [{},{}]
    with h5py.File("output/hdf5/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:4] == "orig":
                res[key[5:]] = hfile[key][()]
            if key[:4] == "samp":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def result_sampled(beta,m0,N_L,E_pi,E_pi_err,E_pipi,E_pipi_err, num_gaussian=2000):
    res = calculations(N_Ls,E_pi,E_pipi)
    res_sample = {}
    for key in res.keys():
        res_sample[key] = []

    for i in tqdm(range(num_gaussian)):
        #print("%1.1f"%(i*100/num_gaussian)+"%")
        E_pi_tmp, E_pipi_tmp = [[],[]]
        for j in range(len(E_pi)):
            E_pi_tmp.append(np.random.normal(E_pi[j], E_pi_err[j]))
            E_pipi_tmp.append(np.random.normal(E_pipi[j], E_pipi_err[j]))
        res_tmp = calculations(N_Ls,E_pi_tmp,E_pipi_tmp)
        for key, val in res_tmp.items():
            res_sample[key].append(val)
    res["beta"],res["m_1"],res["m_2"],res_tmp["beta"],res_tmp["m_1"],res_tmp["m_2"] = [beta,m0,m0,beta,m0,m0]
    return res, res_sample

if __name__ == "__main__":
    beta_arr = []
    m_arr = []

    N_L_arr = []
    E_pis = []
    E_pi_errs = []
    E_pipis = []
    E_pipi_errs = []


    with h5py.File("output/hdf5/fitresults.hdf5") as file:
        for key in file.keys():
            beta = float(file[key+"/pi/beta"][()])
            m_1 = float(file[key+"/pi/m_1"][()])
            N_L = float(file[key+"/pi/N_L"][()])
            E_pi = float(file[key+"/pi/E"][0])
            E_pi_err = float(file[key+"/pi/Delta_E"][0])
            E_pipi = float(file[key+"/pipi/E"][0])
            E_pipi_err = float(file[key+"/pipi/Delta_E"][0])
            if not beta in beta_arr:
                beta_arr.append(beta)
                m_arr.append([])
                N_L_arr.append([])
                E_pis.append([])
                E_pi_errs.append([])
                E_pipis.append([])
                E_pipi_errs.append([])
            if not m_1 in m_arr[beta_arr.index(beta)]:
                m_arr[beta_arr.index(beta)].append(m_1)
                N_L_arr[beta_arr.index(beta)].append([])
                E_pis[beta_arr.index(beta)].append([])
                E_pi_errs[beta_arr.index(beta)].append([])
                E_pipis[beta_arr.index(beta)].append([])
                E_pipi_errs[beta_arr.index(beta)].append([])
            N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m_1)].append(N_L)
            E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m_1)].append(E_pi)
            E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m_1)].append(E_pi_err)
            E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m_1)].append(E_pipi)
            E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m_1)].append(E_pipi_err)
    for beta in beta_arr:
        for m0 in m_arr[beta_arr.index(beta)]:
            ############## default ####################
            N_Ls = N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)]
            E_pis_t = E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)]
            E_pi_errs_t = E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)]
            E_pipis_t = E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)]
            E_pipi_errs_t = E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)]
            print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
            res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
            fn_tmp = "scattering_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
            save_to_hdf(res, res_sam, fn_tmp)
            res, res_sam = read_from_hdf(fn_tmp)
            ##################################
            ############## Fig 5.1 ####################
            N_Ls = []
            E_pis_t = []
            E_pi_errs_t = []
            E_pipis_t = []
            E_pipi_errs_t = []
            for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
                if N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] > 8:
                    if not (beta == 7.2 and m0 == -0.794 and N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] == 12):
                        N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                        E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                        E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                        E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                        E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            if len(N_Ls) > 1:
                print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
                res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
                fn_tmp = "scattering_Fig5.1_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
                save_to_hdf(res, res_sam, fn_tmp)
                res, res_sam = read_from_hdf(fn_tmp)
            print("Ensembles skipped: beta=",beta,", m0=",m0)
            ##################################
            ############## Fig 5.2 ####################
            N_Ls = []
            E_pis_t = []
            E_pi_errs_t = []
            E_pipis_t = []
            E_pipi_errs_t = []
            for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
                if E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] < 0.95:
                    N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                    E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                    E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                    E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                    E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            if len(N_Ls) > 1:
                print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
                res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
                fn_tmp = "scattering_Fig5.2_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
                save_to_hdf(res, res_sam, fn_tmp)
                res, res_sam = read_from_hdf(fn_tmp)
            print("Ensembles skipped: beta=",beta,", m0=",m0)
            ##################################
            ############## Fig 5.3 ####################
            N_Ls = []
            E_pis_t = []
            E_pi_errs_t = []
            E_pipis_t = []
            E_pipi_errs_t = []
            for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
                if N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] > 8:
                    if not (beta == 7.2 and m0 == -0.794 and N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] == 12):
                        if E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] < 0.95:
                            N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                            E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                            E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                            E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                            E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            if len(N_Ls) > 1:
                print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
                res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
                fn_tmp = "scattering_Fig5.3_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
                save_to_hdf(res, res_sam, fn_tmp)
                res, res_sam = read_from_hdf(fn_tmp)
            print("Ensembles skipped: beta=",beta,", m0=",m0)
            ##################################





            # ############# Epipi < 0.95 ####################
            # N_Ls = []
            # E_pis_t = []
            # E_pi_errs_t = []
            # E_pipis_t = []
            # E_pipi_errs_t = []
            # for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
            #     if E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] < 0.95:
            #         N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #         E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #         E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #         E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #         E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            # if len(N_Ls) > 1:
            #     print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
            #     res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
            #     fn_tmp = "scattering_epipi_095_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
            #     save_to_hdf(res, res_sam, fn_tmp)
            #     res, res_sam = read_from_hdf(fn_tmp)
            # print("Ensembles skipped: beta=",beta,", m0=",m0)
            # #################################
            ############ beta = 6.9, m = -0.90, L > x #####################
            for x in [8,10]:
                if beta == 6.9:
                    if m0 == -0.90:
                        N_Ls = []
                        E_pis_t = []
                        E_pi_errs_t = []
                        E_pipis_t = []
                        E_pipi_errs_t = []
                        for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
                            if N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] > x:
                                N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                                E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                                E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                                E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                                E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
                        if len(N_Ls) > 1:
                            print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
                            res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
                            fn_tmp = "scattering_b69_m90_Lg%i_b%1.3f_m%1.3f"%(x,float(res["beta"]),float(res["m_1"]))
                            save_to_hdf(res, res_sam, fn_tmp)
                            res, res_sam = read_from_hdf(fn_tmp)
            print("Ensembles skipped: beta=",beta,", m0=",m0)
            #################################
            # ############ beta = 7.2, m = -0.78, L > x #####################
            # for x in [8,10]:
            #     if beta == 7.2:
            #         if m0 == -0.78:
            #             N_Ls = []
            #             E_pis_t = []
            #             E_pi_errs_t = []
            #             E_pipis_t = []
            #             E_pipi_errs_t = []
            #             for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
            #                 if N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] > x:
            #                     N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                     E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                     E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                     E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                     E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #             if len(N_Ls) > 1:
            #                 print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
            #                 res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
            #                 fn_tmp = "scattering_b72_m78_Lg%i_b%1.3f_m%1.3f"%(x,float(res["beta"]),float(res["m_1"]))
            #                 save_to_hdf(res, res_sam, fn_tmp)
            #                 res, res_sam = read_from_hdf(fn_tmp)
            # print("Ensembles skipped: beta=",beta,", m0=",m0)
            # #################################
            # ############ beta = 7.2, m = -0.794, mpi/mpiinf<1.3 #####################
            # if beta == 7.2:
            #     if m0 == -0.794:
            #         N_Ls = []
            #         E_pis_t = []
            #         E_pi_errs_t = []
            #         E_pipis_t = []
            #         E_pipi_errs_t = []
            #         for i in range(len(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)])):
            #             if N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i] > 13:
            #                 N_Ls.append(N_L_arr[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                 E_pis_t.append(E_pis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                 E_pi_errs_t.append(E_pi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                 E_pipis_t.append(E_pipis[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #                 E_pipi_errs_t.append(E_pipi_errs[beta_arr.index(beta)][m_arr[beta_arr.index(beta)].index(m0)][i])
            #         if len(N_Ls) > 1:
            #             print("Scattering analysis: beta=",beta,", m0=",m0,", resampling .... ")
            #             res, res_sam = result_sampled(beta,m0,N_Ls,E_pis_t,E_pi_errs_t,E_pipis_t,E_pipi_errs_t)
            #             fn_tmp = "scattering_b72_m794_mpir13_b%1.3f_m%1.3f"%(float(res["beta"]),float(res["m_1"]))
            #             save_to_hdf(res, res_sam, fn_tmp)
            #             res, res_sam = read_from_hdf(fn_tmp)
            # print("Ensembles skipped: beta=",beta,", m0=",m0)
            # #################################