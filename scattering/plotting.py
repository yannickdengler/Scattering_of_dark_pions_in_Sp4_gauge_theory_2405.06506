"""
    @Author: Yannick Dengler
    @Date:   2024-Mayy-21
    
    This file creats plots for the phase shift from the results calculated with "scattering.py"
    Execute with: "python3 scattering/plotting.py"
 """

import matplotlib.pyplot as plt
import numpy as np
import scattering as result
import math
import os
import matplotlib
from scipy.optimize import bisect

import generalizedzeta as gz
import h5py
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# nth = 50                                           # the how manyth point should be taken for plotting the lines?
def nth(num):
    if num <= 10:
        return 1
    elif num <= 500:
        return num//10
    else:
        return num//100
num_perc = math.erf(1/np.sqrt(2))
font = {'size'   : 16}
matplotlib.rc('font', **font)

def error_of_array(array):
    tempo = []
    for i in range(len(array[0])):
        tmp = array[:,i]
        tmp.sort()
        tempo.append(tmp)
    num_perc = math.erf(1/np.sqrt(2))
    length = len(array)
    result = np.transpose(tempo)
    return result[length//2], abs(result[length//2]-result[math.floor(length*(1-num_perc)/2)]), abs(result[length//2]-result[math.ceil(length*(1+num_perc)/2)]),result[math.floor(length*(1-num_perc)/2)],result[math.ceil(length*(1+num_perc)/2)]


def error_of_1Darray(array):
    array.sort()
    length = len(array)
    return array[length//2], abs(array[length//2]-array[math.floor(length*(1-num_perc)/2)]), abs(array[length//2]-array[math.ceil(length*(1+num_perc)/2)]),array[math.floor(length*(1-num_perc)/2)],array[math.ceil(length*(1+num_perc)/2)]

def marker_beta(beta):
    if beta == 6.9:
        return "^"
    if beta == 7.05:
        return "s"
    if beta == 7.2:
        return "o"
        
def color_beta(beta):
    if beta == 6.9:
        return "#4363d8"
    if beta == 7.05:
        return "#800000"
    if beta == 7.2:
        return "#f032e6"

def E_of_k(k):
    return 2*np.sqrt(1+k**2)

def k_of_Linv(L_inv, q2):
    return 2*np.pi*L_inv*np.sqrt(q2)

def write_result_file(file):
    res,  res_sample = result.read_from_hdf(file)
    beta = res["beta"]
    m12 = res["m_1"]
    m_pi_inf = np.transpose(res_sample["m_pi_inf"])
    a2 = np.transpose(res_sample["a2"])
    b2 = np.transpose(res_sample["b2"])
    a0 = []
    re0 = []
    for i in range(len(a2[0])):
        a0.append(-1/a2[0][i])
        re0.append(2*b2[0][i])
    m_pi_inf_err = error_of_1Darray(m_pi_inf[0])
    a0_err = error_of_1Darray(a0)
    re0_err = error_of_1Darray(re0)
    with open("output/tables/m_pi_a_r.csv","a") as f:
        f.write("%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n"%(beta,m12,m_pi_inf_err[0],m_pi_inf_err[1],m_pi_inf_err[2],a0_err[0],a0_err[1],a0_err[2],re0_err[0],re0_err[1],re0_err[2]))

def plot_curved_errorbars(ax,xarr,yarr,length=0.2,color="green",label="",ratio=1, offset = 2):    # avarage over offset
    '''
    Ratio consists of scale of the axis times the ratio of the figsize. ratio = dx/dy*ax.bbox.height/ax.bbox.width (for axis)
    '''
    def end_cap(x):
        return -x*(xarr[-1]-xarr[-offset])/(yarr[-1]-yarr[-offset])
    def begin_cap(x):
        return -x*(xarr[offset]-xarr[0])/(yarr[offset]-yarr[0])
    ax.plot(xarr,yarr,color=color,label=label)
    xend = xarr[-1]
    len_x = ratio**2
    len_y = 1
    norm_end   = np.sqrt((end_cap(-len_y)  -  end_cap(1))**2+(2*len_x/ratio)**2)/length
    norm_begin = np.sqrt((begin_cap(-len_y)-begin_cap(1))**2+(2*len_x/ratio)**2)/length

    x1,x2,y1,y2 = [xend-len_x/norm_end, xend+len_x/norm_end, yarr[len(yarr)-1]+end_cap(-len_y)/norm_end, yarr[len(yarr)-1]+end_cap(len_y)/norm_end]
    ax.plot([x1,x2],[y1,y2],color=color)
    x1,x2,y1,y2 = [xarr[0]-len_x/norm_begin, xarr[0]+len_x/norm_begin, yarr[0]+begin_cap(-len_y)/norm_begin, yarr[0]+begin_cap(len_y)/norm_begin]
    ax.plot([x1,x2],[y1,y2],color=color)

def get_Zeta_zeros(num_zero = 6):
    zeros = []
    for i in range(num_zero):
        if i > 7:
            i+=1
        zeros.append(bisect(gz.Zeta, i+1e-9,i+1-1e-9))
    return zeros

def plot_m_inf_with_luscher(file, show = True, save = True):
    zeta_zeros = get_Zeta_zeros()               # q2
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    plt.rcParams['figure.figsize'] = [10, 6] 

    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0.1)   
    w_inch, h_inch = ax2.figure.get_size_inches() * ax2.get_position().size
    axins1 = ax2.inset_axes([1.5/w_inch,0.7/h_inch,3/w_inch,1.3/h_inch]) #
    axins1.grid()
    axins1.set_xlim(0.053, 0.087)
    axins1.set_ylim(0.999, 1.0075)


    axins1.set_xticks([0.06,0.08,], labels = [])
    axins1.set_yticks([1.0,], labels = [])

    res,  res_sample = result.read_from_hdf(file)

    E_pi = error_of_array(res_sample["E_pi_prime"])
    E_pipi = error_of_array(res_sample["E_pipi_prime"])
    m_pi_inf = error_of_array(res_sample["m_pi_inf"])
    mass = m_pi_inf[0][0]
    A_R = error_of_array(res_sample["A_R"])
    N_L = res["N_Ls"]
    N_L_inv = [1./L for L in N_L]
    beta = res["beta"]
    m12 = res["m_1"]

    l1=ax2.errorbar(x=N_L_inv, y=E_pi[0], yerr=[E_pi[1],E_pi[2]], ls = "", capsize=3, color = "blue", marker = "v", markersize = 4, label = "$\pi$")
    l2=ax1.errorbar(x=N_L_inv, y=E_pipi[0], yerr=[E_pipi[1],E_pipi[2]], ls = "", capsize=3, color = "red", marker = "^", markersize = 4, label = "$\pi\pi$")

    xarr = np.linspace(0,0.4,200)
    yarr, yarr_m, yarr_p = [[],[],[]]
    for x in xarr:
        yarr.append(result.inf_mass_fit_Goldstone(1./(x),m_pi_inf[0][0],A_R[0][0])/m_pi_inf[0][0])
        yarr_m.append(result.inf_mass_fit_Goldstone(1./(x),m_pi_inf[3][0],A_R[0][0])/m_pi_inf[0][0])
        yarr_p.append(result.inf_mass_fit_Goldstone(1./(x),m_pi_inf[4][0],A_R[0][0])/m_pi_inf[0][0])
        
    ax2.plot(xarr, yarr, color = "blue", label = "Fit")
    ax2.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")

    axins1.plot(xarr, yarr, color = "blue", label = "Fit")
    axins1.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")    
    axins1.errorbar(x=N_L_inv, y=E_pi[0], yerr=[E_pi[1],E_pi[2]], ls = "", capsize=3, color = "blue", marker = "v", markersize = 4, label = "$\pi$")


    L_inv_lat = np.linspace(0,1.2/min(N_L),200)                                              # in lattice units

    E_triv, E_nontriv, E_nontriv_min, E_nontriv_max = [[],[],[],[]]
    for i in range(6):
        E_triv.append([])
        E_nontriv.append([])
        for j in range(len(L_inv_lat)):
            E_triv[i].append(E_of_k(k_of_Linv(L_inv_lat[j], i)/mass))
            E_nontriv[i].append(E_of_k(k_of_Linv(L_inv_lat[j], zeta_zeros[i])/mass))

    for i in range(5):
        if i == 0:
            l3 = ax1.plot(L_inv_lat, E_triv[i+1], color = "black", ls = "dashed", label = "trivial")
            l4 = ax1.plot(L_inv_lat, E_nontriv[i], color = "black", linewidth = .1, label = "$\cot(\delta_0)=0$")
        else:
            ax1.plot(L_inv_lat, E_triv[i+1], color = "black", ls = "dashed")
            ax1.plot(L_inv_lat, E_nontriv[i], color = "black", linewidth = .1)
        ax1.fill_between(L_inv_lat, y1=E_triv[i], y2=E_nontriv[i],color="lightblue",alpha = 0.5)#,hatch="x")
        ax1.fill_between(L_inv_lat, y1=E_nontriv[i], y2=E_triv[i+1],color="lightblue",alpha = 0.7)#,hatch="\\\\")
    ax1.text(x=0.061, y = 2.308, s="$\delta_0 > 0$")
    ax1.text(x=0.081, y = 2.308, s="$\delta_0 < 0$")

    ax1.annotate("", xy=(0.0565, E_of_k(k_of_Linv(0.0565, zeta_zeros[0])/mass)), xytext=(0.065, 2.15),
            arrowprops=dict(facecolor='black',width = 0.5, headwidth = 3, headlength = 3),fontsize = "small")
    ax1.text(s="$\cot(\delta_0)=0$",x = 0.0655, y = 2.135)

    ax1.text(s="$q^2 \in \mathbb{Z}$",x = 0.0555, y = 2.25)
    ax1.annotate("", xy=(0.0454, E_of_k(k_of_Linv(0.0455, 1)/mass)), xytext=(0.055, 2.26),
            arrowprops=dict(facecolor='black',width = 0.5, headwidth = 3, headlength = 3),fontsize = "small")
    ax1.annotate("", xy=(0.0360, E_of_k(k_of_Linv(0.0361, 2)/mass)), xytext=(0.055, 2.26),
            arrowprops=dict(facecolor='black',width = 0.5, headwidth = 3, headlength = 3),fontsize = "small")

    ax2.set_xlabel("a/L")
    ax2.text(s="$E/m_\pi^\infty$", rotation = "vertical", x=-0.012, y = 1.13, fontsize = 14)

    x1,x2,y1,y2 = [0,0.13,2,2.35]
    ax1.set_xlim([x1,x2])
    ax1.set_ylim([y1,y2])
    ax2.set_ylim([0.995,1.13])

    mark_inset(ax2, axins1, loc1=1, loc2=3, fc="none", ec="0.5", zorder = -5)
    ax2.grid()
    ax1.grid()
    ax2.legend([l1, l2],["$\pi$", "$\pi\pi$"], loc = "upper left")
    if save:
        plt.savefig("output/plots/EpiEpipi_%1.3f_m%1.3f.pdf"%(res["beta"],res["m_1"]), bbox_inches="tight")
    if show:
        plt.show()

def plot_ERT_plus_sigma(file, show=False, save = True, rek_lim = True, vesc_lim = True):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0.1)   
    res,  res_sample = result.read_from_hdf(file)
    num_gaussian = len(res_sample["P2_pipi_prime"])

    x1,x2,y1,y2 = [4,5.4,-3.7,-0.6]
    ax1.set_xlim([x1,x2])
    ax1.set_ylim([y1,y2])
    ratio = (x2-x1)/(y2-y1)*ax1.bbox.height/ax1.bbox.width
    ax2.set_xlim([x1,x2])
    ax2.set_ylim([0,13])

    P_cot_PS_pipi_prime = np.transpose(res_sample["P_cot_PS_pipi_prime"])
    s_prime = np.transpose(res_sample["s_pipi_prime"])
    sigma = np.transpose(res_sample["sigma_pipi_prime"])

    inter_UTE = np.transpose(res_sample["UTE_inter_P_cot_PS"])
    P2_arr = np.logspace(np.log10(1e-4), np.log10(3), len(inter_UTE))
    s_arr_of_P2 = []
    for P2 in P2_arr:
        s_arr_of_P2.append(4+4*P2)
    inter_UTE_plot, inter_UTE_plot_m, inter_UTE_plot_p = [[],[],[]]
    inter_sigma = np.transpose(res_sample["sigma_inter_sigma"])
    s_arr = np.logspace(np.log10(4), np.log10(15), len(inter_sigma))


    for i in range(len(inter_UTE)):
        tmp = inter_UTE[i]
        tmp.sort()
        length = len(tmp)
        inter_UTE_plot.append(tmp[length//2])
        inter_UTE_plot_m.append(tmp[math.floor(length*(1-num_perc)/2)])
        inter_UTE_plot_p.append(tmp[math.ceil(length*(1+num_perc)/2)])
    inter_sigma_plot, inter_sigma_plot_m, inter_sigmaplot_p = [[],[],[]]
    for i in range(len(inter_sigma)):
        tmp = inter_sigma[i]
        tmp.sort()
        length = len(tmp)
        inter_sigma_plot.append(tmp[length//2])
        inter_sigma_plot_m.append(tmp[math.floor(length*(1-num_perc)/2)])
        inter_sigmaplot_p.append(tmp[math.ceil(length*(1+num_perc)/2)])

    for i in range(len(s_prime)):
        length = len(s_prime[i])
        sorted_indices = np.argsort(s_prime[i])  
        av = np.arange(math.floor(length*(1-num_perc)/2),math.ceil(length*(1+num_perc)/2),nth(length))
        s_arr2 = [np.mean(s_prime[i][sorted_indices][av[x]:av[x+1]]) for x in range(len(av)-1)]
        P_cot_arr = [np.mean(P_cot_PS_pipi_prime[i][sorted_indices][av[x]:av[x+1]]) for x in range(len(av)-1)]
        if i == 0:
            plot_curved_errorbars(ax1,s_arr2, P_cot_arr, ratio=ratio, color = "green", label = "pipi", offset=5)
        else:            
            plot_curved_errorbars(ax1,s_arr2, P_cot_arr, ratio=ratio, color = "green", label = "pipi", offset=5)
    ax1.plot(s_arr_of_P2, inter_UTE_plot, label = "Fit")
    ax1.fill_between(s_arr_of_P2, inter_UTE_plot_m, inter_UTE_plot_p, alpha = 0.2)
    
    s_pipi_prime_err = error_of_array(np.transpose(s_prime))
    sigma_prime_err = error_of_array(res_sample["sigma_pipi_prime"])
    ax2.errorbar(x=s_pipi_prime_err[0], xerr=[s_pipi_prime_err[1],s_pipi_prime_err[2]], y=sigma_prime_err[0], yerr=[sigma_prime_err[1],sigma_prime_err[2]], color = "green", ls = "", capsize=5, markersize=10, label = "pipi") # , marker = "^"

    ax2.plot(s_arr, inter_sigma_plot, label = "Fit")
    ax2.fill_between(s_arr, inter_sigma_plot_m, inter_sigmaplot_p, alpha = 0.2)

    ax1.grid()
    ax1.set_ylabel("$P \\cot\delta_0/m_\pi^\infty$")
    ax2.grid()
    ax2.set_xlabel("$s/m_\pi^{\infty\,2}$")
    ax2.set_ylabel("$\sigma m_\pi^{\infty\,2}$")
    if save:
        plt.savefig("output/plots/comb_s_b%1.3f_m%1.3f.pdf"%(res["beta"],res["m_1"]), bbox_inches="tight")
    if show:
        plt.show()
    plt.clf()

def plot_a_0_vs_m_f_pi(show=False, save = True):
    plt.figure(figsize=(8,4.8))
    def chipt(x):
        return x/32
    fpi_full = np.genfromtxt("output/tables/fpi_data.csv",delimiter=",")
    xarr = np.linspace(0,10)
    yarr = [chipt(x**2) for x in xarr]
    plt.plot(xarr, yarr, ls = "dashed", color = "green", label = "LO EFT")
    beta_arr = [[6.9,6.9,6.9,6.9],[7.05,7.05],[7.2,7.2]]                      # those with 3 or more datapoints
    m_arr = [[-0.87,-0.9,-0.91,-0.92],[-0.835,-0.85],[-0.78,-0.794]]
    a0_mpi_total_arr = []
    mpifpi_total_arr = []
    out = []
    for i in range(len(beta_arr)):
        for j in range(len(beta_arr[i])):
            a0mpi_arr = []
            rempi_arr = []
            mpi_arr = []
            mpifpi_arr = []
            res,  res_sample = result.read_from_hdf("scattering_b%1.3f_m%1.3f"%(beta_arr[i][j],m_arr[i][j]))
            for k in range(len(res_sample["a2"])):
                a0mpi_arr.append(float(-1/res_sample["a2"][k][0]))
                rempi_arr.append(float(2*res_sample["b2"][k][0]))
                mpi_arr.append(float(res_sample["m_pi_inf"][k][0]))
                a0_mpi_total_arr.append(float(-1/res_sample["a2"][k][0]))
                for l in range(len(fpi_full)):
                    if beta_arr[i][j] == fpi_full[l][0]:
                        if m_arr[i][j] == fpi_full[l][1]:
                            mpifpi_tmp = res_sample["m_pi_inf"][k][0]/np.random.normal(loc=fpi_full[l][2], scale=fpi_full[l][3])
                            mpifpi_arr.append(mpifpi_tmp)
                            mpifpi_total_arr.append(mpifpi_tmp)
            a0mpi_err = error_of_1Darray(a0mpi_arr)
            rempi_err = error_of_1Darray(rempi_arr)
            mpi_err = error_of_1Darray(mpi_arr)
            mpifpi_err = error_of_1Darray(mpifpi_arr)
            out.append([])
            out[len(out)-1].append(beta_arr[i][j])
            out[len(out)-1].append(m_arr[i][j])
            out[len(out)-1].append(a0mpi_err[0])
            out[len(out)-1].append(rempi_err[0])
            out[len(out)-1].append(mpi_err[0])
            out[len(out)-1].append(mpifpi_err[0])
            plt.errorbar(x=[mpifpi_err[0],],xerr=[[mpifpi_err[1],],[mpifpi_err[2],]],y=[a0mpi_err[0],],yerr=[[a0mpi_err[1],],[a0mpi_err[2],]], marker = marker_beta(beta_arr[i][j]), ls = "", capsize=5, markersize=10, color = color_beta(beta_arr[i][j]))
    with open("output/tables/Sp(4)_data.csv", "w") as f:
        for i in range(len(out)):
            f.write("%e,%e,%e,%e,%e,%e\n"%(out[i][0],out[i][1],out[i][2],out[i][3],out[i][4],out[i][5]))
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(6.9), color = color_beta(6.9), label ="$\\beta=6.90$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.05), color = color_beta(7.05), label ="$\\beta=7.05$", s = 60)
    plt.scatter((-10, -9), y = (0,0), marker = marker_beta(7.2), color = color_beta(7.2), label ="$\\beta=7.20$", s = 60)    
    a0_mpi_total_err = error_of_1Darray(a0_mpi_total_arr)
    mpifpi_total_err = error_of_1Darray(mpifpi_total_arr)
    plt.fill_between(x=(-100,100), y1=(a0_mpi_total_err[3],a0_mpi_total_err[3]), y2=(a0_mpi_total_err[4],a0_mpi_total_err[4]), color = "grey", alpha = 0.25)
    plt.axhline(a0_mpi_total_err[0], color = "grey", alpha = 0.5)
    print("a0_mpi_total_err: ", a0_mpi_total_err[0], a0_mpi_total_err[1], a0_mpi_total_err[2], a0_mpi_total_err[3], a0_mpi_total_err[4])
    print("mpifpi_total_err: ", mpifpi_total_err[0], mpifpi_total_err[1], mpifpi_total_err[2], mpifpi_total_err[3], mpifpi_total_err[4])

    plt.xlim([0,6])
    plt.ylim([0,1])
    plt.xlabel("$m_\pi^\infty / f_\\pi$")
    plt.ylabel("$a_0m_\pi^\infty$")
    plt.grid()
    plt.legend()
    if save:
        plt.savefig("output/plots/chipt_comp.pdf", bbox_inches="tight")
    if show:
        plt.show()
    plt.clf()

def write_fpi_file():
    N_L = {}
    fpi = {}
    betas = {}
    m0s = {}
    fpi_err = {}
    with h5py.File("output/hdf5/fitresults.hdf5", "r") as file:
        for key_tmp in file.keys():
            for key in file[key_tmp]["pi"].keys():
                if key[len(key)-13:] == "Delta_fpi_ren":
                    string = str(file[key_tmp]["pi"]["beta"][()])+str(file[key_tmp]["pi"]["m_1"][()])
                    if string in fpi.keys():
                        if file[key_tmp]["pi"]["N_L"][()] > N_L[string]:
                            betas[string] = float(file[key_tmp]["pi"]["beta"][()])
                            m0s[string] = float(file[key_tmp]["pi"]["m_1"][()])
                            N_L[string] = file[key_tmp]["pi"]["N_L"][()]
                            fpi[string] = float(file[key_tmp]["pi"]["fpi_ren"][()])
                            fpi_err[string] = float(file[key_tmp]["pi"]["Delta_fpi_ren"][()])
                    else:
                        betas[string] = float(file[key_tmp]["pi"]["beta"][()])
                        m0s[string] = float(file[key_tmp]["pi"]["m_1"][()])
                        N_L[string] = file[key_tmp]["pi"]["N_L"][()]
                        fpi[string] = float(file[key_tmp]["pi"]["fpi_ren"][()])
                        fpi_err[string] = float(file[key_tmp]["pi"]["Delta_fpi_ren"][()])
    with open("output/tables/fpi_data.csv", "w") as ofile:
        for key in fpi.keys():
            ofile.write("%f,%f,%f,%f\n"%(betas[key], m0s[key], fpi[key], fpi_err[key]))


if __name__ == "__main__":
    beta_arr = [6.9,6.9,6.9,6.9,7.05,7.05,7.2,7.2]
    m_arr = [-0.87,-0.9,-0.91,-0.92,-0.835,-0.85,-0.78,-0.794] 

    # create directory for plots if it doesn#t exit already
    os.makedirs("output/plots", exist_ok=True)

    write_fpi_file()
    plot_a_0_vs_m_f_pi(show=False,save=True)
    for i in range(len(beta_arr)):
        plot_m_inf_with_luscher("scattering_b%1.3f_m%1.3f"%(beta_arr[i],m_arr[i]), show=False,save=True)
        plot_ERT_plus_sigma("scattering_b%1.3f_m%1.3f"%(beta_arr[i],m_arr[i]), show=False,save=True)