import numpy as np
import matplotlib.pyplot as plt

pdf = False

if __name__ == "__main__":
    plt.figure(figsize=[10,6],dpi=300)

    plt.rcParams["font.size"] = 14

    mass_Sp4 = 100
    mass_chiPT = 100

    sigma_v_data = np.transpose(np.genfromtxt("output/tables/sigma_v_data.csv",delimiter=","))
    DM_halo_data_error = np.transpose(np.genfromtxt("input/sigma_v_data_errors.csv",delimiter=","))

    varr = sigma_v_data[0]
    chiPT_data_in = sigma_v_data[1:3]
    LS_data_in = sigma_v_data[4]
    Sp4_data_in = sigma_v_data[4:]
    LS_data = []
    LS_data.append(0)
    for i in range(len(LS_data_in)):
        LS_data.append(LS_data_in[i])

    varr = np.insert(varr,0,10)

    chiPT_data = []
    Sp4_data = []

    for i in range(len(chiPT_data_in)):
        chiPT_data.append(np.insert(chiPT_data_in[i],0,10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(chiPT_data_in[i][0]))))
    for i in range(len(Sp4_data_in)):
        Sp4_data.append(np.insert(Sp4_data_in[i],0,10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(Sp4_data_in[i][0]))))

    varr[0] = 10

    for i in range(len(chiPT_data)):
        chiPT_data[i][0] = 10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(chiPT_data[i][1]))
    for i in range(len(Sp4_data)):
        Sp4_data[i][0] = 10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(Sp4_data[i][1]))
    LS_data[0]=10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(LS_data[1]))


    for i in range(len(Sp4_data)):
        plt.plot(varr, Sp4_data[i],color="orange")
    plt.errorbar(x=[-1,],y=[-1],xerr=[[1,],[1,]],yerr=[[1,],[1,]],marker="o",color="grey", label ="PRL 116.041302$^{[4]}$",ls ="None", capsize = 3)
    plt.errorbar(x=DM_halo_data_error[0],xerr=[abs(DM_halo_data_error[6]-DM_halo_data_error[0]),abs(DM_halo_data_error[8]-DM_halo_data_error[0])],y=DM_halo_data_error[1],yerr=[abs(DM_halo_data_error[5]-DM_halo_data_error[1]),abs(DM_halo_data_error[3]-DM_halo_data_error[1])],marker="*",markersize=15,color="firebrick",ls ="None", capsize = 3)
    plt.errorbar(x=DM_halo_data_error[10],xerr=[abs(DM_halo_data_error[16]-DM_halo_data_error[10]),abs(DM_halo_data_error[18]-DM_halo_data_error[10])],y=DM_halo_data_error[11],yerr=[abs(DM_halo_data_error[15]-DM_halo_data_error[11]),abs(DM_halo_data_error[13]-DM_halo_data_error[11])],marker="s",markersize=10,color="dodgerblue",ls ="None", capsize = 3)
    plt.errorbar(x=DM_halo_data_error[20],xerr=[abs(DM_halo_data_error[26]-DM_halo_data_error[20]),abs(DM_halo_data_error[28]-DM_halo_data_error[20])],y=DM_halo_data_error[21],yerr=[abs(DM_halo_data_error[25]-DM_halo_data_error[21]),abs(DM_halo_data_error[23]-DM_halo_data_error[21])],marker="p",markersize=10,color="forestgreen",ls ="None", capsize = 3)
    # plt.plot(DM_halo_data_error[30],DM_halo_data_error[31],c="purple", label="ERE Fit ($m_{DM}$=%1.1f$\,$GeV)"%(16.7))                # from old version
    plt.plot(varr,LS_data,c="purple", label="ERE Fit ($m_{DM}$=%1.1f$\,$GeV)"%(16.7))

    for i in range(-5,5):
        plt.plot(varr, varr*10**(i), color="grey",ls = "--", alpha = 0.5)         # lines of constant cross-section

    plt.fill_between(varr,np.min(Sp4_data,axis=0),np.max(Sp4_data,axis=0),color="orange", alpha = 0.7, label = "This work ($m_{DM}$=%i$\,$MeV)"%mass_Sp4)

    if pdf:
        plt.fill_between(varr, chiPT_data[0], chiPT_data[1], color = "green", alpha = 0.3, label = "LO EFT ($m_{DM}$=%i$\,$MeV)"%mass_chiPT)
        plt.fill_between(varr, chiPT_data[0], chiPT_data[1], facecolor="none", edgecolor = "green", hatch="\\\\")
    else:
        plt.fill_between(varr, chiPT_data[0], chiPT_data[1], color = "green", hatch="\\\\", alpha = 0.3, label = "LO EFT ($m_{DM}$=%i$\,$MeV)"%mass_chiPT)


    plt.xlim([2e1,2500])
    plt.ylim([5,4e3])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("$\\left<v\\right>$ in km/s")
    plt.ylabel("$\\left< \sigma v \\right>/m_{DM}$ in $cm^2$/g km/s")
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1,2,0,3]

    legend = plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right')
    legend.get_frame().set_alpha(None)

    plt.tight_layout()
    plt.grid()
    if pdf:
        plt.savefig("output/plots/sigma_v.pdf",bbox_inches = "tight")
    else:
        plt.savefig("output/plots/sigma_v.png",bbox_inches = "tight")
    # plt.show()
    plt.clf()