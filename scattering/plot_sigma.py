import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=[10,6])

plt.rcParams["font.size"] = 14

color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]

mass_Sp4 = 100                                               # For rescaling of our DM mass
mass_chiPT = 100                                               # For rescaling of our DM mass

sigma_v_data = np.transpose(np.genfromtxt("output/sigma_v_data.dat"))
DM_halo_data_error = np.transpose(np.genfromtxt("output/sigma_v_data_errors.dat"))

varr = sigma_v_data[0]
chiPT_data = sigma_v_data[1:3]
Sp4_data = sigma_v_data[3:]

varr[0] = 10

for i in range(len(chiPT_data)):
    chiPT_data[i][0] = 10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(chiPT_data[i][1]))
for i in range(len(Sp4_data)):
    Sp4_data[i][0] = 10**(np.log10(varr[0])-np.log10(varr[1])+np.log10(Sp4_data[i][1]))

plt.errorbar(x=[-1,],y=[-1],xerr=[[1,],[1,]],yerr=[[1,],[1,]],marker="o",color="grey", label ="PRL 116.041302$^{[4]}$",ls ="None", capsize = 3)
plt.errorbar(x=DM_halo_data_error[0],xerr=[abs(DM_halo_data_error[6]-DM_halo_data_error[0]),abs(DM_halo_data_error[8]-DM_halo_data_error[0])],y=DM_halo_data_error[1],yerr=[abs(DM_halo_data_error[5]-DM_halo_data_error[1]),abs(DM_halo_data_error[3]-DM_halo_data_error[1])],marker="o",color="firebrick",ls ="None", capsize = 3)
plt.errorbar(x=DM_halo_data_error[10],xerr=[abs(DM_halo_data_error[16]-DM_halo_data_error[10]),abs(DM_halo_data_error[18]-DM_halo_data_error[10])],y=DM_halo_data_error[11],yerr=[abs(DM_halo_data_error[15]-DM_halo_data_error[11]),abs(DM_halo_data_error[13]-DM_halo_data_error[11])],marker="o",color="dodgerblue",ls ="None", capsize = 3)
plt.errorbar(x=DM_halo_data_error[20],xerr=[abs(DM_halo_data_error[26]-DM_halo_data_error[20]),abs(DM_halo_data_error[28]-DM_halo_data_error[20])],y=DM_halo_data_error[21],yerr=[abs(DM_halo_data_error[25]-DM_halo_data_error[21]),abs(DM_halo_data_error[23]-DM_halo_data_error[21])],marker="o",color="forestgreen",ls ="None", capsize = 3)
plt.plot(DM_halo_data_error[30],DM_halo_data_error[31],c="purple", label="ERE Fit ($m_{DM}$=%1.1f$\,$GeV)"%(16.7))

for i in range(-5,5):
    plt.plot(varr, varr*10**(i), color="grey",ls = "--", alpha = 0.5)         # lines of constant cross-section

plt.fill_between(varr,np.min(Sp4_data,axis=0),np.max(Sp4_data,axis=0),color="orange", alpha = 0.7, label = "This work ($m_{DM}$=%i$\,$MeV)"%mass_Sp4)

plt.fill_between(varr, chiPT_data[0], chiPT_data[1], color = "green", alpha = 0.3, label = "LO EFT ($m_{DM}$=%i$\,$MeV)"%mass_chiPT)#, ls = "--")

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
plt.savefig("plots/sigma_v_MJ.pdf",bbox_inches = "tight")
plt.show()
plt.clf()