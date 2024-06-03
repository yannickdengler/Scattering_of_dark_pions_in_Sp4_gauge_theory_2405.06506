import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 14

color_arr = ["blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick","blue", "green", "red", "purple", "orange", "olive", "skyblue", "lime", "black", "grey", "fuchsia", "peru", "firebrick",]


a_data_linsig = np.transpose(np.genfromtxt("input/a_data_linsig_90.dat"))
re_data_linsig = np.transpose(np.genfromtxt("input/re_data_linsig_90.dat"))

hbarc_GeV = 0.1973
Sp4_data = np.transpose(np.genfromtxt("output/Sp(4)_data.dat"))
Sp4_masses = Sp4_data[4]
Sp4_a = []
Sp4_re = []
for i in range(len(Sp4_masses)):
    Sp4_a.append(-hbarc_GeV*Sp4_masses[i]/(Sp4_data[2][i]))
    Sp4_re.append(2*hbarc_GeV*Sp4_masses[i]*Sp4_data[3][i])


def plot_a():
    mass_plot = np.logspace(-3,2,500)
    a_Sp4_plot = []
    for i in range(len(Sp4_a)):
        a_Sp4_plot.append([])
        for j in range(len(mass_plot)):
            a_Sp4_plot[i].append(Sp4_a[i]/mass_plot[j])

    plt.fill_between(a_data_linsig[1],11.1103,a_data_linsig[0], label='"Kondo et al."', alpha = 0.7, color = "purple", ls = "None")
    plt.fill_between(mass_plot,np.min(a_Sp4_plot,axis=0),np.max(a_Sp4_plot,axis=0), color = "orange", alpha = 0.7, label = "This work", ls = "None")

    plt.xlim([1e0,1e2])
    plt.ylim([1e-3,1e2])
    plt.xlabel("$m_{DM}$ in GeV")
    plt.ylabel("$a_0$ in fm")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="center right")
    plt.grid()
    plt.savefig("plots/a_vs_mass.pdf",bbox_inches = "tight")
    plt.clf()
    # plt.show()

def plot_re():
    mass_plot = np.logspace(0,2,500)
    re_Sp4_plot = []
    for i in range(len(Sp4_a)):
        re_Sp4_plot.append([])
        for j in range(len(mass_plot)):
            re_Sp4_plot[i].append(Sp4_re[i]/mass_plot[j])

    print(len(re_Sp4_plot))
    print(len(re_Sp4_plot[0]))

    plt.fill_between(re_data_linsig[1],0,re_data_linsig[0], color="purple", label='"Kondo et al."', alpha = 0.7, ls = "None")
    plt.fill_between(mass_plot,np.min(re_Sp4_plot,axis=0),np.max(re_Sp4_plot,axis=0), color = "orange", alpha = 0.7, label = "This work", ls = "None")

    plt.xlim([5,30])
    plt.ylim([-4,7.5])
    plt.xlabel("$m_{DM}$ in GeV")
    plt.ylabel("$r_0$ in fm")
    plt.xscale("log")
    plt.legend(loc="upper left")
    plt.grid()
    plt.savefig("plots/re_vs_mass.pdf",bbox_inches = "tight")
    plt.clf()
    # plt.show()

plot_a()
plot_re()