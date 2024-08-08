import numpy as np
import math as m

def sig_dig(x):
    return -int(m.floor(np.log10(abs(x))))

def round_to_sig(x, y=0):
    if y == 0:
        y = x
    return round(x, sig_dig(y))

def get_str(val, errp, errm=None, sd = None):
    if sd == None:
        sd = sig_dig(min(errp, errm))
    if errm == None:
        errm = errp
    res = ""
    if sd == 0:
        res = "$%i^{+%i}_{-%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)))
    else: 
        res = "$(%i^{+%i}_{-%i})\\times 10^{%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)), -sd)
    # res = "$%i^{+%i}_{-%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)))
    return res

def create_latex_table_I2_paper_resuts(pref = ""):
    data = np.genfromtxt("./output/tables/effective_range_parameters"+pref+".csv", delimiter=",")[1:]
    start_str = "\\begin{table}\n\t\centering\n\t\setlength{\\tabcolsep}{2pt}\n\t\\begin{tabular}{|c|c|c|c|c|}\n\t\t\hline\n\t\t"
    after_first_line_str = " \\\\ \hline \hline\n"
    end_str = "\t\end{tabular}\n\t\caption{xxx}\n\t\label{t:results}\n\end{table}\n"
    first_line = "$\\beta$ & $a m_{0}$ & $a m_\pi^\infty\\times 10^4$ & $a_0 m_\pi$ & $r_0 m_\pi$"
    latex_seperator =  " & "

    with open("output/tables/latex_table"+pref+".txt","w") as f:
        f.write(start_str)
        f.write(first_line)
        f.write(after_first_line_str)
        for line in data:
            line_str = "\t\t"
            line_str += "%g"%line[0] + latex_seperator
            line_str += "%g"%line[1] + latex_seperator  

            sd = 0
            line_str += get_str(line[2]*1e4,line[4]*1e4,line[3]*1e4, sd) + latex_seperator                             # m_pi_inf
            sd = 1
            line_str += "$%1.2f^{+%1.2f}_{-%1.2f}$"%(line[5],line[7],line[6]) + latex_seperator                             # a0
            if line[8] > 10:
                line_str += "$%i^{+%i}_{-%i}$"%(line[8],line[10],line[9])                                              # re0
            else:
                line_str += "$%1.1f^{+%1.1f}_{-%1.1f}$"%(line[8],line[10],line[9])                                              # re0

            f.write(line_str)
            f.write("\\\\")
            f.write("\n")
        f.write("\t\\hline\n")
        f.write(end_str)

if __name__ == "__main__":
    for pref in ["","_Fig5.1","_Fig5.2","_Fig5.3","_Fig6"]:
        create_latex_table_I2_paper_resuts(pref)