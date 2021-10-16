from EVS import *
import time

import numpy as np
import matplotlib.pyplot as plt

PATH = 'D:\\Documents\\Python\\Projects\\EVS Cosmic structures\\Derivatives\\Theta = 0.05'

file = open("D:\\Documents\\Python\\Projects\\EVS Cosmic structures\\calculations.dat", "r")
old_text= file.readlines()
file.close()

file = open("D:\\Documents\\Python\\Projects\\EVS Cosmic structures\\calculations.dat", "w")
file.writelines(old_text)     

def main():
    alpha = 0.01
    params = ["Omega_m", "Omega_b_h2", "h", "sigma_8", "n_s", "w_0", 
              "w_a"]
    
    evs = evs_core.EVS()
    z = []
    
    Mmax_plot = []
    legend = []
    
    NM_MAX_arr = [100]

    for k in range(len(NM_MAX_arr)):
        nm_max = NM_MAX_arr[k] # -(NM_MAX_max - NM_MAX_min)/n_it * (0 - k) + NM_MAX_min
        
        for param in params:
            evs.cosmo.loadCosmology()
            defparameters = evs.cosmo.getCosmology()
            
            theta = defparameters[param]
            dtheta = 5e-2
            
            Mmax = []
            
            for i, epsilon in enumerate([-2, -1, 1, 2]):
                parameters = defparameters.copy()
                parameters[param] = theta + epsilon * dtheta
                
                evs.cosmo.setCosmology(**parameters)
                
                res = evs.evs_calculation(NM_MAX=int(nm_max))
                
                Mmax.append(res.Mmax)
                z = res.z_bins
                            
            Mmax = np.array(Mmax, dtype=np.float64)
                
            dF_dtheta = 1/(12 * dtheta) * (Mmax[0] - 8 * Mmax[1] + 8 * Mmax[2] - Mmax[3])
            
            Mmax_plot.append(abs(dF_dtheta))
            legend.append(str(param))
            
            print("results for {0} (NM_MAX={1}) : ".format(param, nm_max), file=file, end='  ')
            np.savetxt(file, dF_dtheta, newline=' ')
            print("", file=file)
        
    
    print(legend)
    print(Mmax_plot)
    
    file.close()

    fig = plt.figure(dpi=300)

    ax = fig.add_axes([0, 0, 0.5, 0.5])
    ax.set_xlabel("z")
    ax.set_xlim(0.0, 1.0)
    ax.set_title("NM_MAX = " + str(nm_max))
    ax.set_ylabel('$dM_{max}/d \Theta$')
    ax.set_yscale("log")
    
    for j in range(len(Mmax_plot)):
        ax.step(np.array([0, *(z)], dtype=float), [Mmax_plot[j][0], *(Mmax_plot[j])], where="mid")
    
    ax.legend(legend, bbox_to_anchor=(1.7, 1))
    
def print_to_file(dtheta, nm_max, z, data):
    time_local = time.localtime()
    filename = "derivatives_{0}_{1}_{2}_{3}:{4}.evsdat".format(time_local.tm_mday, time_local.tm_mon, time_local.tm_year, time_local.tm_hour, time_local.tm_sec) 
    
    file = open(filename, "w+")
    
    print("type : derivatives", file=file)
    print("\ndate : {0}-{1}-{2} {3}:{4}".format(time_local.tm_mday, time_local.tm_mon, time_local.tm_year, time_local.tm_hour, time_local.tm_sec), file=file)
   
    print("\ndtheta = {0:5.2f}".format(dtheta), file=file)
    print("dnm_max = {0:5f}".format(nm_max), file=file)
    
    print("\ndata :", file=file)
        
    for i in range(len(data)):
        print("dnm_max = {0:5f}".format(nm_max), file=file)

if __name__ == '__main__':
    main()