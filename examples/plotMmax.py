"""
    Simple example code used to compute most likely mass M_max^0 of the most massive cluster as a function
    of redshift.
"""

from EVS import *
import time

import numpy as np
import matplotlib.pyplot as plt

def main():
    # setup computation
    evs = evs_core.EVS()
    evs.cosmo.loadCosmology()
    
    # compute
    res = evs.evs_calculation(NM_MAX=500)
    
    # retreive computation results
    Mmax = res.Mmax
    Mmax_s1_inf = res.s1m
    Mmax_s1_sup = res.s1p
    z = res.z_bins
    
    # plot results
    fig = plt.figure(dpi=300)
    ax = fig.subplots()
    
    ax.set_xlabel("z")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylabel('$M_{max}[M_{\odot}.h^{-1}]$')
    ax.set_yscale("log")
    
    ax.step(np.array([0, *z], dtype=float), [Mmax[0], *(Mmax)], "r", where="mid")
    
    ax.step(np.array([0, *z], dtype=float), [Mmax_s1_inf[0], *(Mmax_s1_inf)], "b:", where="mid")
    ax.step(np.array([0, *z], dtype=float), [Mmax_s1_sup[0], *(Mmax_s1_sup)], "b:", where="mid")
    
    ax.grid()
    ax.legend([r"${M_{max}^0}$", r"${\Delta^{1\sigma} M_{max}^0}$"])
   
if __name__ == '__main__':
    main()