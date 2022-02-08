# Run from /ocean/projects/phy200048p/shared to generate plot showing time evolution of <fee> at different dimensionalities

import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator,LogLocator)


base=["N","Fx","Fy","Fz"]
diag_flavor=["00","11","22"]
offdiag_flavor=["01","02","12"]
re=["Re","Im"]
# real/imag
R=0
I=1
    

def offdiagMag(f):
    return np.sqrt(f[:,0,1,R]**2 + f[:,0,1,I]**2 +
                   f[:,0,2,R]**2 + f[:,0,2,I]**2 +
                   f[:,1,2,R]**2 + f[:,1,2,I]**2)


######################
# read averaged data #
######################
def plotdata(filename, t_in):
    
    avgData = h5py.File(filename,"r")
    
    t=np.array(avgData["t"])
    k=np.array(avgData["k"])
    Nee=np.array(avgData["N00_FFT"])
    Nxx=np.array(avgData["N11_FFT"])
    Nex=np.array(avgData["N01_FFT"])
    avgData.close()

    # get time closest to t
    dt = np.abs(t-t_in)
    it = np.argmin(dt)
    trace = Nee[it,0]+Nex[it,0]
    print(it,t[it])
    return k, (Nex/trace)[it, :-1]

################
# plot options #
################
mpl.rcParams['font.size'] = 22
mpl.rcParams['font.family'] = 'serif'
mpl.rc('text', usetex=True)
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2


fig, ax = plt.subplots(1,1, figsize=(6,5))

##############
# formatting #
##############
ax.tick_params(axis='both', which='both', direction='in', right=True,top=True)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.minorticks_on()
#ax.grid(which='both')

#############
# plot data #
#############
tplot = 0.2e-9
filename_emu_2f = "/global/project/projectdirs/m3761/Evan/Fiducial_3D_2F/reduced_data_fft_power.h5"
filename_emu_3f = "/global/project/projectdirs/m3761/Evan/Fiducial_3D_3F_reduced_data_fft_power.h5"
filename_bang = "/global/project/projectdirs/m3761/FLASH/FFI_3D/fid/sim1/reduced_data_fft_power_nov4_test_hdf5_chk.h5"
k1,N1 = plotdata(filename_emu_2f,tplot)
ax.semilogy(k1, N1, 'k-', label=r'${\rm emu\,\,(2f)}$')
k2,N2 = plotdata(filename_emu_3f,tplot)
ax.semilogy(k2, N2, 'k--', label=r'${\rm emu\,\,(3f)}$')
k3,N3 = plotdata(filename_bang,tplot)
ax.semilogy(k3, N3, 'r-', label=r'${\rm FLASH\,\,(2f)}$')
ax.set_xlabel(r"$k,({\rm cm}^{-1})$")
ax.set_ylabel(r"$\widetilde{N}_{ex}$")
ax.legend(loc='upper right')
plt.savefig("N_ex_FFT.pdf", bbox_inches="tight")
print(N1[0]/N3[0])
print(N3[0]/N1[0])
