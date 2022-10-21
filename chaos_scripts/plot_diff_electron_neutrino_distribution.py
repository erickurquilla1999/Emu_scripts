import numpy as np
import matplotlib.pyplot as plt
import h5py

# read data

hf=h5py.File('plt02000.h5', 'r')
hg=h5py.File('../ori/plt02000.h5', 'r')

rhoee_ori=np.array(hf.get('f00_Re'))
rhoee_per=np.array(hg.get('f00_Re'))

N_ori=np.array(hf.get('N'))
N_per=np.array(hg.get('N'))

# chosing the number of partition to the intervar from 0 to 1 for the analysis

n_estratos=100
interval=np.linspace(0,1,n_estratos)

number_neutrino_in_estratos_ori=[]
number_neutrino_in_estratos_per=[]

for i in range(len(interval)-1):
    
    estrato_ini=interval[i]
    estrato_fin=interval[i+1]
    
    number_neutrino_ori=0
    number_neutrino_per=0
    
    for j in range(len(rhoee_ori)):

        if estrato_ini<rhoee_ori[j]<=estrato_fin:
            number_neutrino_ori=number_neutrino_ori+N_ori[j]

        if estrato_ini<rhoee_per[j]<=estrato_fin:
            number_neutrino_per=number_neutrino_per+N_per[j]

    number_neutrino_in_estratos_ori.append(number_neutrino_ori)
    number_neutrino_in_estratos_per.append(number_neutrino_per)
    
rho_=[]
for i in range(len(interval)-1):
    rho_.append(interval[i]+(interval[i+1]-interval[i])/2)

plt.plot(rho_,np.array(number_neutrino_in_estratos_per)-np.array(number_neutrino_in_estratos_ori))
#plt.yscale('log')
plt.ylabel(r'$N_{per}-N_{ori}$')
plt.xlabel(r'$\rho_{ee}$')  
plt.savefig('plots/diff_distribution_electron_neutrinos.pdf')   
plt.clf()  





