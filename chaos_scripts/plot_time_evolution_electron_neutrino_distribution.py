import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob

directories=sorted(glob.glob("plt*/neutrinos"))  
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"

color=['purple','green','orange','black']
co=0

for i in [0,15,80]:

    hf=h5py.File(directories[i]+'.h5', 'r')
    hg=h5py.File('../ori/'+directories[i]+'.h5', 'r')

    time=np.array(hf.get('time'))[0]

    rhoee_ori=np.array(hf.get('f00_Re'))
    rhoee_per=np.array(hg.get('f00_Re'))

    N_ori=np.array(hf.get('N'))
    N_per=np.array(hg.get('N'))
    
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

    plt.plot(rho_,number_neutrino_in_estratos_ori,label='t='+str("{:.2e}".format(time)),linestyle='solid',color=color[co])
    plt.plot(rho_,number_neutrino_in_estratos_per,color=color[co],linestyle='dotted')
    
    co=co+1

plt.yscale('log')
plt.ylabel('number of neutrinos')
plt.xlabel(r'$\rho_{ee}$')  
plt.legend()
plt.savefig('plots/time_evolution_distribution_electron_neutrino.pdf')   
plt.clf()  


