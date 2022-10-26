import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob

n_estratos=100

# reading the names of the files that contains the data

directories=sorted(glob.glob("plt*/neutrinos"))  
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"

# keys to open the data 

keys1=['f00_Re','f11_Re','f22_Re','f00_Rebar','f11_Rebar','f22_Rebar']
keys2=[['f01_Re','f01_Im'],['f02_Re','f02_Im'],['f12_Re','f12_Im'],['f01_Rebar','f01_Imbar'],['f02_Rebar','f02_Imbar'],['f12_Rebar','f12_Imbar']]  
keys3=['N','N','N','Nbar','Nbar','Nbar']

# labels for the plots

labels1=['$\\rho_{ee}$','$\\rho_{\mu\mu} $','$\\rho_{\\tau\\tau}$','$\\bar{\\rho}_{ee} $','$\\bar{\\rho}_{\mu\mu}$','$\\bar{\\rho}_{\\tau\\tau}$']
labels2=['$\left|\\rho_{e\mu}\\right|$','$\left|\\rho_{e\\tau}\\right|$','$\left|\\rho_{\mu\\tau}\\right|$','$\left|\\bar{\\rho}_{e\mu}\\right|$','$\left|\\bar{\\rho}_{e\\tau}\\right|$','$\left|\\bar{\\rho}_{\mu\\tau}\\right|$']

# doing the analysis for the diagonal values of the density matrix

k3=0

for key in keys1:

    # defining the colors for the plots

    color=['purple','green','orange','black']
    co=0

    # looping over three diferent times for plot it

    for i in [0,int(len(directories)/6),len(directories)-1]:

        #opening the hdf5 files hf is the perturbed data and hg is the original data

        hf=h5py.File(directories[i]+'.h5', 'r')
        hg=h5py.File('../ori/'+directories[i]+'.h5', 'r')

        # saving the time values of each plot
        
        time=np.array(hf.get('time'))[0]

        # getting the density matrix components of all particles

        rho_ii_ori=np.array(hf.get(key))
        rho_ii_per=np.array(hg.get(key))

        #getting the number of neutrinos or antineutrinos each particle carries

        N_ori=np.array(hf.get(keys3[k3]))
        N_per=np.array(hg.get(keys3[k3]))
        
        # creating an array with the n_estratos, that is n_estratos division between 0 and 1

        interval=np.linspace(0,1,n_estratos)

        # arrays to save the number of particles between each of the division last created

        number_neutrino_in_estratos_ori=[]
        number_neutrino_in_estratos_per=[]

        # looping over the estratos

        for i in range(len(interval)-1):
            
            estrato_ini=interval[i]
            estrato_fin=interval[i+1]
            
            number_neutrino_ori=0
            number_neutrino_per=0
            
            for j in range(len(rho_ii_ori)):

                if estrato_ini<rho_ii_ori[j]<=estrato_fin:
                    number_neutrino_ori=number_neutrino_ori+N_ori[j]

                if estrato_ini<rho_ii_per[j]<=estrato_fin:
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
    plt.xlabel(labels1[k3])  
    plt.legend()
    plt.savefig('plots/time_evolution_distribution_'+key+'.pdf')   
    plt.clf()  

    k3=k3+1

# doing the analysis for the non-diagonal values of the density matrix

k3=0

for key in keys2:

    # defining the colors for the plots

    color=['purple','green','orange','black']
    co=0

    # looping over three diferent times for plot it

    for i in [0,int(len(directories)/6),len(directories)-1]:

        #opening the hdf5 files hf is the perturbed data and hg is the original data

        hf=h5py.File(directories[i]+'.h5', 'r')
        hg=h5py.File('../ori/'+directories[i]+'.h5', 'r')

        # saving the time values of each plot
        
        time=np.array(hf.get('time'))[0]

        # getting the density matrix components of all particles

        rho_ii_ori=np.sqrt(np.array(hf.get(key[0]))**2+np.array(hf.get(key[1]))**2)
        rho_ii_per=np.sqrt(np.array(hg.get(key[0]))**2+np.array(hg.get(key[1]))**2)

        #getting the number of neutrinos or antineutrinos each particle carries

        N_ori=np.array(hf.get(keys3[k3]))
        N_per=np.array(hg.get(keys3[k3]))
        
        # creating an array with the n_estratos, that is n_estratos division between 0 and 1

        interval=np.linspace(0,1,n_estratos)

        # arrays to save the number of particles between each of the division last created

        number_neutrino_in_estratos_ori=[]
        number_neutrino_in_estratos_per=[]

        # looping over the estratos

        for i in range(len(interval)-1):
            
            estrato_ini=interval[i]
            estrato_fin=interval[i+1]
            
            number_neutrino_ori=0
            number_neutrino_per=0
            
            for j in range(len(rho_ii_ori)):

                if estrato_ini<rho_ii_ori[j]<=estrato_fin:
                    number_neutrino_ori=number_neutrino_ori+N_ori[j]

                if estrato_ini<rho_ii_per[j]<=estrato_fin:
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
    plt.xlabel(labels2[k3])  
    plt.legend()
    plt.savefig('plots/time_evolution_distribution_'+key[0]+key[1]+'.pdf')   
    plt.clf()  

    k3=k3+1


