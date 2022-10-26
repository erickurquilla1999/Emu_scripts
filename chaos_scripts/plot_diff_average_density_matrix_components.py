import numpy as np                                                                                                                                                                                            
import matplotlib.pyplot as plt 
import glob
import h5py

# saving the keys to read the data

keys1=['f00_Re','f11_Re','f22_Re','f00_Rebar','f11_Rebar','f22_Rebar']
keys2=[['f01_Re','f01_Im'],['f02_Re','f02_Im'],['f12_Re','f12_Im'],['f01_Rebar','f01_Imbar'],['f02_Rebar','f02_Imbar'],['f12_Rebar','f12_Imbar']]

# getting the names of the files the contain the data

dir_per=sorted(glob.glob("plt*/neutrinos"))
dir_per=[dir_per[i].split('/')[0] for i in range(len(dir_per))] # remove "neutrinos" 

# Creating the arrays that will save the data 
# 0:time 1:rho_ee 2:rho_uu 3:rho_tt 4:rhobar_ee 5:rhobar_uu 6:rhobar_tt
# 7:rho_eu 8:rho_et 9:rho_ut 10:rhobar_eu 11:rhobar_et 12:rhobar_ut

alldata_per=[[],[],[],[],[],[],[],[],[],[],[],[],[]]
alldata_ori=[[],[],[],[],[],[],[],[],[],[],[],[],[]]

# Looping over the perturbed data files names

for dir in dir_per[1:]:
    
    # opening the original and the perturbed hdf5 files

    h5_per=h5py.File(dir+'.h5', 'r')
    h5_ori=h5py.File('../ori/'+dir+'.h5', 'r')
    
    # saving the time value

    alldata_per[0].append(h5_per.get('time')[0])
    alldata_ori[0].append(h5_ori.get('time')[0])

    # looping on key1 (diagonal components of the density matrix) and saving the average values 

    a=1
    for key in keys1:
        alldata_per[a].append(np.average(h5_per.get(key)))
        alldata_ori[a].append(np.average(h5_ori.get(key)))
        a=a+1

    # looping over key2 (nondiagonal components of the density matrix) computing the magnitud and saving the average

    b=7
    for key in keys2:
        alldata_per[b].append(np.average(np.sqrt(np.array(h5_per.get(key[0]))**2+np.array(h5_per.get(key[1]))**2)))
        alldata_ori[b].append(np.average(np.sqrt(np.array(h5_ori.get(key[0]))**2+np.array(h5_ori.get(key[1]))**2)))
        b=b+1
    
    ## closing hdf5 files

    h5_per.close() 
    h5_ori.close()

# labels for the plots

labels=['$\Delta\left<\\rho_{ee}\\right>$','$\Delta\left<\\rho_{\mu\mu}\\right>$','$\Delta\left<\\rho_{\\tau\\tau}\\right>$','$\Delta\left<\\bar{\\rho}_{ee}\\right>$','$\Delta\left<\\bar{\\rho}_{\mu\mu}\\right>$','$\Delta\left<\\bar{\\rho}_{\\tau\\tau}\\right>$','$\Delta\left<\left|\\rho_{e\mu}\\right|\\right>$','$\Delta\left<\left|\\rho_{e\\tau}\\right|\\right>$','$\Delta\left<\left|\\rho_{\mu\\tau}\\right|\\right>$','$\Delta\left<\left|\\bar{\\rho}_{e\mu}\\right|\\right>$','$\Delta\left<\left|\\bar{\\rho}_{e\\tau}\\right|\\right>$','$\Delta\left<\left|\\bar{\\rho}_{\mu\\tau}\\right|\\right>$']

# names for the pdf files that will be saved

names=['average_rho_ee.pdf','average_rho_uu.pdf','average_rho_tt.pdf','average_rhobar_ee.pdf','average_rhobar_uu.pdf','average_rhobar_tt.pdf','average_rho_eu.pdf','average_rho_et.pdf','average_rho_ut.pdf','average_rhobar_eu.pdf','average_rhobar_et.pdf','average_rhobar_ut.pdf']

# looping over the data for plot each components

for i in range(1,len(alldata_ori)):
    plt.plot(alldata_ori[0],np.absolute(np.array(alldata_ori[i])-np.array(alldata_per[i]))) 
    plt.ylabel(labels[i-1])
    plt.xlabel(r'Time (s)') 
    plt.savefig('plots/difference_'+names[i-1])  
    plt.clf() 

# looping ever the neutrino data for plot each components in same plot

for i in [1,2,3,7,8,9]:
    plt.plot(alldata_ori[0],np.absolute(np.array(alldata_ori[i])-np.array(alldata_per[i])),label=labels[i-1]) 
    
plt.ylabel(r'$ \Delta \left < \left | \rho_{ii}  \right | \right > $')
plt.xlabel(r'Time (s)') 
plt.yscale('log')
plt.legend()
plt.savefig('plots/difference_average_rho_all.pdf')  
plt.clf() 

# looping ever the antineutrino data for plot each components in same plot

for i in [4,5,6,10,11,12]:
    plt.plot(alldata_ori[0],np.absolute(np.array(alldata_ori[i])-np.array(alldata_per[i])),label=labels[i-1]) 
    
plt.ylabel(r'$ \Delta \left < \left | \bar{\rho}_{ii} \right | \right > $')
plt.xlabel(r'Time (s)') 
plt.yscale('log')
plt.legend()
plt.savefig('plots/difference_average_rhobar_all.pdf')  
plt.clf() 



