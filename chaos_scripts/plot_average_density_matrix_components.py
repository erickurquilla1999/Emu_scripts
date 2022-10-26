import numpy as np                                                                                                                                                                                            
import matplotlib.pyplot as plt 
import glob
import h5py

# keys to read and analysis the data

keys1=['f00_Re','f11_Re','f22_Re','f00_Rebar','f11_Rebar','f22_Rebar']
keys2=[['f01_Re','f01_Im'],['f02_Re','f02_Im'],['f12_Re','f12_Im'],['f01_Rebar','f01_Imbar'],['f02_Rebar','f02_Imbar'],['f12_Rebar','f12_Imbar']]

# reading the original and perturbed data files 

dir_per=sorted(glob.glob("plt*/neutrinos"))
dir_per=[dir_per[i].split('/')[0] for i in range(len(dir_per))] # remove "neutrinos" 

dir_ori=sorted(glob.glob("../ori/plt*/neutrinos"))
dir_ori=[dir_ori[i].split('/')[2] for i in range(len(dir_ori))] # remove "neutrinos" 

# Creating the arrays that will save the data
# 0:time 1:rho_ee 2:rho_uu 3:rho_tt 4:rhobar_ee 5:rhobar_uu 6:rhobar_tt
# 7:rho_eu 8:rho_et 9:rho_ut 10:rhobar_eu 11:rhobar_et 12:rhobar_ut

alldata_per=[[],[],[],[],[],[],[],[],[],[],[],[],[]] 
alldata_ori=[[],[],[],[],[],[],[],[],[],[],[],[],[]]

# looping over the perturbed files for save the average data in the alldata_per list

for dir in dir_per:
     
    hf = h5py.File(dir+'.h5', 'r')

    alldata_per[0].append(hf.get('time')[0])

    a=1
    for key in keys1:
        alldata_per[a].append(np.average(hf.get(key)))
        a=a+1

    b=7
    for key in keys2:
        alldata_per[b].append(np.average(np.sqrt(np.array(hf.get(key[0]))**2+np.array(hf.get(key[1]))**2)))
        b=b+1
    
    hf.close()

# looping over the original files for save the average data in the alldata_ori list

for dir in dir_ori:
     
    hf = h5py.File('../ori/'+dir+'.h5', 'r')

    alldata_ori[0].append(hf.get('time')[0])

    a=1
    for key in keys1:
        alldata_ori[a].append(np.average(hf.get(key)))
        a=a+1

    b=7
    for key in keys2:
        alldata_ori[b].append(np.average(np.sqrt(np.array(hf.get(key[0]))**2+np.array(hf.get(key[1]))**2)))
        b=b+1
    
    hf.close()
    
# defining the labels for the plots

labels=['$\left<\\rho_{ee}\\right>$','$\left<\\rho_{\mu\mu}\\right>$','$\left<\\rho_{\\tau\\tau}\\right>$','$\left<\\bar{\\rho}_{ee}\\right>$','$\left<\\bar{\\rho}_{\mu\mu}\\right>$','$\left<\\bar{\\rho}_{\\tau\\tau}\\right>$','$\left<\left|\\rho_{e\mu}\\right|\\right>$','$\left<\left|\\rho_{e\\tau}\\right|\\right>$','$\left<\left|\\rho_{\mu\\tau}\\right|\\right>$','$\left<\left|\\bar{\\rho}_{e\mu}\\right|\\right>$','$\left<\left|\\bar{\\rho}_{e\\tau}\\right|\\right>$','$\left<\left|\\bar{\\rho}_{\mu\\tau}\\right|\\right>$']

# saving the names of the pdf files that will be saved

names=['average_rho_ee.pdf','average_rho_uu.pdf','average_rho_tt.pdf','average_rhobar_ee.pdf','average_rhobar_uu.pdf','average_rhobar_tt.pdf','average_rho_eu.pdf','average_rho_et.pdf','average_rho_ut.pdf','average_rhobar_eu.pdf','average_rhobar_et.pdf','average_rhobar_ut.pdf']

# looping over the original and perturbed data for plot them together

for i in range(1,len(alldata_per)):
    
    plt.plot(alldata_ori[0],alldata_ori[i],label='original') 
    plt.plot(alldata_per[0],alldata_per[i],label='perturbed')
    plt.ylabel(labels[i-1])
    plt.xlabel(r'Time (s)') 
    plt.legend()
    plt.savefig('plots/'+names[i-1])  
    plt.clf() 

# looping over the neutrino original data for plot each components in same plot

for i in [1,2,3,7,8,9]:
    plt.plot(alldata_ori[0],np.array(alldata_ori[i]),label=labels[i-1])

plt.ylabel(r'$ \left < \left | \rho_{ii} \right | \right >$')
plt.xlabel(r'Time (s)')
#plt.yscale('log')
plt.legend()
plt.savefig('../ori/plots/average_rho_all.pdf')
plt.clf()

# looping over the antineutrino original data for plot each components in same plot

for i in [4,5,6,10,11,12]:
    plt.plot(alldata_ori[0],np.array(alldata_ori[i]),label=labels[i-1])

plt.ylabel(r'$ \left < \left | \bar{\rho}_{ii} \right | \right >$')
plt.xlabel(r'Time (s)')
#plt.yscale('log')
plt.legend()
plt.savefig('../ori/plots/average_rhobar_all.pdf')
plt.clf()

# looping over the neutrino perturbed data for plot each components in same plot

for i in [1,2,3,7,8,9]:
    plt.plot(alldata_per[0],np.array(alldata_per[i]),label=labels[i-1])

plt.ylabel(r'$ \left < \left | \rho_{ii} \right | \right >$')
plt.xlabel(r'Time (s)')
#plt.yscale('log')
plt.legend()
plt.savefig('plots/average_rho_all.pdf')
plt.clf()

# looping over the antineutrino original data for plot each components in same plot

for i in [4,5,6,10,11,12]:
    plt.plot(alldata_per[0],np.array(alldata_per[i]),label=labels[i-1])

plt.ylabel(r'$ \left < \left | \bar{\rho}_{ii} \right | \right >$')
plt.xlabel(r'Time (s)')
#plt.yscale('log')
plt.legend()
plt.savefig('plots/average_rhobar_all.pdf')
plt.clf()

