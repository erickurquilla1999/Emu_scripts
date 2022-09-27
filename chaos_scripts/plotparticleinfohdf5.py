import numpy as np                                                                                                                                                                                            
import matplotlib.pyplot as plt 
import glob
import h5py

directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos" 

keys=['pos_x','pos_y','pos_z','time','x', 'y', 'z', 'pupx', 'pupy', 'pupz', 'pupt', 'N', 'L', 'f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im' ,'f22_Re', 'Nbar' ,'Lbar', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar' ,'f12_Imbar', 'f22_Rebar']

alldata=[[],[],[],[],[],[],[],[],[],[],[],[],[]]

for dir in directories:
     
    print('.................................')
    print(dir)  
    print('.................................') 

    hf = h5py.File(dir+'.h5', 'r')
    
    for i in keys:
        print(len(np.array(hf.get(i))))

    alldata[0].append(hf.get('time')[0])

    for i in range(13,22):
        alldata[i-12].append(np.average(hf.get(keys[i])))

    alldata[10].append(np.average(np.sqrt(np.array(hf.get(keys[14]))**2+np.array(hf.get(keys[15]))**2)))
    alldata[11].append(np.average(np.sqrt(np.array(hf.get(keys[16]))**2+np.array(hf.get(keys[17]))**2))) 
    alldata[12].append(np.average(np.sqrt(np.array(hf.get(keys[19]))**2+np.array(hf.get(keys[20]))**2)))

labels=['$<\\rho_{ee}>$','$Re<\\rho_{e\mu}>$','$Im<\\rho_{e\mu}>$','$Re<\\rho_{e\\tau}>$','$Im<\\rho_{e\\tau}>$','$<\\rho_{\mu\mu}>$','$Re<\\rho_{\mu\\tau}>$','$Im<\\rho_{\mu\\tau}>$','$<\\rho_{\tau\tau}>$','$|<\\rho_{e\mu}>|$','$|<\\rho_{e\\tau}>|$','$|<\\rho_{\mu\\tau}>$|']
names=['h5_rho_ee_av.pdf','h5_Re_rho_eu_av.pdf','h5_Im_rho_eu_av.pdf','h5_Re_rho_et_av.pdf','h5_Im_rho_et_av.pdf','h5_rho_uu_av.pdf','h5_Re_rho_ut_av.pdf','h5_Im_rho_ut_av.pdf','h5_rho_tt_av.pdf','h5_mag_rho_eu_av.pdf','h5_mag_rho_et_av.pdf','h5_mag_rho_mt_av.pdf']

for i in range(1,len(alldata)):

    plt.plot(alldata[0],alldata[i]) 
    plt.ylabel(labels[i-1])
    plt.xlabel(r'Time (s)')  
    plt.savefig(names[i-1])  
    plt.clf() 

    
    
     
    
    
    
