import numpy as np                                                                                                                                                                                            
import matplotlib.pyplot as plt 
import glob

directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos" 

#labels=['pos_x','pos_y','pos_z','time','x', 'y', 'z', 'pupx', 'pupy', 'pupz', 'pupt', 'N', 'L', 'f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im' ,'f22_Re', 'Nbar' ,'Lbar', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar' ,'f12_Imbar', 'f22_Rebar']

alldata=[[],[],[],[],[],[],[],[],[],[],[],[],[]]

for dir in directories:
    
    print('.................................')
    print(dir)
    print('.................................')
   
    data=np.loadtxt(dir+'.txt', unpack = True, skiprows=1)
    
    for i in data:
        print(len(i))

    alldata[0].append(data[3][0])

    for i in range(13,22):
        alldata[i-12].append(np.average(data[i]))

    alldata[10].append(np.average(np.sqrt(data[14]**2+data[15]**2)))
    alldata[11].append(np.average(np.sqrt(data[16]**2+data[17]**2))) 
    alldata[12].append(np.average(np.sqrt(data[19]**2+data[20]**2)))


labels=['$<\\rho_{ee}>$','$Re<\\rho_{e\mu}>$','$Im<\\rho_{e\mu}>$','$Re<\\rho_{e\\tau}>$','$Im<\\rho_{e\\tau}>$','$<\\rho_{\mu\mu}>$','$Re<\\rho_{\mu\\tau}>$','$Im<\\rho_{\mu\\tau}>$','$<\\rho_{\\tau\\tau}>$','$|<\\rho_{e\mu}>|$','$|<\\rho_{e\\tau}>|$','$|<\\rho_{\mu\\tau}>$|']
names=['txt_rho_ee_av.pdf','txt_Re_rho_eu_av.pdf','txt_Im_rho_eu_av.pdf','txt_Re_rho_et_av.pdf','txt_Im_rho_et_av.pdf','txt_rho_uu_av.pdf','txt_Re_rho_ut_av.pdf','txt_Im_rho_ut_av.pdf','txt_rho_tt_av.pdf','txt_mag_rho_eu_av.pdf','txt_mag_rho_et_av.pdf','txt_mag_rho_mt_av.pdf']

for i in range(1,len(alldata)):

    plt.plot(alldata[0],alldata[i]) 
    plt.ylabel(labels[i-1])
    plt.xlabel(r'Time (s)')  
    plt.savefig(names[i-1])  
    plt.clf() 

    
    
     
    
    
    
