import matplotlib.pyplot as plt                                                                                                                  
import numpy as np                                                                                                                               
import h5py             
import glob

directories = sorted(glob.glob("plt*/neutrinos"))        
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"  

#labels=['pos_x','pos_y','pos_z', 'time', 'x', 'y', 'z', 'pupx', 'pupy', 'pupz', 'pupt', 'N', 'L', 'f00_Re', 'f01_Re', 'f01_Im', 'f11_Re', 'Nbar', 'Lbar', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f11_Rebar']
#labels=['f00_Re', 'f01_Re', 'f01_Im', 'f11_Re','f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f11_Rebar']
labels=['f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im', 'f11_Re','f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar', 'f12_Imbar', 'f11_Rebar']

time=[]
ssvecdiff=[]

for dir in directories[1:]:

    hper = h5py.File(dir+'.h5', 'r') 
    hori = h5py.File('../ori/'+dir+'.h5', 'r')
    
    time.append(hori.get('time')[0])

    ssdiff=0

    for label in labels:
        ssdiff=ssdiff+np.sum((np.array(hper.get(label))-np.array(hori.get(label)))**2)

    ssvecdiff.append(np.sqrt(ssdiff))

    hper.close()
    hori.close()

time=np.array(time)

initialvalue=3
lastvalue=10

lyapunov=(1/(time[lastvalue]-time[initialvalue]))*np.log(ssvecdiff[lastvalue]/ssvecdiff[initialvalue])

fitdata=ssvecdiff[initialvalue]*np.exp(lyapunov*(time[initialvalue:lastvalue+1]-time[initialvalue]))

plt.plot(time,ssvecdiff) 
plt.plot(time[initialvalue:lastvalue+1],fitdata,linestyle='dashed',label='$\lambda=$'+str(lyapunov))
plt.legend()
plt.yscale('log')
plt.ylabel(r'$\Delta \vec{r}_{ss}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/fit_state_space_vec_diff.pdf')   
plt.clf()


