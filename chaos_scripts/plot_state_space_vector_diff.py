import matplotlib.pyplot as plt                                                                                                                  
import numpy as np                                                                                                                               
import h5py             
import glob


directories = sorted(glob.glob("plt*/neutrinos"))  
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"                                                                                                              

# plotting magnitud of the state space vectors differences

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

plt.plot(time,ssvecdiff) 
plt.yscale('log')  
plt.ylabel(r'$\Delta \vec{r}_{ss}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/state_space_vector_difference.pdf')   
plt.clf() 



