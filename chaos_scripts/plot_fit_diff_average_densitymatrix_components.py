import matplotlib.pyplot as plt                                                                                                                  
import numpy as np                                                                                                                               
import h5py             
import glob

#########################################################################################
# plot the difference between rho_ee original and rho_ee perturbed both average over all particles
#########################################################################################

directories = sorted(glob.glob("plt*/neutrinos"))        
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"  

time=[]
rhoee_per=[]
rhoee_ori=[]

for dir in directories:

    hf = h5py.File(dir+'.h5', 'r')
    rhoee_per.append(np.average(np.array(hf.get('f00_Re'))))
    time.append(hf.get('time')[0])
    hf.close() 

    hf = h5py.File('../ori/'+dir+'.h5', 'r')
    rhoee_ori.append(np.average(np.array(hf.get('f00_Re'))))
    hf.close() 

diff_rho_ee=np.absolute(np.array(rhoee_per)[1:]-np.array(rhoee_ori)[1:])
time=time[1:]

plt.plot(time,diff_rho_ee) 

initialvalue=5
lastvalue=8

lyapunov=(1/(time[lastvalue]-time[initialvalue]))*np.log(diff_rho_ee[lastvalue]/diff_rho_ee[initialvalue])
fitdata=diff_rho_ee[initialvalue]*np.exp(lyapunov*(time[initialvalue:lastvalue+1]-time[initialvalue]))

plt.plot(time[initialvalue:lastvalue+1],fitdata,linestyle='dashed',label='$\lambda=$'+str(lyapunov))

plt.legend()
plt.yscale('log')  
plt.ylabel(r'$\Delta \rho_{ee}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/fit_difference_average_rho_ee.pdf')   
plt.clf()




#########################################################################################
# plot the difference between rho_uu original and rho_uu perturbed both average over all particles
#########################################################################################

directories = sorted(glob.glob("plt*/neutrinos"))        
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"  

time=[]
rhouu_per=[]
rhouu_ori=[]

for dir in directories:

    hf = h5py.File(dir+'.h5', 'r')
    rhouu_per.append(np.average(np.array(hf.get('f11_Re'))))
    time.append(hf.get('time')[0])
    hf.close() 

    hf = h5py.File('../ori/'+dir+'.h5', 'r')
    rhouu_ori.append(np.average(np.array(hf.get('f11_Re'))))
    hf.close() 


diff_rho_uu=np.absolute(np.array(rhouu_per)[1:]-np.array(rhouu_ori)[1:])
time=time[1:]

plt.plot(time,diff_rho_uu) 

initialvalue=5
lastvalue=8

lyapunov=(1/(time[lastvalue]-time[initialvalue]))*np.log(diff_rho_uu[lastvalue]/diff_rho_uu[initialvalue])
fitdata=diff_rho_uu[initialvalue]*np.exp(lyapunov*(time[initialvalue:lastvalue+1]-time[initialvalue]))

plt.plot(time[initialvalue:lastvalue+1],fitdata,linestyle='dashed',label='$\lambda=$'+str(lyapunov))

plt.legend()
plt.yscale('log')  
plt.ylabel(r'$\Delta \rho_{\mu\mu}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/fit_difference_average_rho_uu.pdf')   
plt.clf()

#########################################################################################
# plot the difference between rho_tt original and rho_tt perturbed both average over all particles
#########################################################################################

directories = sorted(glob.glob("plt*/neutrinos"))        
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"  

time=[]
rhott_per=[]
rhott_ori=[]

for dir in directories:

    hf = h5py.File(dir+'.h5', 'r')
    rhott_per.append(np.average(np.array(hf.get('f22_Re'))))
    time.append(hf.get('time')[0])
    hf.close() 

    hf = h5py.File('../ori/'+dir+'.h5', 'r')
    rhott_ori.append(np.average(np.array(hf.get('f22_Re'))))
    hf.close() 

diff_rho_tt=np.absolute(np.array(rhott_per)[1:]-np.array(rhott_ori)[1:])
time=time[1:]

plt.plot(time,diff_rho_tt) 

initialvalue=5
lastvalue=8

lyapunov=(1/(time[lastvalue]-time[initialvalue]))*np.log(diff_rho_tt[lastvalue]/diff_rho_tt[initialvalue])
fitdata=diff_rho_tt[initialvalue]*np.exp(lyapunov*(time[initialvalue:lastvalue+1]-time[initialvalue]))

plt.plot(time[initialvalue:lastvalue+1],fitdata,linestyle='dashed',label='$\lambda=$'+str(lyapunov))

plt.legend()
plt.yscale('log')  
plt.ylabel(r'$\Delta \rho_{\tau\tau}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/fit_difference_average_rho_tt.pdf')   
plt.clf()






#########################################################################################
# plot the difference between rho_eu original and rho_eu perturbed both average over all particles
#########################################################################################

directories = sorted(glob.glob("plt*/neutrinos"))        
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"  

time=[]
rhoeu_per=[]
rhoeu_ori=[]

for dir in directories:

    hf = h5py.File(dir+'.h5', 'r')
    rhoeu_per.append(np.average(np.sqrt(np.array(hf.get('f01_Re'))**2+np.array(hf.get('f01_Im'))**2)))
    time.append(hf.get('time')[0])
    hf.close() 

    hf = h5py.File('../ori/'+dir+'.h5', 'r')
    rhoeu_ori.append(np.average(np.sqrt(np.array(hf.get('f01_Re'))**2+np.array(hf.get('f01_Im'))**2)))
    hf.close() 

diff_rho_eu=np.absolute(np.array(rhoeu_per)[1:]-np.array(rhoeu_ori)[1:])
time=time[1:]

plt.plot(time,diff_rho_eu) 

initialvalue=0
lastvalue=3

lyapunov=(1/(time[lastvalue]-time[initialvalue]))*np.log(diff_rho_eu[lastvalue]/diff_rho_eu[initialvalue])
fitdata=diff_rho_eu[initialvalue]*np.exp(lyapunov*(time[initialvalue:lastvalue+1]-time[initialvalue]))

plt.plot(time[initialvalue:lastvalue+1],fitdata,linestyle='dashed',label='$\lambda=$'+str(lyapunov))

plt.legend()
plt.yscale('log')  
plt.ylabel(r'$\Delta \rho_{e\mu}$')
plt.xlabel(r'Time (s)')
plt.savefig('plots/fit_difference_average_rho_eu.pdf')   
plt.clf()





