import matplotlib.pyplot as plt                                                                                                                  
import numpy as np                                                                                                                               
import h5py             
import glob

# read the names of perturbed data files

directories = sorted(glob.glob("plt*/neutrinos"))  
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"                                                                                                              

# saving the keys for read and analyze the data

keys=['f00_Re','f01_Re','f01_Im','f02_Re','f02_Im','f11_Re','f12_Re','f12_Im','f22_Re','f00_Rebar','f01_Rebar','f01_Imbar','f02_Rebar','f02_Imbar','f11_Rebar','f12_Rebar','f12_Imbar','f22_Rebar']
keys_test=['pos_x','pos_y','pos_z','pupx','pupy','pupz','time']

# creating the variables for save the time and the magnitud of the difference of the state space vector

time=[]
ssvecdiff=[]

# loopping over the perturbed data

for directory in directories[1:]:
    
    # opening the hdf5 files

    hper = h5py.File(directory+'.h5', 'r') 
    hori = h5py.File('../ori/'+directory+'.h5', 'r')

    # generate a unique identifier for each particle for the original and perturbed data
    
    identifier_ori=np.asarray(hori.get('pos_x'), dtype=str)
    identifier_per=np.asarray(hper.get('pos_x'), dtype=str)
    
    for keyt in keys_test[1:]:
        identifier_ori=np.char.add(identifier_ori,np.asarray(hori.get(keyt), dtype=str))
        identifier_per=np.char.add(identifier_per,np.asarray(hper.get(keyt), dtype=str))
    
    # test if the particles in the original and perturbed file match time, position and momentum

    if not(np.all(identifier_ori==identifier_per)):

        print(directory+' - error: particles do not match')
        
        index_ori=np.argsort(identifier_ori)
        index_per=np.argsort(identifier_per)
        
        # test again if the particles match, if not end the execution

        assert np.all(np.array(hori.get('pos_x'))[index_ori]==np.array(hper.get('pos_x'))[index_per]), '\n \n '+directory+' ---> ending execution: x position do not match \n \n' 
        assert np.all(np.array(hori.get('pos_y'))[index_ori]==np.array(hper.get('pos_y'))[index_per]), '\n \n '+directory+' ---> ending execution: y position do not match \n \n'   
        assert np.all(np.array(hori.get('pos_z'))[index_ori]==np.array(hper.get('pos_z'))[index_per]), '\n \n '+directory+' ---> ending execution: z position do not match \n \n'
        assert np.all(np.array(hori.get('pupx'))[index_ori]==np.array(hper.get('pupx'))[index_per]), '\n \n '+directory+' ---> ending execution: x momentum do not match \n \n'  
        assert np.all(np.array(hori.get('pupy'))[index_ori]==np.array(hper.get('pupy'))[index_per]), '\n \n '+directory+' ---> ending execution: y momentum do not match \n \n' 
        assert np.all(np.array(hori.get('pupz'))[index_ori]==np.array(hper.get('pupz'))[index_per]), '\n \n '+directory+' ---> ending execution: z momentum do not match \n \n' 
        assert np.all(np.array(hori.get('time'))[index_ori]==np.array(hper.get('time'))[index_per]), '\n \n '+directory+' ---> ending execution: t momentum do not match \n \n'   

        # getting the time value

        time.append(hori.get('time')[0])

        # computing the magnitud of the difference of the state space vector

        ssdiff=0

        for key in keys:
            ssdiff=ssdiff+np.sum(np.square(np.array(hper.get(key))[index_per]-np.array(hori.get(key))[index_ori]))

        ssvecdiff.append(np.sqrt(ssdiff))

    else:

        print(directory)

        # getting the time value

        time.append(hori.get('time')[0])

        # computing the magnitud of the difference of the state space vector

        ssdiff=0

        for key in keys:
            ssdiff=ssdiff+np.sum(np.square(np.array(hper.get(key))-np.array(hori.get(key))))

        ssvecdiff.append(np.sqrt(ssdiff))

    # closing the hdf5 files

    hper.close()
    hori.close()


# plotting magnitud of the state space vectors differences

plt.plot(time,ssvecdiff) 
plt.yscale('log')  
plt.ylabel(r'$\left | \Delta \vec{r}_{ss} \right | $')
plt.xlabel(r'Time (s)')
plt.savefig('plots/difference_state_space_vector.pdf')   
plt.clf() 




