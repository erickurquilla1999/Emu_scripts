import matplotlib.pyplot as plt                                                                                                                  
import numpy as np                                                                                                                               
import h5py             
import glob
from scipy.signal import find_peaks
from multiprocessing import Pool

nproc = 16

# loopping over the perturbed data

def cssdiff(directory):

    # saving the keys for read and analyze the data

    keys=['f00_Re','f01_Re','f01_Im','f02_Re','f02_Im','f11_Re','f12_Re','f12_Im','f22_Re','f00_Rebar','f01_Rebar','f01_Imbar','f02_Rebar','f02_Imbar','f11_Rebar','f12_Rebar','f12_Imbar','f22_Rebar']
    keys_test=['pos_x','pos_y','pos_z','pupx','pupy','pupz','time']

    # opening the hdf5 files

    hper = h5py.File(directory+'.h5', 'r') 
    hori = h5py.File('../'+directory+'.h5', 'r')

    # generate a unique identifier for each particle for the original and perturbed data
    
    identifier_ori=np.asarray(hori.get('pos_x'), dtype=str)
    identifier_per=np.asarray(hper.get('pos_x'), dtype=str)
    
    for keyt in keys_test[1:]:
        identifier_ori=np.char.add(identifier_ori,np.asarray(hori.get(keyt), dtype=str))
        identifier_per=np.char.add(identifier_per,np.asarray(hper.get(keyt), dtype=str))
    
    # test if the particles in the original and perturbed file match time, position and momentum

    if not(np.all(identifier_ori==identifier_per)):
        
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

        time=hori.get('time')[0]

        # computing the magnitud of the difference of the state space vector

        ssmag_per=0
        ssmag_ori=0
        ssdiff=0

        for key in keys:
            ssmag_per=ssmag_per+np.sum(np.square(np.array(hper.get(key))[index_per]))
            ssmag_ori=ssmag_ori+np.sum(np.square(np.array(hori.get(key))[index_ori]))            
            ssdiff=ssdiff+np.sum(np.square(np.array(hper.get(key))[index_per]-np.array(hori.get(key))[index_ori]))

        ssvecdiff=np.sqrt(ssdiff)
        ssmagper=np.sqrt(ssmag_per)
        ssmagori=np.sqrt(ssmag_ori)
        
        print()
        print(directory+' - error: particles do not match')
        print(ssdiff)
        print(str(np.sqrt(ssmag_per))+' ---> '+str(100*np.sqrt(ssdiff)/np.sqrt(ssmag_per))+' % ')
        print(str(np.sqrt(ssmag_ori))+' ---> '+str(100*np.sqrt(ssdiff)/np.sqrt(ssmag_ori))+' % ')
        print()

        # closing the hdf5 files

        hper.close()
        hori.close()

        return [time,ssvecdiff,ssmagper,ssmagori]
        
    else:

        # getting the time value

        time=hori.get('time')[0]

        # computing the magnitud of the difference of the state space vector

        ssmag_per=0
        ssmag_ori=0
        ssdiff=0
        
        for key in keys:
            ssmag_per=ssmag_per+np.sum(np.square(np.array(hper.get(key))))
            ssmag_ori=ssmag_ori+np.sum(np.square(np.array(hori.get(key))))            
            ssdiff=ssdiff+np.sum(np.square(np.array(hper.get(key))-np.array(hori.get(key))))

        ssvecdiff=np.sqrt(ssdiff)
        ssmagper=np.sqrt(ssmag_per)
        ssmagori=np.sqrt(ssmag_ori)
        
        print()
        print(directory)
        print(ssdiff)
        print(str(np.sqrt(ssmag_per))+' ---> '+str(100*np.sqrt(ssdiff)/np.sqrt(ssmag_per))+' % ')
        print(str(np.sqrt(ssmag_ori))+' ---> '+str(100*np.sqrt(ssdiff)/np.sqrt(ssmag_ori))+' % ')
        print()

        # closing the hdf5 files

        hper.close()
        hori.close()

        return [time,ssvecdiff,ssmagper,ssmagori]


# read the names of perturbed data files

directories = sorted(glob.glob("plt*/neutrinos"))  
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"                                                                                                              

# run the cssdiff files function in parallel
if __name__ == '__main__':
    pool = Pool(nproc)
    finalresult=pool.map(cssdiff,directories)

t=[]
ssv=[]
ssmp=[]
ssmo=[]

for f in finalresult:
    t.append(f[0])
    ssv.append(f[1])
    ssmp.append(f[2])
    ssmo.append(f[3])

ind=np.argsort(np.array(t))

# writting a hdf5 file for save the generated data
hf = h5py.File("state_space_difference.h5", 'w')
hf.create_dataset("time",data=np.array(t)[ind])
hf.create_dataset("statespacediff",data=np.array(ssv)[ind])
hf.create_dataset("ssmagper",data=np.array(ssmp)[ind])
hf.create_dataset("ssmagori",data=np.array(ssmo)[ind])
hf.close()