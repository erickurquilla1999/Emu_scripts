##########################################################
# this script save a .h5 file with the magnitud of the 
# difference of state vectors for a given time for two 
# differente simulations
##########################################################

import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/data_reduction')
import numpy as np
import matplotlib.pyplot as plt
import yt
import glob
import multiprocessing as mp
import h5py
import amrex_plot_tools as amrex
import emu_yt_module as emu
from multiprocessing import Pool
import scipy.special
import time

##########
# INPUTS #
##########
nproc = 10         

#start counting time
st = time.time()

# get the plt directories of the original simulation
directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/n')[0] for i in range(len(directories))] # remove "neutrinos"

# get number of flavors ---> NF
eds = emu.EmuDataset(directories[0])
NF = eds.get_num_flavors()
if NF==2:
    rkey, ikey = amrex.get_particle_keys()
if NF==3:
    rkey, ikey = amrex.get_3flavor_particle_keys()

class GridData(object):
    def __init__(self, ad):
        x = ad['index','x'].d
        y = ad['index','y'].d
        z = ad['index','z'].d
        dx = ad['index','dx'].d
        dy = ad['index','dy'].d
        dz = ad['index','dz'].d
        self.ad = ad
        self.dx = dx[0]
        self.dy = dy[0]
        self.dz = dz[0]
        self.xmin = np.min(x-dx/2.)
        self.ymin = np.min(y-dy/2.)
        self.zmin = np.min(z-dz/2.)
        self.xmax = np.max(x+dx/2.)
        self.ymax = np.max(y+dy/2.)
        self.zmax = np.max(z+dz/2.)
        self.nx = int((self.xmax - self.xmin) / self.dx + 0.5)
        self.ny = int((self.ymax - self.ymin) / self.dy + 0.5)
        self.nz = int((self.zmax - self.zmin) / self.dz + 0.5)
        #print(self.nx, self.ny, self.nz)
       
    def get_particle_cell_ids(self,rdata):
        # get coordinates
        x = rdata[:,rkey["x"]]
        y = rdata[:,rkey["y"]]
        z = rdata[:,rkey["z"]]
        ix = (x/self.dx).astype(int)
        iy = (y/self.dy).astype(int)
        iz = (z/self.dz).astype(int)

        # HACK - get this grid's bounds using particle locations
        ix -= np.min(ix)
        iy -= np.min(iy)
        iz -= np.min(iz)
        nx = np.max(ix)+1
        ny = np.max(iy)+1
        nz = np.max(iz)+1
        idlist = (iz + nz*iy + nz*ny*ix).astype(int)

        return idlist

#function that given a plt 2 files (original and perturbed) for the same time, compute the diference between state spaces vectors
def compute_ssdiff(dire):
    
    #start counting time
    cpt = time.time()
   
    #directories for original and perturbed data
    dirori_='../ori/'+str(dire)
    dirper_=dire

    eds = emu.EmuDataset(dirori_)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    # get the number of grid cells 
    header = amrex.AMReXParticleHeader(dirori_+"/neutrinos/Header")
    grid_data = GridData(ad)
    nlevels = len(header.grids)
    assert nlevels==1
    level = 0
    ngrids = len(header.grids[level])
    
    # sumss will save the sum the squeres of the difference between each state space vector components 
    sumss=0

    for gridID in range(ngrids):
        
        # read particle data on a single grid
        idataori, rdataori = amrex.read_particle_data(dirori_, ptype="neutrinos", level_gridID=(level,gridID))
        idataper, rdataper = amrex.read_particle_data(dirper_, ptype="neutrinos", level_gridID=(level,gridID))
        
        # sort the rdataori 
        idori=[]
        for p in rdataori:
            iden=''
            for pdata in p:
                iden=iden+str(pdata)
            idori.append(iden)

        idx_ori=np.argsort(idori) 
        rdataori=np.array(rdataori)[idx_ori]
        
        # sort the rdataper
        idper=[]
        for p in rdataper:
            iden=''
            for pdata in p:
                iden=iden+str(pdata)
            idper.append(iden)

        idx_per=np.argsort(idper) 
        rdataper=np.array(rdataper)[idx_per]

        # test if the particles in the original and perturbed file match position and momentum
        assert np.all(rdataori[:,rkey["pos_x"]] == rdataper[:,rkey["pos_x"]]), '\n \n '+dire+' ---> x position do not match \n \n'
        assert np.all(rdataori[:,rkey["pos_y"]] == rdataper[:,rkey["pos_y"]]), '\n \n '+dire+' ---> y position do not match \n \n'
        assert np.all(rdataori[:,rkey["pos_z"]] == rdataper[:,rkey["pos_z"]]), '\n \n '+dire+' ---> z position do not match \n \n'
        assert np.all(rdataori[:,rkey["pupx"]] == rdataper[:,rkey["pupx"]]), '\n \n '+dire+' ---> x momentum do not match \n \n'
        assert np.all(rdataori[:,rkey["pupy"]] == rdataper[:,rkey["pupy"]]), '\n \n '+dire+' ---> y momentum do not match \n \n'
        assert np.all(rdataori[:,rkey["pupz"]] == rdataper[:,rkey["pupz"]]), '\n \n '+dire+' ---> z momentum do not match \n \n'
        assert np.all(rdataori[:,rkey["time"]] == rdataper[:,rkey["time"]]), '\n \n '+dire+' ---> t momentum do not match \n \n' 
    
        keys=['f00_Re','f01_Re','f01_Im','f02_Re','f02_Im','f11_Re','f12_Re','f12_Im','f22_Re','f00_Rebar','f01_Rebar','f01_Imbar','f02_Rebar','f02_Imbar','f11_Rebar','f12_Rebar','f12_Imbar','f22_Rebar']

        for key in keys:
            sumss=sumss+np.sum(np.square(rdataori[:,rkey[key]]-rdataper[:,rkey[key]]))

        '''
        #loop over each components of the density matrix and sum the squere of the difference
        #of each componentes for the original and perturbed simulation
        for a in range(NF):
            for b in range(a,NF):
                nodiagonal=1
                im="_Im"
                if(a==b):
                    nodiagonal=0
                    im="_Re"
                if(b>a):coefficient = 2
                else:coefficient = 1
                label = "f"+str(a)+str(b)+"_Re"
                sumss+=coefficient*np.sum((rdataori[:,rkey[label]]-rdataper[:,rkey[label]])**2)
                label = "f"+str(a)+str(b)+im
                sumss+=nodiagonal*coefficient*np.sum((rdataori[:,rkey[label]]-rdataper[:,rkey[label]])**2)
        '''

    #print the total execution time of this function
    print(str(dire)+" time: "+str(time.time() - cpt))

    return rdataori[:,rkey["time"]][1],sumss**0.5

if __name__ == '__main__':

    # run compute_ssdiff in  paralell
    pool = Pool(nproc)   
    finalresult=pool.map(compute_ssdiff,directories)  
    
    t=[]
    sumss=[]
    for results in finalresult:
        t.append(results[0])
        sumss.append(results[1])

    # writting a hdf5 file for save the generated data
    hf = h5py.File("statespacediff.h5", 'w')
    hf.create_dataset("time",data=t)
    hf.create_dataset("statespacediff",data=sumss)
    hf.close()

    print('\n \n')
    for i in directories:
        print('Done '+i)

# print total execution time of the program
print("\n \n Total time: "+str(time.time() - st)+" \n \n")
