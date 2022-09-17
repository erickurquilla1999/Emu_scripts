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

#start counting time
st = time.time()

# get the plt directories of the original simulation
directories = sorted(glob.glob("ori/plt*/neutrinos"))
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
        

    # particle cell id ON THE CURRENT GRID
    # the x, y, and z values are assumed to be relative to the
    # lower boundary of the grid
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
    dirori_='ori/'+str(dire)
    dirper_='per/'+str(dire)

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
        
        # test if the particles in the original and perturbed file match position and momentum
        assert np.all(rdataori[:,rkey["pos_x"]] == rdataper[:,rkey["pos_x"]]), 'x position do not match'
        assert np.all(rdataori[:,rkey["pos_y"]] == rdataper[:,rkey["pos_y"]]), 'y position do not match'
        assert np.all(rdataori[:,rkey["pos_z"]] == rdataper[:,rkey["pos_z"]]), 'z position do not match'
        assert np.all(rdataori[:,rkey["pupx"]] == rdataper[:,rkey["pupx"]]), 'x momentum do not match'
        assert np.all(rdataori[:,rkey["pupy"]] == rdataper[:,rkey["pupy"]]), 'y momentum do not match'
        assert np.all(rdataori[:,rkey["pupz"]] == rdataper[:,rkey["pupz"]]), 'z momentum do not match'

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

    #print the total execution time of this function
    print(str(dire)+" time: "+str(time.time() - cpt))
    
    return rdataori[:,rkey["time"]][1],sumss**0.5

directories=[directories[i].split('/')[1] for i in range(len(directories))]
directories=directories[1:-1]

nproc = 2

if __name__ == '__main__':

    pool = Pool(nproc)
    
    finalresult=pool.map(compute_ssdiff,directories)  
    
    # writting a text file for save the generated data
    f = open("ssreduction.txt", "w")
    f.write("time ssdiff \n")
    for e in finalresult:
        f.write(str(e[0])+" "+str(e[1])+" \n")
    f.close()

# print total execution time of the program
print("\n \n Total time: "+str(time.time() - st)+" \n \n")
