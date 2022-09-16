import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
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

##########
# INPUTS #
##########

nproc = 2

#########################
# loop over directories #
#########################

dirori = sorted(glob.glob("ori/plt*/neutrinos"))
dirori = [dirori[i].split('/n')[0] for i in range(len(dirori))] # remove "neutrinos"

dirper = sorted(glob.glob("per/plt*/neutrinos"))
dirper = [dirper[i].split('/n')[0] for i in range(len(dirper))] # remove "neutrinos"

# get NF
eds = emu.EmuDataset(dirori[0])
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


#function that given a plt file (original and perturbed for the same time), compute the diference between state spaces vectors

def compute_ssdiff(dire):
   
    print('.............................................')
    print(dire)
    print('.............................................')
    
    dirori_='ori/'+str(dire)
    dirper_='per/'+str(dire)

    ##################################
    #oridata
    ##################################

    eds = emu.EmuDataset(dirori_)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    ################
    # angular work #
    ################
    header = amrex.AMReXParticleHeader(dirori_+"/neutrinos/Header")
    grid_data = GridData(ad)
    nlevels = len(header.grids)
    assert nlevels==1
    level = 0
    ngrids = len(header.grids[level])

    # average the angular power spectrum over many cells
    # loop over all cells within each grid
    

    oridata=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    for gridID in range(ngrids):
        #print("grid",gridID+1,"/",ngrids)
        
        # read particle data on a single grid
        idata, rdata = amrex.read_particle_data(dirori_, ptype="neutrinos", level_gridID=(level,gridID))
        
        for i in rdata:

            oridata[0].append(str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[7])+str(i[8])+str(i[9]))
            oridata[1].append(i[13])
            oridata[2].append(i[14])
            oridata[3].append(i[15])
            oridata[4].append(i[16])
            oridata[5].append(i[17])
            oridata[6].append(i[18])
            oridata[7].append(i[19])
            oridata[8].append(i[20])
            oridata[9].append(i[21])
            oridata[10].append(i[24])
            oridata[11].append(i[25])
            oridata[12].append(i[26])
            oridata[13].append(i[27])
            oridata[14].append(i[28])
            oridata[15].append(i[29])
            oridata[16].append(i[30])
            oridata[17].append(i[31])
            oridata[18].append(i[32])
    
    idx=np.argsort(oridata[0])
    for m in range(0,len(oridata)):    
        oridata[m]=np.array(oridata[m])[idx]
    
    ##################################
    #perdata
    ##################################

    eds = emu.EmuDataset(dirper_)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    ################
    # angular work #
    ################
    header = amrex.AMReXParticleHeader(dirper_+"/neutrinos/Header")
    grid_data = GridData(ad)
    nlevels = len(header.grids)
    assert nlevels==1
    level = 0
    ngrids = len(header.grids[level])

    perdata=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    # average the angular power spectrum over many cells
    # loop over all cells within each grid
    for gridID in range(ngrids):
        #print("grid",gridID+1,"/",ngrids)
        
        # read particle data on a single grid
        idata, rdata = amrex.read_particle_data(dirper_, ptype="neutrinos", level_gridID=(level,gridID))
        
        for i in rdata:

            perdata[0].append(str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[7])+str(i[8])+str(i[9]))
            perdata[1].append(i[13])
            perdata[2].append(i[14])
            perdata[3].append(i[15])
            perdata[4].append(i[16])
            perdata[5].append(i[17])
            perdata[6].append(i[18])
            perdata[7].append(i[19])
            perdata[8].append(i[20])
            perdata[9].append(i[21])
            perdata[10].append(i[24])
            perdata[11].append(i[25])
            perdata[12].append(i[26])
            perdata[13].append(i[27])
            perdata[14].append(i[28])
            perdata[15].append(i[29])
            perdata[16].append(i[30])
            perdata[17].append(i[31])
            perdata[18].append(i[32])

    idx=np.argsort(perdata[0])
    for m in range(0,len(perdata)):
        perdata[m]=np.array(perdata[m])[idx]
    
    tiempo=rdata[0][3]

    # test leave the code in case the perturbed original particles positions and momentum do not match

    for n in range(0,len(oridata[0])): 
        assert oridata[0][n]==perdata[0][n],'the particles do not match'

    sum_state_space=0.0
    
    numbers=[1,2,2,3,3,4,4,5,5,6,7,7,8,8,9]

    for j in range(0,len(oridata[0])):
        for k in numbers:
            sum_state_space=(oridata[k][j]-perdata[k][j])**2+sum_state_space

    sum_state_space=np.sqrt(sum_state_space)
    
    return tiempo,sum_state_space

dir_ = [dirori[i].split('/')[1] for i in range(len(dirori))]

if __name__ == '__main__':
    pool = Pool(nproc)
    #finalresult=pool.map(compute_ssdiff,dir_)  
    finalresult=pool.map(compute_ssdiff,['plt00000','plt00025']) 
    f = open("ssreduction.txt", "w")
    f.write("time ssdiff")
    for e in finalresult:
        f.write(str(e[0])+" "+str(e[1])+" \n")
    f.close()
        
