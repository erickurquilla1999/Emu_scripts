##########################################################
#This script write all the particle information in the plt* directories into .txt files
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

##########
# INPUTS #
##########
nproc = 2

directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"

# get NF
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
        print(self.nx, self.ny, self.nz)
        
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

def writetxtfiles(dire):

    eds = emu.EmuDataset(dire)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    header = amrex.AMReXParticleHeader(dire+"/neutrinos/Header")
    grid_data = GridData(ad)
    nlevels = len(header.grids)
    assert nlevels==1
    level = 0
    ngrids = len(header.grids[level])

    #creating the file to save the particle data
    file1 = open(str(dire)+".txt","w")
    file1.write('pos_x pos_y pos_z time x y z pupx pupy pupz pupt N L f00_Re f01_Re f01_Im f02_Re f02_Im f11_Re f12_Re f12_Im f22_Re Nbar Lbar f00_Rebar f01_Rebar f01_Imbar f02_Rebar f02_Imbar f11_Rebar f12_Rebar f12_Imbar f22_Rebar \n') 

    # loop over all cells within each grid
    for gridID in range(ngrids):
        
        # read particle data on a single grid
        idata, rdata = amrex.read_particle_data(dire, ptype="neutrinos", level_gridID=(level,gridID))
        # writing the particle data 
        for i in rdata:
            file1.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+" "+str(i[3])+" "+str(i[4])+" "+str(i[5])+" "+str(i[6])+" "+str(i[7])+" "+str(i[8])+" "+str(i[9])+" "+str(i[10])+" "+str(i[11])+" "+str(i[12])+" "+str(i[13])+" "+str(i[14])+" "+str(i[15])+" "+str(i[16])+" "+str(i[17])+" "+str(i[18])+" "+str(i[19])+" "+str(i[20])+" "+str(i[21])+" "+str(i[22])+" "+str(i[23])+" "+str(i[24])+" "+str(i[25])+" "+str(i[26])+" "+str(i[27])+" "+str(i[28])+" "+str(i[29])+" "+str(i[30])+" "+str(i[31])+" "+str(i[31])+" \n")

    return dire

# run the writetxtfiles function in parallel
if __name__ == '__main__':
    pool = Pool(nproc)
    finalresult=pool.map(writetxtfiles,directories)
    for i in finalresult: print("completed ---> "+i)
    