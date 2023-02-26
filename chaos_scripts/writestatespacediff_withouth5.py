##########################################################
#This script write all the particle information in the plt* directories into hdf5 format files
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
nproc = 1

directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"

# get NF
eds = emu.EmuDataset(directories[0])
NF = eds.get_num_flavors()
if NF==2:
    rkey, ikey = amrex.get_particle_keys()
    labels=['pos_x','pos_y','pos_z', 'time', 'x', 'y', 'z', 'pupx', 'pupy', 'pupz', 'pupt', 'N', 'L', 'f00_Re', 'f01_Re', 'f01_Im', 'f11_Re', 'Nbar', 'Lbar', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f11_Rebar']
if NF==3:
    rkey, ikey = amrex.get_3flavor_particle_keys()
    labels=['pos_x','pos_y','pos_z','time','x', 'y', 'z', 'pupx', 'pupy', 'pupz', 'pupt', 'N', 'L', 'f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im' ,'f22_Re', 'Nbar' ,'Lbar', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar' ,'f12_Imbar', 'f22_Rebar']
    
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

def writehdf5files(dire):

    eds = emu.EmuDataset(dire)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    header = amrex.AMReXParticleHeader(dire+"/neutrinos/Header")
    grid_data = GridData(ad)
    nlevels = len(header.grids)
    assert nlevels==1
    level = 0
    ngrids = len(header.grids[level])                                                                                                                        
    
    data_given=[]
    data_per=[]
    
    number_of_particles=0
    
    for gridID in range(ngrids):
        
        idata_given, rdata_given = amrex.read_particle_data('../'+dire, ptype="neutrinos", level_gridID=(level,gridID))
        idata_per, rdata_per = amrex.read_particle_data(dire, ptype="neutrinos", level_gridID=(level,gridID))
        
        number_of_particles=number_of_particles+len(rdata_given)
        
        data_given.append(rdata_given)
        data_per.append(rdata_per)
    
    data_given=np.array(data_given)
    data_per=np.array(data_per)

    number_of_particles_variables=len(data_given[0][0])
    
    data_given=np.reshape(data_given,(number_of_particles,number_of_particles_variables))
    data_per=np.reshape(data_per,(number_of_particles,number_of_particles_variables))
    
    labels_to_compare=[rkey['pos_x'],rkey['pos_y'],rkey['pos_z'],rkey['time'],rkey['pupx'],rkey['pupy'],rkey['pupz'],rkey['pupt']]

    identifier_given=[]
    identifier_per=[]

    for i in range(len(data_given)):
        
        given_string=''
        per_string=''
        
        for j in labels_to_compare:
            given_string=given_string+str(data_given[i][j])
            per_string=per_string+str(data_per[i][j])

        identifier_given.append(given_string)
        identifier_per.append(per_string)
    
    identifier_given=np.array(identifier_given)
    identifier_per=np.array(identifier_per)
    
    index_given=np.argsort(identifier_given)
    index_per=np.argsort(identifier_per)

    assert np.all(identifier_given[index_given]==identifier_per[index_per]),'\n \n ---> ending execution: particles do not match \n \n'

    data_given=data_given[index_given]
    data_per=data_per[index_per]

    return dire

# run the write hdf5 files function in parallel
if __name__ == '__main__':
    pool = Pool(nproc)
    finalresult=pool.map(writehdf5files,directories)
    for i in finalresult: print("completed ---> "+i)
