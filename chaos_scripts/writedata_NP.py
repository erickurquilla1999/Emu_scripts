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
do_average=1
do_ssdiff=1
particle_index=1

##########
# creating a list of directories of the data to be read #
##########

directories_h5 = sorted(glob.glob("*.h5"))
directories_h5 = [directories_h5[i].split('.h5')[0] for i in range(len(directories_h5))]
directories_all = sorted(glob.glob("plt*/neutrinos"))
directories_all = [directories_all[i].split('/')[0] for i in range(len(directories_all))]

directories=[]

for dir1 in directories_all:
    thereis=0
    for dir2 in directories_h5:
        if dir1==dir2:
            thereis=1
            break
    if thereis==0:
        directories.append(dir1)

assert len(directories)!=0,'\n \n ---> ending execution: plt*.h5 already computed \n \n'

##########
# for read emu data stuff #
##########

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

##########
# function that do the work #
##########

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

    #########################################################################
    # do average and ss diff
    #########################################################################

    if do_average==1 and do_ssdiff==1:

        #these list will save the particles data
        data_given=[]
        data_per=[]

        #variable to count the number of particles
        number_of_particles=0

        #loop over the grid cells and save al the particles data 
        for gridID in range(ngrids):
            
            idata_given, rdata_given = amrex.read_particle_data('../'+dire.split('.')[0], ptype="neutrinos", level_gridID=(level,gridID))
            idata_per, rdata_per = amrex.read_particle_data(dire, ptype="neutrinos", level_gridID=(level,gridID))
            
            number_of_particles=number_of_particles+len(rdata_given)
            
            data_given.append(rdata_given)
            data_per.append(rdata_per)

            #delete the data to save memory 
            del idata_given
            del rdata_given
            del idata_per
            del rdata_per

        data_given=np.array(data_given)
        data_per=np.array(data_per)

        #save the number of variables that were stores for each particle
        number_of_particles_variables=len(data_given[0][0])

        #reshape the data array in order to be a single array (not arrays of arrays)
        data_given=np.reshape(data_given,(number_of_particles,number_of_particles_variables))
        data_per=np.reshape(data_per,(number_of_particles,number_of_particles_variables))

        #delete some varibles to save memory 
        del number_of_particles
        del number_of_particles_variables

        #these labels will be used for create a unique identifies for each particle in the simulation
        labels_to_compare=[rkey['N'],rkey['Nbar'],rkey['time'],rkey['pupx'],rkey['pupy'],rkey['pupz'],rkey['pupt'],rkey['pos_x'],rkey['pos_y'],rkey['pos_z']]

        #list that will store each particle identifier
        identifier_given=[]
        identifier_per=[]

        #creating the indentifiers
        for i in range(len(data_given)):
            
            given_string=''
            per_string=''
            
            for j in labels_to_compare:
                given_string=given_string+str(data_given[i][j])+' '
                per_string=per_string+str(data_per[i][j])+' '

            identifier_given.append(given_string)
            identifier_per.append(per_string)

        identifier_given=np.array(identifier_given)
        identifier_per=np.array(identifier_per)

        #creating an index to sorting the data according to the identifiers
        index_given=np.argsort(identifier_given)
        index_per=np.argsort(identifier_per)
        
        #sorting the identifiers
        identifier_given=identifier_given[index_given]
        identifier_per=identifier_per[index_per]
        
        #sorting the data
        data_given=data_given[index_given]
        data_per=data_per[index_per]

        #stop the code is the sort data does match each other (given data and perturbed data)
        if not np.all(identifier_given[index_given]==identifier_per[index_per]):
            cond=[]
            for ind in labels_to_compare:
                cond.append(np.allclose(data_per[:,ind],data_given[:,ind], rtol=1e-8, atol=1e-8, equal_nan=False))
            if not np.all(cond):
                assert True,'\n \n ---> ending execution: particles do not match in file '+dire+' \n \n'

        #deleting some varibles to save memory
        del identifier_given
        del identifier_per
        del index_given
        del index_per
        del labels_to_compare

        #creating a hdf5 diles to save the particles data
        hf = h5py.File(str(dire)+".h5", 'w')

        #saving time
        hf.create_dataset('time', data=t)

        #saving all of the information of a single particle
        hf.create_dataset('single_particle_given', data=data_given[particle_index])
        hf.create_dataset('single_particle_per', data=data_per[particle_index])

        #this keys will be used to compute the average of all the components of the density matrices
        keys_for_average=[['f00_Re'], ['f01_Re', 'f01_Im'], ['f02_Re', 'f02_Im'], ['f11_Re'], ['f12_Re', 'f12_Im'] ,['f22_Re'], ['f00_Rebar'], ['f01_Rebar', 'f01_Imbar'], ['f02_Rebar', 'f02_Imbar'],['f11_Rebar'],['f12_Rebar' ,'f12_Imbar'],['f22_Rebar']]
        #keys_for_average=[['f00_Re'], ['f01_Re', 'f01_Im'], ['f11_Re'], ['f00_Rebar'], ['f01_Rebar', 'f01_Imbar'], ['f11_Rebar']]

        #count the index of the array keys_for_average in the next loop
        countforNorNbar=0

        #loogping over all the keys
        for key in keys_for_average:

            # determine what key use if N for number of neutrino or Nbar for number of antineutrinos
            if countforNorNbar<=5:
                NorNbarKey='N'
            if countforNorNbar>5:
                NorNbarKey='Nbar'

            #computed the average value over the entire domaim of the diagonal components of the density matrix
            if len(key)==1: 
                hf.create_dataset(key[0]+'_given', data=np.average(data_given[:,rkey[NorNbarKey]]*data_given[:,rkey[key[0]]]))
                hf.create_dataset(key[0]+'_per', data=np.average(data_per[:,rkey[NorNbarKey]]*data_per[:,rkey[key[0]]]))

            #computed the average value over the entire domaim of the non-diagonal components of the density matrix
            if len(key)==2: 
                hf.create_dataset(key[0]+'_given', data=np.average(data_given[:,rkey[NorNbarKey]]*np.sqrt(np.square(data_given[:,rkey[key[0]]])+np.square(data_given[:,rkey[key[1]]]))))
                hf.create_dataset(key[0]+'_per', data=np.average(data_per[:,rkey[NorNbarKey]]*np.sqrt(np.square(data_per[:,rkey[key[0]]])+np.square(data_per[:,rkey[key[1]]]))))

            countforNorNbar=countforNorNbar+1

        # defining gell-man matrices
        gm1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        gm2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])
        gm3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]])
        gm4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
        gm5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]])
        gm6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
        gm7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]])
        gm8 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / np.sqrt(3)

        GM=np.array([gm1,gm2,gm3,gm4,gm5,gm6,gm7,gm8])

        rho_given=np.array([[data_given[:,rkey['f00_Re']],data_given[:,rkey['f01_Re']]+np.array(+1j)*data_given[:,rkey['f01_Im']],data_given[:,rkey['f02_Re']]+np.array(+1j)*data_given[:,rkey['f02_Im']]],[data_given[:,rkey['f01_Re']]+np.array(-1j)*data_given[:,rkey['f01_Im']],data_given[:,rkey['f11_Re']],data_given[:,rkey['f12_Re']]+np.array(+1j)*data_given[:,rkey['f12_Im']]],[data_given[:,rkey['f02_Re']]+np.array(-1j)*data_given[:,rkey['f02_Im']],data_given[:,rkey['f12_Re']]+np.array(-1j)*data_given[:,rkey['f12_Im']],data_given[:,rkey['f22_Re']]]])
        rhobar_given=np.array([[data_given[:,rkey['f00_Rebar']],data_given[:,rkey['f01_Rebar']]+np.array(+1j)*data_given[:,rkey['f01_Imbar']],data_given[:,rkey['f02_Rebar']]+np.array(+1j)*data_given[:,rkey['f02_Imbar']]],[data_given[:,rkey['f01_Rebar']]+np.array(-1j)*data_given[:,rkey['f01_Imbar']],data_given[:,rkey['f11_Rebar']],data_given[:,rkey['f12_Rebar']]+np.array(+1j)*data_given[:,rkey['f12_Imbar']]],[data_given[:,rkey['f02_Rebar']]+np.array(-1j)*data_given[:,rkey['f02_Imbar']],data_given[:,rkey['f12_Rebar']]+np.array(-1j)*data_given[:,rkey['f12_Imbar']],data_given[:,rkey['f22_Rebar']]]])
        rho_per=np.array([[data_per[:,rkey['f00_Re']],data_per[:,rkey['f01_Re']]+np.array(+1j)*data_per[:,rkey['f01_Im']],data_per[:,rkey['f02_Re']]+np.array(+1j)*data_per[:,rkey['f02_Im']]],[data_per[:,rkey['f01_Re']]+np.array(-1j)*data_per[:,rkey['f01_Im']],data_per[:,rkey['f11_Re']],data_per[:,rkey['f12_Re']]+np.array(+1j)*data_per[:,rkey['f12_Im']]],[data_per[:,rkey['f02_Re']]+np.array(-1j)*data_per[:,rkey['f02_Im']],data_per[:,rkey['f12_Re']]+np.array(-1j)*data_per[:,rkey['f12_Im']],data_per[:,rkey['f22_Re']]]])
        rhobar_per=np.array([[data_per[:,rkey['f00_Rebar']],data_per[:,rkey['f01_Rebar']]+np.array(+1j)*data_per[:,rkey['f01_Imbar']],data_per[:,rkey['f02_Rebar']]+np.array(+1j)*data_per[:,rkey['f02_Imbar']]],[data_per[:,rkey['f01_Rebar']]+np.array(-1j)*data_per[:,rkey['f01_Imbar']],data_per[:,rkey['f11_Rebar']],data_per[:,rkey['f12_Rebar']]+np.array(+1j)*data_per[:,rkey['f12_Imbar']]],[data_per[:,rkey['f02_Rebar']]+np.array(-1j)*data_per[:,rkey['f02_Imbar']],data_per[:,rkey['f12_Rebar']]+np.array(-1j)*data_per[:,rkey['f12_Imbar']],data_per[:,rkey['f22_Rebar']]]])

        # computing polarization vectors
        P_given=[]
        Pbar_given=[]
        P_per=[]
        Pbar_per=[]

        for gm in GM:

            Pi_given=0
            Pbari_given=0
            Pi_per=0
            Pbari_per=0            

            for m in [0,1,2]:
                for n in [0,1,2]:
                    Pi_given=Pi_given+rho_given[m][n]*gm[m][n]
                    Pbari_given=Pbari_given+rhobar_given[m][n]*gm[m][n]
                    Pi_per=Pi_per+rho_per[m][n]*gm[m][n]
                    Pbari_per=Pbari_per+rhobar_per[m][n]*gm[m][n]

                P_given.append(np.real(Pi_given))
                Pbar_given.append(np.real(Pbari_given))
                P_per.append(np.real(Pi_per))
                Pbar_per.append(np.real(Pbari_per))

        del rho_given
        del rhobar_given
        del rho_per
        del rhobar_per

        #this keys will be used to compute the state space diference vector magnitud
        # keys=['f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im' ,'f22_Re', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar' ,'f12_Imbar', 'f22_Rebar']
        #keys=['f00_Re', 'f01_Re', 'f01_Im', 'f11_Re', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f11_Rebar']
        ssmag_per=0
        ssmag_ori=0
        ssdiff=0

        #count the index of the array keys in the next loop
        # countforNorNbar=0

        #computing the state space diference vector magnitud
        for h in range(8):
            
            ssmag_ori=ssmag_ori+np.sum(np.square(data_given[:,rkey['N']]*P_given[h]))+np.sum(np.square(data_given[:,rkey['Nbar']]*Pbar_given[h]))
            ssmag_per=ssmag_per+np.sum(np.square(data_per[:,rkey['N']]*P_per[h]))+np.sum(np.square(data_per[:,rkey['Nbar']]*Pbar_per[h]))
            ssdiff=ssdiff+np.sum(np.square(data_per[:,rkey['N']]*P_per[h]-data_given[:,rkey['N']]*P_given[h]))+np.sum(np.square(data_per[:,rkey['Nbar']]*Pbar_per[h]-data_given[:,rkey['Nbar']]*Pbar_given[h]))

        hf.create_dataset('difference_state_space_vector_magnitud', data=np.sqrt(ssdiff))
        hf.create_dataset('state_space_vector_magnitud_per', data=np.sqrt(ssmag_per))
        hf.create_dataset('state_space_vector_magnitud_given', data=np.sqrt(ssmag_ori))

        hf.close()
        
        # deleting the data to save memory
        del data_per
        del data_given
        del ssmag_per
        del ssmag_ori
        del ssdiff
        
        del P_given
        del Pbar_given
        del P_per
        del Pbar_per

    #########################################################################
    # do average
    #########################################################################

    if do_average==1 and do_ssdiff==0:

        #these list will save the particles data
        data_given=[]

        #variable to count the number of particles
        number_of_particles=0

        #loop over the grid cells and save al the particles data 
        for gridID in range(ngrids):
            
            idata_given, rdata_given = amrex.read_particle_data(dire, ptype="neutrinos", level_gridID=(level,gridID))
            
            number_of_particles=number_of_particles+len(rdata_given)
            
            data_given.append(rdata_given)

            #delete the data to save memory 
            del idata_given
            del rdata_given

        data_given=np.array(data_given)

        #save the number of variables that were stores for each particle
        number_of_particles_variables=len(data_given[0][0])

        #reshape the data array in order to be a single array (not arrays of arrays)
        data_given=np.reshape(data_given,(number_of_particles,number_of_particles_variables))

        #delete some varibles to save memory 
        del number_of_particles
        del number_of_particles_variables

        #creating a hdf5 diles to save the particles data
        hf = h5py.File(str(dire)+".h5", 'w')

        #saving time
        hf.create_dataset('time', data=t)

        #saving all of the information of a single particle
        hf.create_dataset('single_particle_given', data=data_given[particle_index])

        #this keys will be used to compute the average of all the components of the density matrices
        keys_for_average=[['f00_Re'], ['f01_Re', 'f01_Im'], ['f02_Re', 'f02_Im'], ['f11_Re'], ['f12_Re', 'f12_Im'] ,['f22_Re'], ['f00_Rebar'], ['f01_Rebar', 'f01_Imbar'], ['f02_Rebar', 'f02_Imbar'],['f11_Rebar'],['f12_Rebar' ,'f12_Imbar'],['f22_Rebar']]

        #loogping over all the keys
        for key in keys_for_average:
            
            #computed the average value over the entire domaim of the diagonal components of the density matrix
            if len(key)==1: 
                hf.create_dataset(key[0]+'_given', data=np.average(data_given[:,rkey[key[0]]]))

            #computed the average value over the entire domaim of the non-diagonal components of the density matrix
            if len(key)==2: 
                hf.create_dataset(key[0]+'_given', data=np.average(np.sqrt(np.square(data_given[:,rkey[key[0]]])+np.square(data_given[:,rkey[key[1]]]))))

        hf.close()
        
        # deleting the data to save memory
        del data_given

    # #########################################################################
    # # do ss diff
    # #########################################################################

    # if do_average==0 and do_ssdiff==1:

    #     data_given=[]
    #     data_per=[]

    #     number_of_particles=0

    #     for gridID in range(ngrids):
            
    #         idata_given, rdata_given = amrex.read_particle_data('../'+dire, ptype="neutrinos", level_gridID=(level,gridID))
    #         idata_per, rdata_per = amrex.read_particle_data(dire, ptype="neutrinos", level_gridID=(level,gridID))
            
    #         number_of_particles=number_of_particles+len(rdata_given)
            
    #         data_given.append(rdata_given)
    #         data_per.append(rdata_per)

    #         del idata_given
    #         del rdata_given
    #         del idata_per
    #         del rdata_per

    #     data_given=np.array(data_given)
    #     data_per=np.array(data_per)

    #     number_of_particles_variables=len(data_given[0][0])

    #     data_given=np.reshape(data_given,(number_of_particles,number_of_particles_variables))
    #     data_per=np.reshape(data_per,(number_of_particles,number_of_particles_variables))

    #     del number_of_particles
    #     del number_of_particles_variables

    #     labels_to_compare=[rkey['N'],rkey['Nbar'],rkey['time'],rkey['pupx'],rkey['pupy'],rkey['pupz'],rkey['pupt'],rkey['pos_x'],rkey['pos_y'],rkey['pos_z']]

    #     identifier_given=[]
    #     identifier_per=[]

    #     for i in range(len(data_given)):
            
    #         given_string=''
    #         per_string=''
            
    #         for j in labels_to_compare:
    #             given_string=given_string+str(data_given[i][j])
    #             per_string=per_string+str(data_per[i][j])

    #         identifier_given.append(given_string)
    #         identifier_per.append(per_string)

    #     identifier_given=np.array(identifier_given)
    #     identifier_per=np.array(identifier_per)

    #     index_given=np.argsort(identifier_given)
    #     index_per=np.argsort(identifier_per)

    #     assert np.all(identifier_given[index_given]==identifier_per[index_per]),'\n \n ---> ending execution: particles do not match \n \n'

    #     data_given=data_given[index_given]
    #     data_per=data_per[index_per]

    #     del identifier_given
    #     del identifier_per
    #     del index_given
    #     del index_per
    #     del labels_to_compare

    #     keys=['f00_Re', 'f01_Re', 'f01_Im', 'f02_Re', 'f02_Im', 'f11_Re', 'f12_Re', 'f12_Im' ,'f22_Re', 'f00_Rebar', 'f01_Rebar', 'f01_Imbar', 'f02_Rebar', 'f02_Imbar', 'f11_Rebar', 'f12_Rebar' ,'f12_Imbar', 'f22_Rebar']

    #     ssmag_per=0
    #     ssmag_ori=0
    #     ssdiff=0

    #     for key in keys:
            
    #         ssmag_per=ssmag_per+np.sum(np.square(data_per[:,rkey[key]]))
    #         ssmag_ori=ssmag_ori+np.sum(np.square(data_given[:,rkey[key]]))
    #         ssdiff=ssdiff+np.sum(np.square(data_per[:,rkey[key]]-data_given[:,rkey[key]]))

    #     hf = h5py.File(str(dire)+".h5", 'w')

    #     hf.create_dataset('time', data=t)
    #     hf.create_dataset('single_particle_given', data=data_given[particle_index])
    #     hf.create_dataset('single_particle_per', data=data_per[particle_index])
    #     hf.create_dataset('difference_state_space_vector_magnitud', data=np.sqrt(ssdiff))
    #     hf.create_dataset('state_space_vector_magnitud_per', data=np.sqrt(ssmag_per))
    #     hf.create_dataset('state_space_vector_magnitud_given', data=np.sqrt(ssmag_ori))

    #     hf.close()

    #     del data_per
    #     del data_given
    #     del ssmag_per
    #     del ssmag_ori
    #     del ssdiff

    return dire

# run the write hdf5 files function in parallel
if __name__ == '__main__':
    pool = Pool(nproc)
    finalresult=pool.map(writehdf5files,directories)
    for i in finalresult: print("completed ---> "+i)
