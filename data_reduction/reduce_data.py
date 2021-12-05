# used to make plots but now just generates a hdf5 file with domain-averaged data.
# Run in the directory of the simulation the data should be generated for.
# Still has functionality for per-snapshot plots, but the line is commented out.
# This version averages the magnitudes of off-diagonal components rather than the real/imaginary parts
# also normalizes fluxes by sumtrace of N rather than F.
# This data is used for the growth plot.
# Note - also tried a version using maxima rather than averages, and it did not make the growth plot look any better.

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
NF = 2
nproc = 8
do_average = False
do_fft     = False
do_angular = True

do_MPI = False
output_base = "plt"

#####################
# FFT preliminaries #
#####################
def get_kmid(fft):
    if fft.kx is not None:
        kmid = fft.kx[np.where(fft.kx>=0)]
    if fft.ky is not None:
        kmid = fft.ky[np.where(fft.ky>=0)]
    if fft.kz is not None:
        kmid = fft.kz[np.where(fft.kz>=0)]
    return kmid

def fft_coefficients(fft):
    # add another point to the end of the k grid for interpolation
    # MAKES POWER SPECTRUM HAVE SIZE ONE LARGER THAN KTEMPLATE
    kmid = get_kmid(fft)
    dk = kmid[1]-kmid[0]
    kmid = np.append(kmid, kmid[-1]+dk)
    
    # compute the magnitude of the wavevector for every point
    kmag = 0
    if fft.kx is not None:
        kmag = kmag + fft.kx[:,np.newaxis,np.newaxis]**2
    if fft.ky is not None:
        kmag = kmag + fft.ky[np.newaxis,:,np.newaxis]**2
    if fft.kz is not None:
        kmag = kmag + fft.kz[np.newaxis,np.newaxis,:]**2
    kmag = np.sqrt(np.squeeze(kmag))
    kmag[np.where(kmag>=kmid[-1])] = kmid[-1]
    
 
    # compute left index for interpolation
    ileft = (kmag/dk).astype(int)
    iright = ileft+1
    iright[np.where(iright>=len(kmid)-1)] = len(kmid)-1

    # compute the fraction of the power that goes toward the left and right k point
    cleft = (kmid[iright]-kmag)/dk
    cright = 1.0-cleft

    return cleft, cright, ileft, iright, kmid

def fft_power(fft, cleft, cright, ileft, iright, kmid):

    # compute power contributions to left and right indices
    power = fft.magnitude**2
    powerLeft = power*cleft
    powerRight = power*cright

    # accumulate onto spectrum
    spectrum = np.array( [ 
        np.sum( powerLeft*(ileft ==i) + powerRight*(iright==i) )
        for i in range(len(kmid))] )

    return spectrum

#########################
# average preliminaries #
#########################
def get_matrix(base,suffix,fmt="Emu"):
    if(fmt=="Emu"):
        f00  = ad['boxlib',base+"00_Re"+suffix]
        f01  = ad['boxlib',base+"01_Re"+suffix]
        f01I = ad['boxlib',base+"01_Im"+suffix]
        f11  = ad['boxlib',base+"11_Re"+suffix]
        if(NF>=3):
            f02  = ad['boxlib',base+"02_Re"+suffix]
            f02I = ad['boxlib',base+"02_Im"+suffix]
            f12  = ad['boxlib',base+"12_Re"+suffix]
            f12I = ad['boxlib',base+"12_Im"+suffix]
            f22  = ad['boxlib',base+"22_Re"+suffix]

    if(fmt=="FLASH"):
        assert NF==2
        # need to translate Emu dataset names to FLASH ones
        # WARNING - we are calculating energy densities instead of number densities
        energyGroup = "01"

        if base=="N":
            baseFlash = "e"
        if base=="Fx":
            baseFlash = "f"
        if base=="Fy":
            baseFlash = "g"
        if base=="Fz":
            baseFlash = "h"
        
        if suffix=="":
            suffixFlash = ["e","m","r","i"]
        if suffix=="bar":
            suffixFlash = ["a","n","s","j"]

        f00  = ad['flash',baseFlash+suffixFlash[0]+energyGroup]
        f11  = ad['flash',baseFlash+suffixFlash[1]+energyGroup]
        f01  = ad['flash',baseFlash+suffixFlash[2]+energyGroup]
        f01I = ad['flash',baseFlash+suffixFlash[3]+energyGroup]
            
    zero = np.zeros(np.shape(f00))

    if(NF==2):
        fR = [[f00 , f01 ], [ f01 ,f11 ]]
        fI = [[zero, f01I], [-f01I,zero]]
    if(NF==3):
        fR = [[f00 , f01 , f02 ], [ f01 ,f11 ,f12 ], [ f02 , f12 ,f22 ]]
        fI = [[zero, f01I, f02I], [-f01I,zero,f12I], [-f02I,-f12I,zero]]
    return fR, fI

def sumtrace_N(N):
    sumtrace = 0
    for fi in range(NF):
        sumtrace += np.sum(N[fi][fi])
    return sumtrace

def averaged_N(N, NI, sumtrace):
    R=0
    I=1
    
    # do the averaging
    # f1, f2, R/I
    Nout = np.zeros((NF,NF))
    for i in range(NF):
        for j in range(NF):
            Nout[i][j] = float(np.sum(np.sqrt(N[i][j]**2 + NI[i][j]**2)) / sumtrace)
    return np.array(Nout)

def averaged_F(F, FI, sumtrace):
    R=0
    I=1
    
    # do the averaging
    # direction, f1, f2, R/I
    Fout = np.zeros((3,NF,NF))
    for i in range(3):
        for j in range(NF):
            for k in range(NF):
                Fout[i][j][k] = float(np.sum(np.sqrt( F[i][j][k]**2 + FI[i][j][k]**2))/sumtrace)

    return Fout

def offdiagMag(f):
    R = 0
    I = 1
    result = 0
    for f0 in range(NF):
        for f1 in range(f0+1,NF):
            result += f[:,f0,f1,R]**2 + f[:,f0,f1,I]**2
    return np.sqrt(result)



#########################
# angular preliminaries #
#########################
rkey, ikey = amrex.get_3flavor_particle_keys()
paramfile = open("inputs","r")
for line in paramfile:
    line_without_comments = line.split("#")[0]
    if "nphi_equator" in line_without_comments:
        nl = int(line_without_comments.split("=")[1]) // 2
paramfile.close()
nl += 1

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

# get the number of particles per cell
def get_nppc(d):
    eds = emu.EmuDataset(d)
    t = eds.ds.current_time
    ad = eds.ds.all_data()
    grid_data = GridData(ad)
    level = 0
    gridID = 0
    idata, rdata = amrex.read_particle_data(d, ptype="neutrinos", level_gridID=(level,gridID))
    idlist = grid_data.get_particle_cell_ids(rdata)
    ncells = np.max(idlist)+1
    nppc = len(idlist) // ncells
    return nppc

    
# input list of particle data separated into grid cells
# output the same array, but sorted by zenith angle, then azimuthal angle
# also output the grid of directions in each cell (assumed to be the same)
def sort_rdata_chunk(p):
    # sort first in theta
    sorted_indices = p[:,rkey["pupz"]].argsort()
    p = p[sorted_indices,:]

    # loop over unique values of theta
    costheta = p[:,rkey["pupz"]] / p[:,rkey["pupt"]]
    for unique_costheta in np.unique(costheta):
        # get the array of particles with the same costheta
        costheta_locs = np.where(costheta == unique_costheta)[0]
        p_theta = p[costheta_locs,:]
        
        # sort these particles by the azimuthal angle
        phi = np.arctan2(p_theta[:,rkey["pupy"]] , p_theta[:,rkey["pupx"]] )
        sorted_indices = phi.argsort()
        p_theta = p_theta[sorted_indices,:]
        
        # put the sorted data back into p
        p[costheta_locs,:] = p_theta
        
    # return the sorted array
    return p

def Ylm_indices(l):
    start = l**2
    stop = (l+1)**2
    return start,stop

def create_shared_Ylm_star(rdata):
    nparticles = len(rdata)

    # create local arrays that use this memory for easy modification
    Ylm_star = np.frombuffer(Ylm_star_shared, dtype='complex').reshape( ( (nl+1)**2, nparticles) )
    
    # get direction coordinates
    pupx = rdata[:,rkey["pupx"]]
    pupy = rdata[:,rkey["pupy"]]
    pupz = rdata[:,rkey["pupz"]]
    pupt = rdata[:,rkey["pupt"]]
    xhat = pupx/pupt
    yhat = pupy/pupt
    zhat = pupz/pupt            
    theta = np.arccos(zhat)
    phi = np.arctan2(yhat,xhat)
                
    # evaluate spherical harmonic amplitudes
    for l in range(nl):
        start,stop = Ylm_indices(l)
        nm = stop-start
        mlist = np.array(range(nm))-l
        Ylm_star_thisl = [np.conj(scipy.special.sph_harm(m, l, phi, theta)) for m in mlist]
        Ylm_star[start:stop] = np.array( Ylm_star_thisl )

def get_shared_Ylm_star(l,nparticles):
    start, stop = Ylm_indices(l)
    Ylm_star = np.frombuffer(Ylm_star_shared, dtype='complex').reshape( ( (nl+1)**2, nparticles) )
    return Ylm_star[start:stop]

# use scipy.special.sph_harm(m, l, azimuthal_angle, polar_angle)
# np.arctan2(y,x)
def spherical_harmonic_power_spectrum_singlel(l, Nrho):
    nparticles = np.shape(Nrho)[2]
    Ylm_star = get_shared_Ylm_star(l,nparticles)
    Nrholm_integrand = np.array([Nrho*Ylm_star[im,:] for im in range(len(Ylm_star))])
    Nrholm = np.sum(Nrholm_integrand, axis=3)
    result = np.sum(np.abs(Nrholm)**2, axis=0)
    return result

def get_Nrho(p):
    # build Nrho complex values
    nparticles = len(p)
    if NF==2:
        Nrho = np.zeros((2,3,nparticles))*1j
        Nrho[0,0,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f00_Re"   ]] + 1j*0                      )
        Nrho[0,1,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f01_Re"   ]] + 1j*p[:,rkey["f01_Im"   ]] )
        Nrho[0,2,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f11_Re"   ]] + 1j*0                      )
        Nrho[1,0,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f00_Rebar"]] + 1j*0                      )
        Nrho[1,1,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f01_Rebar"]] + 1j*p[:,rkey["f01_Imbar"]] )
        Nrho[1,2,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f11_Rebar"]] + 1j*0                      )
    if NF==3:
        Nrho = np.zeros((2,6,nparticles))*1j
        Nrho[0,0,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f00_Re"   ]] + 1j*0                      )
        Nrho[0,1,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f01_Re"   ]] + 1j*p[:,rkey["f01_Im"   ]] )
        Nrho[0,2,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f02_Re"   ]] + 1j*p[:,rkey["f02_Im"   ]] )
        Nrho[0,3,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f11_Re"   ]] + 1j*0                      )
        Nrho[0,4,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f12_Re"   ]] + 1j*p[:,rkey["f12_Im"   ]] )
        Nrho[0,5,:] = p[:,rkey["N"   ]] * ( p[:,rkey["f22_Re"   ]] + 1j*0                      )
        Nrho[1,0,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f00_Rebar"]] + 1j*0                      )
        Nrho[1,1,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f01_Rebar"]] + 1j*p[:,rkey["f01_Imbar"]] )
        Nrho[1,2,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f02_Rebar"]] + 1j*p[:,rkey["f02_Imbar"]] )
        Nrho[1,3,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f11_Rebar"]] + 1j*0                      )
        Nrho[1,4,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f12_Rebar"]] + 1j*p[:,rkey["f12_Imbar"]] )
        Nrho[1,5,:] = p[:,rkey["Nbar"]] * ( p[:,rkey["f22_Rebar"]] + 1j*0                      )
    return Nrho

def spherical_harmonic_power_spectrum(Nrho):
    spectrum = np.array([spherical_harmonic_power_spectrum_singlel(l, Nrho) for l in range(nl)])
    return spectrum

#########################
# loop over directories #
#########################
if do_MPI:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()
else:
    mpi_rank = 0
    mpi_size = 1
directories = sorted(glob.glob(output_base+"*"))
if( (not do_average) and (not do_fft)):
    directories = []
for d in directories[mpi_rank::mpi_size]:
    print("# rank",mpi_rank,"is working on", d)
    eds = emu.EmuDataset(d)
    t = eds.ds.current_time
    ad = eds.ds.all_data()

    ################
    # average work #
    ################
    # write averaged data
    outputfilename = d+"/reduced_data.h5"
    already_done = len(glob.glob(outputfilename))>0
    if do_average and not already_done:
        thisN, thisNI = get_matrix("N",""   )
        sumtrace = sumtrace_N(thisN)
        trace = sumtrace
        N = averaged_N(thisN,thisNI,sumtrace)

        thisN, thisNI = get_matrix("N","bar")
        sumtrace = sumtrace_N(thisN)
        tracebar = sumtrace
        Nbar = averaged_N(thisN,thisNI,sumtrace)

        thisFx, thisFxI = get_matrix("Fx","")
        thisFy, thisFyI = get_matrix("Fy","")
        thisFz, thisFzI = get_matrix("Fz","")
        Ftmp  = np.array([thisFx , thisFy , thisFz ])
        FtmpI = np.array([thisFxI, thisFyI, thisFzI])
        F = averaged_F(Ftmp, FtmpI,sumtrace)
    
        thisFx, thisFxI = get_matrix("Fx","bar") 
        thisFy, thisFyI = get_matrix("Fy","bar") 
        thisFz, thisFzI = get_matrix("Fz","bar") 
        Ftmp  = np.array([thisFx , thisFy , thisFz ])
        FtmpI = np.array([thisFxI, thisFyI, thisFzI])
        Fbar = averaged_F(Ftmp, FtmpI,sumtrace)

        print("# rank",mpi_rank,"writing",outputfilename)
        avgData = h5py.File(outputfilename,"w")
        avgData["N_avg_mag"] = [N,]
        avgData["Nbar_avg_mag"] = [Nbar,]
        avgData["F_avg_mag"] = [F,]
        avgData["Fbar_avg_mag"] = [Fbar,]
        avgData["t"] = [t,]
        avgData.close()

    ############
    # FFT work #
    ############
    outputfilename = d+"/reduced_data_fft_power.h5"
    already_done = len(glob.glob(outputfilename))>0
    if do_fft and not already_done:

        print("# rank",mpi_rank,"writing",outputfilename)
        fout = h5py.File(outputfilename,"w")
        fout["t"] = [np.array(t),]

        fft = eds.fourier("N00_Re",nproc=nproc)
        fout["k"] = get_kmid(fft)
        cleft, cright, ileft, iright, kmid = fft_coefficients(fft)
        N00_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
        fft = eds.fourier("N11_Re",nproc=nproc)
        N11_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
        fft = eds.fourier("N01_Re","N01_Im",nproc=nproc)
        N01_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
        fout["N00_FFT"] = [np.array(N00_FFT),]
        fout["N11_FFT"] = [np.array(N11_FFT),]
        fout["N01_FFT"] = [np.array(N01_FFT),]
        if NF>2:
            fft = eds.fourier("N22_Re",nproc=nproc)
            N22_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
            fft = eds.fourier("N02_Re","N02_Im",nproc=nproc)
            N02_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
            fft = eds.fourier("N12_Re","N12_Im",nproc=nproc)
            N12_FFT = fft_power(fft, cleft, cright, ileft, iright, kmid)
            fout["N22_FFT"] = [np.array(N22_FFT),]
            fout["N02_FFT"] = [np.array(N02_FFT),]
            fout["N12_FFT"] = [np.array(N12_FFT),]
        
        #fft = eds.fourier("Fx00_Re")
        #Fx00_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fx11_Re")
        #Fx11_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fx01_Re","Fx01_Im")
        #Fx01_FFT.append(fft_power(fft))
        #fout["Fx00_FFT"] = [np.array(Fx00_FFT),]
        #fout["Fx11_FFT"] = [np.array(Fx11_FFT),]
        #fout["Fx01_FFT"] = [np.array(Fx01_FFT),]
        #if NF>2:
            #fft = eds.fourier("Fx22_Re")
            #Fx22_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fx02_Re","Fx02_Im")
            #Fx02_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fx12_Re","Fx12_Im")
            #Fx12_FFT.append(fft_power(fft))
            #fout["Fx22_FFT"] = [np.array(Fx22_FFT),]
            #fout["Fx02_FFT"] = [np.array(Fx02_FFT),]
            #fout["Fx12_FFT"] = [np.array(Fx12_FFT),]
        
        #fft = eds.fourier("Fy00_Re")
        #Fy00_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fy11_Re")
        #Fy11_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fy01_Re","Fy01_Im")
        #Fy01_FFT.append(fft_power(fft))
        #fout["Fy00_FFT"] = [np.array(Fy00_FFT),]
        #fout["Fy11_FFT"] = [np.array(Fy11_FFT),]
        #fout["Fy01_FFT"] = [np.array(Fy01_FFT),]
        #if NF>2:
            #fft = eds.fourier("Fy22_Re")
            #Fy22_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fy02_Re","Fy02_Im")
            #Fy02_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fy12_Re","Fy12_Im")
            #Fy12_FFT.append(fft_power(fft))
            #fout["Fy22_FFT"] = [np.array(Fy22_FFT),]
            #fout["Fy02_FFT"] = [np.array(Fy02_FFT),]
            #fout["Fy12_FFT"] = [np.array(Fy12_FFT),]
        
        #fft = eds.fourier("Fz00_Re")
        #Fz00_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fz11_Re")
        #Fz11_FFT.append(fft_power(fft))
        #fft = eds.fourier("Fz01_Re","Fz01_Im")
        #Fz01_FFT.append(fft_power(fft))
        #fout["Fz00_FFT"] = [np.array(Fz00_FFT),]
        #fout["Fz11_FFT"] = [np.array(Fz11_FFT),]
        #fout["Fz01_FFT"] = [np.array(Fz01_FFT),]
        #if NF>2:
            #fft = eds.fourier("Fz22_Re")
            #Fz22_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fz02_Re","Fz02_Im")
            #Fz02_FFT.append(fft_power(fft))
            #fft = eds.fourier("Fz12_Re","Fz12_Im")
            #Fz12_FFT.append(fft_power(fft))
            #fout["Fz22_FFT"] = [np.array(Fz22_FFT),]
            #fout["Fz02_FFT"] = [np.array(Fz02_FFT),]
            #fout["Fz12_FFT"] = [np.array(Fz12_FFT),]

        fout.close()

# separate loop for angular spectra so there is no aliasing and better load balancing
directories = sorted(glob.glob("plt*/neutrinos"))
directories = [directories[i].split('/')[0] for i in range(len(directories))] # remove "neutrinos"

# get number of particles to be able to construct 
nppc = get_nppc(directories[-1])

# create shared object
# double the size for real+imaginary parts
Ylm_star_shared = mp.RawArray('d', int(2 * (nl+1)**2 * nppc))


if __name__ == '__main__':
    pool = Pool(nproc)
    for d in directories:
        if mpi_rank==0:
            print("# working on", d)
        eds = emu.EmuDataset(d)
        t = eds.ds.current_time
        ad = eds.ds.all_data()

        ################
        # angular work #
        ################
        outputfilename = d+"/reduced_data_angular_power_spectrum.h5"
        already_done = len(glob.glob(outputfilename))>0
        if do_angular and not already_done:

            if mpi_rank==0:
                print("Computing up to l =",nl-1)

            header = amrex.AMReXParticleHeader(d+"/neutrinos/Header")
            grid_data = GridData(ad)
            nlevels = len(header.grids)
            assert nlevels==1
            level = 0
            ngrids = len(header.grids[level])
        
            # average the angular power spectrum over many cells
            # loop over all cells within each grid
            if NF==2:
                ncomps = 3
            if NF==3:
                ncomps = 6
            spectrum = np.zeros((nl,2,ncomps))
            Nrho_avg = np.zeros((2,ncomps,nppc))*1j
            total_ncells = 0
            for gridID in range(mpi_rank,ngrids,mpi_size):
                print("    rank",mpi_rank,"grid",gridID+1,"/",ngrids)
            
                # read particle data on a single grid
                idata, rdata = amrex.read_particle_data(d, ptype="neutrinos", level_gridID=(level,gridID))
                
                # get list of cell ids
                idlist = grid_data.get_particle_cell_ids(rdata)
                
                # sort rdata based on id list
                sorted_indices = idlist.argsort()
                rdata = rdata[sorted_indices]
                idlist = idlist[sorted_indices]
                
                # split up the data into cell chunks
                ncells = np.max(idlist)+1
                nppc = len(idlist) // ncells
                rdata  = [ rdata[icell*nppc:(icell+1)*nppc,:] for icell in range(ncells)]
                chunksize = ncells//nproc
                if ncells % nproc != 0:
                    chunksize += 1
                
                # sort particles in each chunk
                rdata = pool.map(sort_rdata_chunk, rdata, chunksize=chunksize)
                
                # create the spherical harmonic grid
                if gridID == mpi_rank:
                    create_shared_Ylm_star(rdata[0])
                    phat = rdata[0][:,rkey["pupx"]:rkey["pupz"]+1] / rdata[0][:,rkey["pupt"]][:,np.newaxis]

                # accumulate the spatial average of the angular distribution
                Nrho = pool.map(get_Nrho,rdata, chunksize=chunksize)
                Nrho_avg += np.sum(Nrho, axis=0)
                
                # accumulate a spectrum from each cell
                spectrum_each_cell = pool.map(spherical_harmonic_power_spectrum, Nrho, chunksize=chunksize)
                spectrum += np.sum(spectrum_each_cell, axis=0)

                # count the total number of cells
                total_ncells += ncells
            
            if do_MPI:
                comm.Barrier()
                spectrum     = comm.reduce(spectrum    , op=MPI.SUM, root=0)
                Nrho_avg     = comm.reduce(Nrho_avg    , op=MPI.SUM, root=0)
                total_ncells = comm.reduce(total_ncells, op=MPI.SUM, root=0)

            # write averaged data
            if mpi_rank==0:
                spectrum /= total_ncells*ad['index',"cell_volume"][0]
                Nrho_avg /= total_ncells*ad['index',"cell_volume"][0]
            
                print("# writing",outputfilename)
                avgData = h5py.File(outputfilename,"w")
                avgData["angular_spectrum"] = [spectrum,]
                avgData["Nrho (1|ccm)"] = [Nrho_avg,]
                avgData["phat"] = phat
                avgData["t"] = [t,]
                avgData.close()

    
