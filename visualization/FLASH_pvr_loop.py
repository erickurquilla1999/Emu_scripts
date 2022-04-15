import os
import glob

output_base = "/project/projectdirs/m3761/FLASH/FFI_3D/NSM/sim1/from_payne/v0/NSM_sim_hdf5_chk_"
store_base = "/project/projectdirs/m3761/FLASH/FFI_3D/NSM/sim1/from_payne/v0/volume_rendering/"

directories = sorted(glob.glob(output_base+"*"))
start_ind = 0
end_ind = 100

for d in directories[start_ind:end_ind]:
    os.system('python3 ~/Emu_scripts/visualization/FLASH_phase_volume_render.py {name} -lo 0 0 0 -hi 7.86524303 7.86524303 7.86524303 -f all -o {storage}'.format(name=d, storage=store_base))

