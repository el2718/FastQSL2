"""
Download a synoptic HMI magnetogram and calculate the squashing factor Q by calling FastQSL directly from Python
This demo is modified from HMI_Example.py of https://github.com/Valentin-Aslanyan/UFiT
"""
import pfsspy, sunpy, wget, os, pickle
import astropy.units
import numpy as np
from fastqsl import fastqsl
# ------------------------------------------------------------
def pfss4fastqsl(num_CR, *, \
                 Rss=2.5, n_lon=361, n_lat=181, n_r=76, \
                 data_type="hmi.Synoptic_Mr_polfil", \
                 data_dir=os.getcwd()+os.sep, fname=None, \
                 ):
    if fname is None: fname = 'pfss_'+data_type+'.'+str(num_CR)

    Bfile = data_dir+fname+'.pkl'
    if fname+'.pkl' not in os.listdir(data_dir):
        
        if num_CR<2096: raise Exception("CR "+str(num_CR)+" too early for HMI")
        HMI_file=data_type+'.'+str(num_CR)+'.fits'
        if HMI_file not in os.listdir(data_dir):
            wget.download('http://jsoc.stanford.edu/data/hmi/synoptic/'+HMI_file, out=data_dir)
            if HMI_file not in os.listdir(data_dir): 
                raise Exception("Could not get HMI map for CR "+str(num_CR))

        HMI_map = sunpy.map.Map(HMI_file)
        print(dir(HMI_map))
        # Downsample and remove NaNs as required by pfsspy
        HMI_map = HMI_map.resample([n_lon-1, n_lat-1] * astropy.units.pix)
        HMI_map.data[np.isnan(HMI_map.data)]=0.0    # NaNs set to zero
        pfss_in  = pfsspy.Input(HMI_map, n_r-1, Rss)
        pfss_out = pfsspy.pfss(pfss_in)



        # exchange the index order of R and phi (longitude)
        Bvec= pfss_out.bg.transpose(2,1,0,3) 

        # # b_lat = -b_theta
        Bvec[:,:,:,1]= - Bvec[:,:,:,1]

        lon_rad= pfss_out.grid.pg
        lat_rad= np.arcsin(pfss_out.grid.sg)
        radius = np.exp(pfss_out.grid.rg)
    
        with open(Bfile, 'wb') as file: 
            pickle.dump((Bvec, lon_rad, lat_rad, radius), file)

        print(Bfile+' is saved')
    return Bfile
# ------------------------------------------------------------
# Carrington rotation
num_CR= 2284

data_type="hmi.Synoptic_Mr_polfil" 
# data_type="hmi.Synoptic_Mr" 
# data_type="hmi.Synoptic_Ml"
# data_type="hmi.Synoptic_Ml_small"
# data_type="hmi.B_synoptic_Br"
# data_type="hmi.B_synoptic_Br_small"

Bfile= pfss4fastqsl(num_CR, data_type=data_type)
with open(Bfile, "rb") as file:
    Bvec, lon_rad, lat_rad, radius = pickle.load(file)
# ------------------------------------------------------------
r_cut=2
fname = 'pfss_'+data_type+'.'+str(num_CR)

# # # compute Q at bottom
fastqsl(Bvec, xa=lon_rad, ya=lat_rad, za=radius, spherical=True, \
fname= fname+'_orig', preview=True, keep_tmp=True)

# # # remove first two layers to remove small scale structure
fastqsl(Bvec[r_cut:,:,:,:], xa=lon_rad, ya=lat_rad, za=radius[r_cut:], spherical=True, \
fname= fname+'_rcut2', scottFlag=False, preview=True, keep_tmp=True)

# # trace field lines from two points
# # Since keep_tmp=True was set in the command above, bfield.bin has already been saved in tmp_dir; 
# # therefore, the input magnetic field is unnecessary here
qsl=fastqsl(\
# Bvec[r_cut:,:,:,:], xa=lon_rad, ya=lat_rad, za=radius[r_cut:], spherical=True, \
fname= fname+'_rcut2_seed_path', preview=True, \
length_out=True, \
seed=[[np.pi*0.85, 0.1, 1.7], [np.pi*1.5, -0.2, 1.2]], \
path_out=True, loopB_out=True)

# # git clone https://github.com/el2718/par2solarwind
# from par2solarwind import par2solarwind
# # compute two parameters for solar wind modeling at bottom
# par2solarwind(Bvec[r_cut:,:,:, :], lon_rad, lat_rad, radius[r_cut:], \
#               bottomFlag=True, fname=fname+'_rcut2_solarwind', preview=True)
