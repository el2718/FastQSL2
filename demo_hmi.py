"""
Download a synoptic HMI magnetogram and calculate the squashing factor Q by calling FastQSL directly from Python
This demo is modified from HMI_Example.py of https://github.com/Valentin-Aslanyan/UFiT
"""
import pfsspy, sunpy, wget, os, pickle
import astropy.units
import numpy as np
from fastqsl import fastqsl
# ------------------------------------------------------------
def HMI_pfss4fastqsl(num_CR, num_r, num_t, num_p, Rss, data_dir):
    Bfile='B4fastqsl'+str(num_CR)+'.pkl'
    if Bfile not in os.listdir(data_dir):
        
        if num_CR<2096: raise Exception("CR "+str(num_CR)+" too early for HMI")
        
        HMI_file='hmi.Synoptic_Mr.'+str(num_CR)+'.fits'
        
        if HMI_file not in os.listdir(data_dir):
            wget.download('http://jsoc.stanford.edu/data/hmi/synoptic/'+HMI_file, out=data_dir)
            if HMI_file not in os.listdir(data_dir): 
                raise Exception("Could not get HMI map for CR "+str(num_CR))
            
        HMI_map = sunpy.map.Map(HMI_file)
        
        # Downsample and remove NaNs as required by pfsspy
        HMI_map = HMI_map.resample([num_p, num_t] * astropy.units.pix)
        HMI_map.data[np.isnan(HMI_map.data)]=0.0    # NaNs set to zero

        pfss_in  = pfsspy.Input(HMI_map, num_r, Rss)
        pfss_out = pfsspy.pfss(pfss_in)

        b_lon= pfss_out.bg[:,:,:,0].transpose(2,1,0)
        b_lat=-pfss_out.bg[:,:,:,1].transpose(2,1,0)
        b_r  = pfss_out.bg[:,:,:,2].transpose(2,1,0)

        lon_rad= pfss_out.grid.pg
        lat_rad= np.arcsin(pfss_out.grid.sg)
        radius = np.exp(pfss_out.grid.rg)
    
        with open(Bfile, 'wb') as file: 
            pickle.dump((b_lon, b_lat, b_r, lon_rad, lat_rad, radius), file)
    return Bfile
# ------------------------------------------------------------
num_CR=2281    # Carrington rotation

num_r=60    # dimensions of B grid
num_t=180
num_p=360

Rss=2.5        # the radius of Source surface

data_dir = os.getcwd()+os.sep
# ------------------------------------------------------------
Bfile= HMI_pfss4fastqsl(num_CR, num_r, num_t, num_p, Rss, data_dir)
with open(data_dir+Bfile, "rb") as file:
    b_lon, b_lat, b_r, lon_rad, lat_rad, radius = pickle.load(file)
# ------------------------------------------------------------
# compute Q at bottom
fastqsl(b_lon, b_lat, b_r, \
xa=lon_rad, ya=lat_rad, za=radius, spherical=True, \
fname='pfss_orig', preview=True, keep_tmp=True)

# remove first two layers to remove small scale structure
r_cut=2
fastqsl(b_lon[r_cut:,:,:], b_lat[r_cut:,:,:], b_r[r_cut:,:,:], \
xa=lon_rad, ya=lat_rad, za=radius[r_cut:], spherical=True, \
fname='pfss_rcut2', preview=True, keep_tmp=True)

# trace field lines from two points
# Since keep_tmp=True was set in the command above, bfield.bin has already been saved in tmp_dir; 
# therefore, the input magnetic field is unnecessary here
qsl=fastqsl(\
# b_lo[r_cut:,:,:], b_lat[r_cut:,:,:], b_r[r_cut:,:,:], \
# xa=lon_rad, ya=lat_rad, za=radius[r_cut:], spherical=True, \
fname='pfss_rcut2_seed_path', preview=True, \
length_out=True, \
seed=[[np.pi*0.85, 0.1, 1.7], [np.pi*1.5, -0.2, 1.2]], \
path_out=True, loopB_out=True)