"""
Download a synoptic HMI magnetogram and calculate the squashing factor Q by calling FastQSL directly from Python
This demo is modified from HMI_Example.py of https://github.com/Valentin-Aslanyan/UFiT
"""
import pfsspy, sunpy, wget, os
import astropy.units
import numpy as np
from fastqsl import fastqsl
# ------------------------------------------------------------
def get_HMI_magnetogram(cr_target,data_dir):
	if cr_target<2096:
		print("CR "+str(cr_target)+" too early for HMI")
		return ""
	target_filename='hmi.Synoptic_Mr.'+str(cr_target)+'.fits'
	
	if target_filename not in os.listdir(data_dir):
		wget.download('http://jsoc.stanford.edu/data/hmi/synoptic/'+target_filename, out=data_dir)
	if target_filename in os.listdir(data_dir):
		return os.path.join(data_dir,target_filename)
	else:
		print("Could not get HMI map for CR "+str(cr_target))
		return ""
# ------------------------------------------------------------
target_CR=2284	#Carrington rotation

num_r=60	#dimensions of B grid
num_t=180
num_p=360

Rss=2.5		#Source surface
cdir = os.getcwd()+os.sep 
HMI_data_dir=cdir
# ------------------------------------------------------------
#Download magnetogram fits file
HMI_filename=get_HMI_magnetogram(target_CR, HMI_data_dir)
HMI_map = sunpy.map.Map(HMI_filename)

#Downsample and remove NaNs as required by pfsspy
HMI_map = HMI_map.resample([num_p, num_t] * astropy.units.pix)
HMI_map.data[np.isnan(HMI_map.data)]=0.0	#NaNs set to zero

#Compute B from PFSS and transform the field and coordinates into required form
pfss_in = pfsspy.Input(HMI_map, num_r, Rss)
pfss_out = pfsspy.pfss(pfss_in)
# ------------------------------------------------------------
# compute Q with FastQSL

b_r  = pfss_out.bg[:,:,:,2].transpose(2,1,0)
b_lat=-pfss_out.bg[:,:,:,1].transpose(2,1,0)
b_lon= pfss_out.bg[:,:,:,0].transpose(2,1,0)

radius=np.exp(pfss_out.grid.rg)
lat_rad=np.arcsin(pfss_out.grid.sg)
lon_rad=pfss_out.grid.pg

fastqsl(b_lon, b_lat, b_r, \
xa=lon_rad, ya=lat_rad, za=radius, spherical=True, \
preview=True)

r_cut=2
fastqsl(b_lon[r_cut:,:,:], b_lat[r_cut:,:,:], b_r[r_cut:,:,:], \
xa=lon_rad, ya=lat_rad, za=radius[r_cut:], spherical=True, \
fname='r_cut2', preview=True)
