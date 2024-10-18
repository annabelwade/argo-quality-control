import os
import shutil
import warnings
import re
import numpy as np
import xarray as xr
from tqdm import tqdm
import gsw
from scipy.interpolate import PchipInterpolator
from scipy import integrate
import dask
from datetime  import datetime

from QC_functions import find_nc_filenames,make_dir
warnings.filterwarnings("ignore")


@dask.delayed
def compute_dh_level(interp_func,upper_pres,lower_pres):
	"""
	Creates delayed dask item to integrate the specific volume anomaly interpolated 
	function between two levels.

	Input:
		interp_func: interpolated specific volume function, scipy object
		upper_pres: upper bound of integration
		lower_pres: lower bound of integration

	Returns:
		dynamic height for a given level
	"""
	dh_lev,err = integrate.quad(interp_func,upper_pres,lower_pres)
	return dh_lev

# Input paramters: 
# Change this to True if you want to reload everything - diagnostic use only
# Change this to False if you want to only update missing files
reload_data = False

# Create the levels - these can be changed as needed
level_ar = [0,5,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,
           900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200]

# Select a threshold around the levels
pres_th = 0.1

# This is the directory for everything to get DH
# You should be working in this directory
main_dir = os.getcwd()

# Input directories
data_dir = main_dir + '/NETCDF/'
flag_array_dir = main_dir + '/flag_array_files/'

# Output directory
dynamic_height_dir = main_dir +  '/dynamic_height/'
make_dir(dynamic_height_dir,clean=False)

# Lists of nc files for all directories
datalist = sorted(find_nc_filenames(data_dir))
flaglist = sorted(find_nc_filenames(flag_array_dir))
dhlist = sorted(find_nc_filenames(dynamic_height_dir))

# Check to keep files in the output directory
if reload_data:
	datanum = re.findall(r'\d+', ''.join(datalist))
	load_data = datanum
else:
	datanum = re.findall(r'\d+', ''.join(datalist))
	dhnum = re.findall(r'\d+', ''.join(dhlist))
	load_data = datanum.copy()
	for index in range(0,len(dhnum)):
		if (dhnum[index] in datanum) and (len(datalist) == len(flaglist)) :
			load_data.remove(dhnum[index])

print('Loading files:',load_data)
input('Press enter to proceed...')

for indie in range(0,len(load_data)):

	filenum = int(load_data[indie]) 
	nc_d = data_dir+'global_profiles.{:02d}.nc'.format(filenum)
	nc_f = flag_array_dir+'flag_arrays.{:02d}.nc'.format(filenum)

	# Open the data and flags
	ds = xr.open_dataset(nc_d)
	flag = xr.open_dataset(nc_f)

	print('\n',nc_d)
	print(nc_f)

	# Select only data that is good or can be fixed
	ds['TEMP'] = ds['TEMP'].where(flag['bad_pnts']==0)
	ds['PRES'] = ds['PRES'].where(flag['bad_pnts']==0)
	ds['PSAL'] = ds['PSAL'].where(flag['bad_pnts']==0)
	ds = ds.where(flag['bad_prof']<4,drop=True)
	flag=flag.where(flag['bad_prof']<4,drop=True)

	# Choose variables with sample as a dimension, excluding those with dtype object
	sample_vars = [var for var in ds.data_vars 
				if ('sample' in ds[var].dims 
				and ds[var].dtype != 'object')] 

	# Profile flag 1 and 3: Out of order or repeated pressure
	# --------------------------------------------------------------------------
	if 1 in flag['bad_prof'] or 3 in flag['bad_prof']:
		flag_13_cycs = np.where((flag['bad_prof'] == 1) | (flag['bad_prof'] == 3))
		for cyc_i in flag_13_cycs[0]:
			### Correcting for out of order PRES ###
			cycle_data = ds.isel(cycle=int(cyc_i)).sortby('PRES')
			newdP = cycle_data['PRES'].dropna(dim='sample').diff(dim='sample')
			ds.loc[dict(cycle=int(cyc_i))] = cycle_data
			### Correcting for repeated PRES ###
			if np.sum(newdP == 0)>0:
				sample_cycle_data = cycle_data[sample_vars]
				# Average to have unique pressures
				averaged_data = sample_cycle_data.groupby('PRES').mean()
				averaged_data = averaged_data.swap_dims({'PRES':'sample'})
				averaged_data = averaged_data.reset_coords('PRES')
				# Pad data with NaNs to match original sample length
				pad_after = len(cycle_data['PRES']) - len(averaged_data['PRES'])
				averaged_data = averaged_data.pad(pad_width={'sample': (0, pad_after)})
				ds[sample_vars].loc[dict(cycle=int(cyc_i))] = averaged_data


	# Profile flag 2 and 3: Salinity gradient larger than 0.002 below 1000 dbar
	# --------------------------------------------------------------------------
	if 2 in flag['bad_prof'] or 3 in flag['bad_prof']:
		flag_23_cycs = np.where((flag['bad_prof'] == 2) | (flag['bad_prof'] == 3))
		for cyc_i in flag_23_cycs[0]:
			### Correcting for last sample causing salinity gradient ###
			data_vars = ['PRES', 'TEMP', 'PSAL'] # variables to correct
			cycle_data = ds.isel(cycle=int(cyc_i))[data_vars]
			last_sample_i = np.where(~np.isnan(cycle_data['PRES']))[0][-1]
			for var in data_vars:
				cycle_data[var].loc[dict(sample=last_sample_i)] = np.nan
			ds[data_vars].loc[dict(cycle=int(cyc_i))] = cycle_data

	# Find where there are pressure measurements in each level bin, with a threshold around the levels
	val_mask = np.zeros((len(level_ar)-1,len(ds['cycle'])))
	for levs in range(0,len(level_ar)-1):
		val_mask[levs,:] = sum((ds['PRES']>= (level_ar[levs]-pres_th)) & (ds['PRES']<= (level_ar[levs+1]+pres_th))) > 0

	# Create a mask of levels where pressures are in one bin above and two bins below
	level_maskarray = (val_mask[0:-2,:] + val_mask[1:-1,:] + val_mask[2:,:]) == 3

	# Put the data into an xarray dataset for easy export
	level_data = xr.Dataset(data_vars=dict(
							level_pres=(['level'], level_ar[1:-2]),
							level_mask=(['level','cycle'], level_maskarray)))

	level_data['level_pres'].attrs = {'description': 'level pressure','units':'dbar'}

	# Get the specific volume anomaly from the TEOS-10 package
	sa = gsw.SA_from_SP(ds['PSAL'].transpose(),ds['PRES'].transpose(),ds['LON'],ds['LAT'])
	ct = gsw.CT_from_t(sa,ds['TEMP'].transpose(),ds['PRES'].transpose())
	rho = gsw.sigma0(sa,ct)
	svan = gsw.specvol_anom_standard(sa,ct,ds['PRES'].transpose())

	dh_all = np.ones((len(level_data['level_pres']),len(level_data['cycle'])))*np.nan
	dh_rho = np.ones((len(level_data['level_pres']),len(level_data['cycle'])))*np.nan

	# Go through each cycle to interpolate the specific volume anomaly
	for cycnum in tqdm(range(0,len(level_data['cycle']))):
		dh_future = []
		ind_vals = [] 
		cyc_pres = (ds['PRES']).isel(cycle=cycnum).dropna(dim='sample')
		cyc_mask = level_data['level_mask'].isel(cycle=cycnum)
		interp_func = PchipInterpolator(cyc_pres*1e4,svan.isel(cycle=cycnum).dropna(dim='sample'),extrapolate=False)
		interp_rho = PchipInterpolator(cyc_pres,rho.isel(cycle=cycnum).dropna(dim='sample'),extrapolate=False)
		# Integrate the specific volume anomaly for every viable level using Dask
		for int_num in range(0,len(level_data['level_pres'])):
			if cyc_mask.isel(level=int_num):
				dh_lev = compute_dh_level(interp_func,level_ar[int_num+1]*1e4,level_ar[int_num+2]*1e4)
				dh_future.append(dh_lev)
				ind_vals.append(int_num)
		dh_comp = np.array(dask.compute(*dh_future))
		dh_all[ind_vals,cycnum] = dh_comp
		dh_rho[ind_vals,cycnum] = interp_rho(level_data['level_pres'].isel(level=ind_vals))

	# Compile DataArrays
	dh_array = xr.DataArray(dh_all,dims=['level','cycle'],
		attrs={'description':'dynamic height','units':'m^2/s^2','FillValue':'NaN'})

	dh_dense = xr.DataArray(dh_rho,dims=['level','cycle'],
		attrs={'description':'density','units':'kg/m^2','FillValue':'NaN'})
	# # Take only the levels that are interpolated
	# dh_array = dh_array.where(level_data['level_mask'])

	dh_ds = xr.Dataset(data_vars={'DH':dh_array,
								 'RHO':dh_dense,
	                   			 'LEVEL_PRES':level_data['level_pres'],
	                             'LAT':ds['LAT'],
	                             'LON':ds['LON'],
	                             'DATENUM':ds['DATENUM'],
	                             'FLOATNUM':ds['FLOATNUM'],
	                             'CYCLENUM':ds['CYCLENUM'],
	                             'SPECIAL_AREAS':flag['special_areas']
	                             },
						attrs={'creators':'Katy Christensen & Annabel Wade',
								'date_created':datetime.now().strftime("%Y-%m-%d"),
								'description':'Dynamic height computed from quality controlled '}) 

	# Drop the cycles that are only nans
	# This happens if there isn't adequate vertical resolution
	dh_ds = dh_ds.where(dh_ds['DH'].sum(dim='level')!=0,drop=True)

	# Save .nc file
	dh_ds.to_netcdf(dynamic_height_dir+'dynamic_height.{:02d}.nc'.format(filenum))

# Save all of 
dhlist = sorted(find_nc_filenames(dynamic_height_dir))
alldh = xr.open_mfdataset(dhlist,combine='nested',concat_dim='cycle')
alldh.to_netcdf(main_dir+'/global_DH_'+datetime.strftime(datetime.today(),format='%Y%m%d')+'.nc')



# alldh_loaded = xr.open_dataset(main_dir+'global_DH_20231005.nc')
