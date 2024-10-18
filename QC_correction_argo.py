### Last updated January 11, 2024 by Annabel Wade ###

### External imports ###
import xarray as xr, os, numpy as np
from os.path import join
from QC_functions import *
main_dir = os.getcwd()  # Main / current directory
os.chdir(main_dir)

### Program parameters ###
flag_QC = True 					# Flag bad points and profiles
plot_QC = False 				# Plot raw vs. QC data
save_QC = True	 				# Save flag files 
overwrite_flag_files = False	# If True, overwrites the flag array files

### Directories and file management ###
subdirs = ["data", "figures", "flag_array_files"]

# Make any subdirectories that don't already exist
for subdir in subdirs:
    make_dir(join(main_dir, subdir))
data_dir, figures_dir, flag_array_dir = [join(main_dir, subdir) for subdir in subdirs]

maxnum = 80
# Loop through each data file
for filenum in range(1,maxnum+1):

	# Get the file names
	data_nc = join( data_dir, 'global_profiles.{:02d}.nc'.format(filenum) )
	flag_arrays_nc = join( flag_array_dir, 'flag_arrays.{:02d}.nc'.format(filenum) )
	print(data_nc.split('\\')[-1])

	# Create copies of the raw data and flag_arrays array file
	data = xr.open_dataset(data_nc)
	flag_arrays = xr.open_dataset(flag_arrays_nc)

	ds = data.copy()
	# Implement flag_arrays arrays onto raw data
	ds[['TEMP','PRES','PSAL']] = ds[['TEMP','PRES','PSAL']].where(flag_arrays['bad_pnts']==0)
	ds = ds.where(flag_arrays['bad_prof']<4)

	# Choose variables with sample as a dimension, excluding those with dtype object
	sample_vars = [var for var in ds.data_vars 
				if ('sample' in ds[var].dims 
				and ds[var].dtype != 'object')] 

	# Profile flag_arrays 1 and 3: Out of order or repeated pressure
	# --------------------------------------------------------------------------
	if 1 in flag_arrays['bad_prof'] or 3 in flag_arrays['bad_prof']:
		# get indicies of cycles with flag_arrays 1 or 3
		flag_13_cycs = np.where((flag_arrays['bad_prof'] == 1) | (flag_arrays['bad_prof'] == 3))[0]
		for cyc_i in flag_13_cycs:
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

	# Profile flag_arrays 2 and 3: Salinity gradient larger than 0.002 below 1000 dbar
	# --------------------------------------------------------------------------
	data_vars = ['PRES', 'TEMP', 'PSAL'] # variables to correct
	
	if 2 in flag_arrays['bad_prof'] or 3 in flag_arrays['bad_prof']:
		flag_23_cycs = np.where((flag_arrays['bad_prof'] == 2) | (flag_arrays['bad_prof'] == 3))[0]
		for cyc_i in flag_23_cycs:
			### Correcting for last sample causing salinity gradient ###
			cycle_data = ds.isel(cycle=int(cyc_i))[data_vars]
			last_sample_i = np.where(~np.isnan(cycle_data['PRES']))[0][-1]
			for var in data_vars:
				cycle_data[var].loc[dict(sample=last_sample_i)] = np.nan
			ds[data_vars].loc[dict(cycle=int(cyc_i))] = cycle_data
