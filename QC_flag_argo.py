### Last updated April 1, 2024 by Annabel Wade ###

### External imports ###
import xarray as xr, os, numpy as np, pandas as pd, re, matplotlib.pyplot as plt
from datetime import datetime; from tqdm import trange; from QC_functions import *
from shapely import contains_xy
from os.path import join
main_dir = os.getcwd()  # Main / current directory
if not os.path.isfile(join(main_dir, 'QC_hand_pick_pnts.py')):
	raise FileNotFoundError("QC_hand_pick_pnts.py is not in the current directory.")
from QC_hand_pick_pnts import hp_TEMP, hp_PSAL

### Program parameters ###
flag_QC = True 					# Flag bad points and profiles
plot_QC = False 				# Plot raw vs. QC data
save_QC = True	 				# Save flag files 
overwrite_flag_files = False	# If True, overwrites the flag array files

### Directories and file management ###
subdirs = ["data_dir", "figures_dir", "flag_array_dir"]

# Make any subdirectories that don't already exist
for subdir in subdirs:
    make_dir(join(main_dir, subdir))
data_dir, figures_dir, flag_array_dir = [join(main_dir, subdir) for subdir in subdirs]

# Get list of data files and flag array files to run through program
datalist = find_nc_filenames(data_dir)
if (len(datalist) == 0):
	raise FileNotFoundError("There are no data files present in the data directory.")
if overwrite_flag_files:
	flaglist = []
else: 
	flaglist = find_nc_filenames(flag_array_dir) 
datanums = re.findall(r'\d+', ''.join(datalist)) # Extract file numbers from data filepaths
flagnums = re.findall(r'\d+', ''.join(flaglist)) # Extract file numbers from flag array filepaths

# Overwrite filenums that don't yet have a corresponding flag array file
load_data = [datanums[index] for index in range(len(datanums)) if datanums[index] not in flagnums]

# Load special regions polygons
if not os.path.isfile(join(main_dir,'special_areas_polygons.shp')):
	raise FileNotFoundError("The polygons file does not exist in the current directory.")
polygons = load_polygons(join(main_dir,'special_areas_polygons.shp'))

### Iterate through each data file ###
for index in trange(len(load_data)):

	# Get the file name
	filenum = int(load_data[index])
	nc = join(data_dir, 'global_profiles.{:02d}.nc'.format(filenum))
	print("\n"+nc.split('\\')[-1])

	# Read the raw data and create a copy for manipulating
	data = xr.open_dataset(nc)

	### Quality control flagging of bad points and bad profiles ###
	if flag_QC:
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# QC Point Flags
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# bad_pnts:
		# 0 - Good Data
		# 1 - Missing or out of range pressure 
		# 2 - Missing or out of range temperature
		# 3 - Missing or out of range salinity
	 	# 4 - Bad temperature points hand picked by Annabel Wade 
	 	# 5 - Bad salinity points hand picked by Annabel Wade

		# Create points flag array
		bad_pnts = np.zeros([len(data['sample']), len(data['cycle'])])

		# bad_pnts flag 1-3: missing and out of range data 
		variables = ['PRES', 'TEMP', 'PSAL']
		for v_name in reversed(variables):
			vars()[v_name] = data[v_name]
			if v_name == 'PSAL':  # Consider PSAL < 10 as OOR
				bad_pnts[np.where(data['PSAL'] < 10)] = 3
			curr_flag = np.logical_or(np.isnan(vars()[v_name]).values, 		# missing data
							~np.isnan(data[v_name+'_OUTOFRANGE'].values))   # out of range data
			bad_pnts[curr_flag] = variables.index(v_name) + 1

		# bad_pnts flags 4-5: hand picking bad points
		# Loop through array of index tuples (sample, cyc) in hand picked dictionaries
		if filenum in hp_TEMP.keys():
			for sample_cyc in hp_TEMP[filenum]:
				bad_pnts[sample_cyc[0], sample_cyc[1]] = 4
		if filenum in hp_PSAL.keys():
			for sample_cyc in hp_PSAL[filenum]:
				bad_pnts[sample_cyc[0], sample_cyc[1]] = 5


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# QC Profile Flags
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# bad_prof:
		# 0 - Good Data
		# 1 - Repeated or out of order pressure
		# 2 - Salinity gradient (larger than 0.002 below 1000 dbar)
		# 3 - Both flag 1 and 2
		# 4 - 4 or fewer measurements
		# 5 - Out of bounds time range (everything before 2022) or no time data
		# 6 - No geo-coordinates
		#
		# special_areas:
		# 0 - Non-special
		# 1 - Mediterranean Sea
		# 2 - Baltic Sea
		# 3 - Indian Ocean flow-through; coastal West Pacific Ocean
		# 4 - Gulf of Mexico; Caribbean Sea
		# 5 - Coastal Pakistan

		# Create profile flag arrays
		bad_prof = np.zeros(len(data['cycle']))
		special_loc = np.zeros(bad_prof.shape)

		# Select just the pressure and salinity data that is not pre-flagged
		PRES = PRES.where(bad_pnts==0)
		PSAL = PSAL.where(bad_pnts==0)
		TEMP = TEMP.where(bad_pnts==0)

		# Get the number of measurements in each profile
		cyc_len = np.logical_not(np.isnan(PRES)).sum(axis=0)
	
		# Get the P and S gradients below 1000 dbar
		dP_1000 = PRES.where(PRES>1000).diff('sample')
		dS_1000 = PSAL.where(PRES>1000).diff('sample')
		meandP_1000 = dP_1000.mean(axis=0,skipna=True)
		sgrad = abs(dS_1000/dP_1000)
		
		for cyc_i in range(0,len(cyc_len)):
			# bad_prof flags 1-3: pressure order and salinity gradient is larger than 0.002
			# Only flag profiles where the last dP is < 10 and mean dP > 20
			dP = PRES.sel(cycle=cyc_i).dropna(dim='sample').diff(dim='sample')
			if sum(dP <= 0) > 0:
				bad_prof[cyc_i] = 1
			cyc = dP_1000.sel(cycle=cyc_i)
			cyc = cyc.where(~np.isnan(cyc), drop=True)
			if len(cyc) > 0:
				if cyc[-1] < 10 and meandP_1000[cyc_i] > 20:
					grade = sgrad.sel(cycle=cyc_i).values
					grade = grade[~np.isnan(grade)]
					if grade[-1] > 0.002 and np.nanmean(grade) < 0.002:
						if bad_prof[cyc_i] == 1:  # both flag 1 and 2
							bad_prof[cyc_i] = 3
						else:
							bad_prof[cyc_i] = 2

			# Special areas flag 1-5
			lat = data['LAT'].sel(cycle=cyc_i).values
			lon = data['LON'].sel(cycle=cyc_i).values
			for polygon_i in range(len(polygons)):
				if contains_xy(polygons[polygon_i], lon, lat):
					special_loc[cyc_i] = polygon_i+1
					continue

		# bad_prof flag 4: 4 or fewer measurements
		bad_prof[cyc_len.values < 5] = 4

		# bad_prof flag 5: date isn't before 2022 or missing time data
		datenums = data['DATENUM'].values
		timestamps = pd.to_datetime(datenums-719529, unit='D')
		bad_prof[np.logical_or(timestamps>datetime(2022,1,1,0,0,0,0),np.isnan(timestamps))] = 5

		# bad_prof flag 6: locations that were flagged bad or missing
		lat = data['LAT']
		lon = data['LON']
		PQ = data['POS_QC']
		bad_prof[np.where(PQ.values > 2)[0]] = 6
		bad_prof[np.logical_or(np.isnan(lat.values),np.isnan(lon.values))] = 6
		

	### Plot raw vs. QC data ###
	if plot_QC:
		# Temperature vs Pressure profiles
		fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(15,9))

		# RAW DATA
		p1 = ax[0].plot(data['TEMP'],data['PRES'],'-',lw=0.5,alpha=0.5)
		ax[0].set_xlabel('Temperature (degC)'); ax[0].set_ylabel('Pressure (dbar)')
		ax[0].set_title('Raw Data'); ax[0].invert_yaxis(); ax[0].grid()

		# QC DATA
		p2 = ax[1].plot(data['TEMP'].where(bad_pnts==0).sel(cycle=bad_prof==0),
					    data['PRES'].where(bad_pnts==0).sel(cycle=bad_prof==0),
					    '-',lw=0.5,alpha=0.5)
		ax[1].set_xlabel('Temperature (degC)'); ax[1].set_ylabel('Pressure (dbar)')
		ax[1].set_title('QC Data'); ax[1].invert_yaxis(); ax[1].grid()

		plt.savefig(join(figures_dir, 'QC_TEMP_figs', 'QCprofs_TEMP.{:02d}.pdf'.format(filenum)))

		# Salinity vs Pressure profiles
		fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(15,9))

		# RAW DATA
		p1 = ax[0].plot(data['PSAL'],data['PRES'],'-',lw=0.5,alpha=0.5)
		ax[0].set_xlabel('Salinity (PSU)'); ax[0].set_ylabel('Pressure (dbar)')
		ax[0].set_title('Raw Data'); ax[0].invert_yaxis(); ax[0].grid()

		# QC Data
		p2 = ax[1].plot(data['PSAL'].where(bad_pnts==0).sel(cycle=bad_prof==0),
						data['PRES'].where(bad_pnts==0).sel(cycle=bad_prof==0),
						'-',lw=0.5,alpha=0.5)
		ax[1].set_xlabel('Salinity (PSU)'); ax[1].set_ylabel('Pressure (dbar)')
		ax[1].set_title('QC Data'); ax[1].invert_yaxis(); ax[1].grid()

		plt.savefig(join(figures_dir, 'QC_PSAL_figs', 'QCprofs_PSAL.{:02d}.pdf'.format(filenum)))


	### Save flag arrays to NetCDF file ###
	if save_QC:
		# Convert Numpy flag arrays to xarray DataArrays
		bad_pnts_array = xr.DataArray(bad_pnts,dims=['sample','cycle'],
			attrs={'description':'Flagged points',
			 	   '0':'Good Data',
				   '1':'Missing or out of range pressure',
				   '2':'Missing or out of range temperature',
				   '3':'Missing or out of range salinity',
				   '4':'Bad temperature points hand picked by Annabel Wade',
				   '5':'Bad salinity points hand picked by Annabel Wade'})

		bad_prof_array = xr.DataArray(bad_prof,dims=['cycle'],
			attrs={'description':'Flagged profiles',
				   '0':'Good Data',
				   '1':'Repeated or out of order pressure', 
				   '2':'Salinity gradient (larger than 0.002 below 1000 dbar)', 
				   '3':'Both flag 1 and 2', 
		   		   '4':'4 or fewer measurements', 
				   '5':'Out of bounds time range (everything before 2022) or no time data' ,
				   '6':'No geo-coordinates'})

		special_areas_array = xr.DataArray(special_loc, dims=['cycle'], 
			attrs={'description':'Special regions flag',
					'0':'Non-special' ,
					'1':'Mediterranean Sea',
					'2':'Baltic Sea',
					'3':'Indian Ocean flow-through; coastal West Pacific Ocean',
					'4':'Gulf of Mexico; Caribbean Sea',
					'5':'Coastal Pakistan'})
		
		# Compile DataArrays
		ds = xr.Dataset(data_vars={'bad_pnts':bad_pnts_array,
								   'bad_prof':bad_prof_array,
								   'FLOATNUM':data['FLOATNUM'],
								   'CYCLENUM':data['CYCLENUM'],
								   'LAT':data['LAT'],
								   'LON':data['LON'], 
								   'special_areas':special_areas_array},
			attrs={'creators':'Katy Christensen & Annabel Wade',
				   'date_created':datetime.now().strftime("%Y-%m-%d"),
				   'description':'Flag arrays intended to be used with global Argo ' + 
				   				 'float dataset obtained by Robert Drucker on 26-Aug-2022'}) 

		# Save current .nc flag array file
		ds.to_netcdf(join(flag_array_dir, 'flag_arrays.{:02d}.nc'.format(filenum)))
