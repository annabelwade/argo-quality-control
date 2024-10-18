### Last updated April 1, 2024 by Annabel Wade

# External imports
import numpy as np, xarray as xr, matplotlib.pyplot as plt, os
plt.rcParams['figure.dpi'] = 300

# Directories
main_dir = os.getcwd()  									# main / current directory
data_dir = os.path.join(main_dir,"data_dir")				# data files 
os.chdir(main_dir)

def pnt_tuples(filenum,cyc_i,sample_i_arr=None,all_pnts=False,start_sample=0, end_sample=-1):
	"""
	Generates an array of tuples in the format (sample_index, cycle_index) where the indices
	correspond to a bad point to be flagged in main program. 
	
	Arguments:
	filenum -- number of the data file
	cyc_i -- index of the cycle to flag
	sample_i_arr -- An array of sample indices can be provided
	start_sample & end_sample -- Starting and ending (exclusive) sample indicies can be provided 
	all_pnts -- True if all points of the file need to be flagged

	Returns: NumPy array of tuples in the format (sample_i, cyc_i)
	"""

	result = []

	if all_pnts: # all points in profile are bad
		nc = os.path.join( data_dir, 'global_profiles.{:02d}.nc'.format(filenum) )
		data = xr.open_dataset(nc)
		sample_num = len(data.isel(cycle=cyc_i).sample)
		for i in range(sample_num):
			result.append((i,cyc_i))
	elif not sample_i_arr == None: # sample indicies to flag were specified 
		for sample_i in sample_i_arr:
			result.append((sample_i,cyc_i))
	else: # start and end index were specified
		for sample_i in range(start_sample, end_sample):
			result.append((sample_i,cyc_i))

	return np.array(result)

# Dictionary mapping filenum to array of tuples in the form (sample, cycle)
hp_PSAL = dict()
hp_TEMP = dict()

# File 2
hp_PSAL[2] = np.array([(44,22958)])

# File 11
hp_PSAL[11] = np.array([(96,7144)])

# File 15
file15 = np.vstack((pnt_tuples(15,4878,all_pnts=True),
					pnt_tuples(15,5806,all_pnts=True),
					np.array([(58,10980)]),
					np.array([(58,10991)])))
hp_PSAL[15] = file15

# File 16
hp_PSAL[16] = np.array([(54,534)])

# File 20
file20 = np.vstack((pnt_tuples(20,9984,start_sample=196,end_sample=204),
					pnt_tuples(20,13837,start_sample=419,end_sample=422),
					pnt_tuples(20,14320,start_sample=264,end_sample=270),
					pnt_tuples(20,15994,start_sample=386,end_sample=392)))
hp_PSAL[20] = file20

# File 37
hp_PSAL[37] = np.vstack((pnt_tuples(37,2827,start_sample=483,end_sample=485),
							pnt_tuples(37,7596,start_sample=49,end_sample=53)))

# File 62
hp_PSAL[62] = pnt_tuples(62,416,sample_i_arr=[483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 494])

# File 67
filenum=67
nc = os.path.join( data_dir, 'global_profiles.{:02d}.nc'.format(filenum) ) 
data = xr.open_dataset(nc)
sample_num = len(data.isel(cycle=24886).sample)
hp_PSAL[67] = pnt_tuples(67,24886,start_sample=410,end_sample=sample_num)

# File 12
hp_PSAL[12] = np.vstack((pnt_tuples(12,19073,all_pnts=True),
						pnt_tuples(12,1542,sample_i_arr=[16])))

# File 13
hp_PSAL[13] = pnt_tuples(13,2609,sample_i_arr=[35,39,41,58])

# File 14
hp_PSAL[14] = np.vstack((pnt_tuples(14,27807,sample_i_arr=[0]),
						pnt_tuples(14, 27930, all_pnts=True), # bad DH profile due to unreasonable anomalousness. added 1/6/24
					))

# File 17
hp_PSAL[17] = np.vstack((pnt_tuples(17,6388, all_pnts=True), # Southern Ocean, instrument failure . added 1/6/24
						pnt_tuples(17,22993, all_pnts=True), # near Arabian Sea, unreasonable anomalousness. added 1/6/24 
						pnt_tuples(17,23425, all_pnts=True), # near Arabian Sea, unreasonable anomalousness. added 1/6/24
						))

# File 40
hp_PSAL[40] = pnt_tuples(40,5653,all_pnts=True)

# file 46
hp_PSAL[46] = np.vstack((pnt_tuples(46, 21031,all_pnts=True),
						pnt_tuples(46, 21032,all_pnts=True),
						pnt_tuples(46, 21037,all_pnts=True)))

# File 48
hp_PSAL[48] = pnt_tuples(48,4015,sample_i_arr=[0])

# File 78
hp_PSAL[78] = pnt_tuples(78,12,sample_i_arr=[1]) # PSAL>40 single pnt

# File 16
hp_TEMP[16] = np.vstack((pnt_tuples(16,10969,sample_i_arr=[2]),
					pnt_tuples(16,10232,sample_i_arr=[5])))

# File 17
hp_TEMP[17] = np.vstack((pnt_tuples(17,21978,sample_i_arr=[2,4]),
						pnt_tuples(17,22119,sample_i_arr=[2,5]),
						pnt_tuples(17,22259,sample_i_arr=[2,3,4]),
						pnt_tuples(17,21260,sample_i_arr=[3]),
						pnt_tuples(17,22118,sample_i_arr=[5])))

# File 32
hp_TEMP[32] = np.array([(69,2770)])
