# Global Dynamic Height from Quality-Controlled Argo float profiles
Last updated April 1, 2024 by Annabel Wade and Katy Christensen.

## Description
These scripts carry out quality control and analysis of the global Argo data sets obtained by Robert Drucker on August 29th, 2022 in the form of 80 netcdf files (`global_profiles.nn.nc`, where `nn` denotes the file number), and calculate dynamic height with these data.  


## Files and order of usage:
1. `QC_flag_argo.py`
	- dependent on:
		- `QC_hand_pick_pnts.py`
		- `QC_functions.py`
		- `global_profiles.nn.nc` data files in *data_dir*
		- `special_areas_polygons.shp` file in *main_dir*
2. `dynamic_height_calculation.py`
	- dependent on: 
		- `QC_functions.py`
		- `global_profiles.nn.nc` data files in *data_dir*
		- `flag_array.nn.nc` files in *flag_array_dir*
3. `QC_analysis.ipynb`
	- dependent on:
		- `global_profiles.nn.nc` data files in *data_dir*
		- `flag_array.nn.nc` files in *flag_array_dir*
		- `hp_pnts` dictionary generated in `QC_hand_pick_pnts.py`
		- `global_dh.nc` files in *dynamic_height_dir*
		- `special_areas_polygons.shp` file in *main_dir*

## Directories
Note: if a directory does not exist, it is created. 

+ *main_dir*  
	- contains Python scripts and notebooks, *data_dir*, *figures_dir*, *flag_array_dir*, and *dynamic_height_dir*  
+ *data_dir*
	- contains data files  
+ *flag_array_dir*
	- contains flag array files  
+ *figures_dir* 
	- contains *QC_TEMP_figs*, *QC_PSAL_figs*, *hp_pnts_figs* 
+ *dynamic_height_dir*  
	- contains calculations of dynamic height for each data file at set levels  


## Python scripts/notebooks
+ `QC_flag_argo.py` - Main QC program that generates flag array files:
   - Inputs:
     - Directories - Uses current directory as *main_dir* by default, but path of *main_dir* can be inputted. Other directories are defined based on *main_dir*.
     - Data files - Contained in *data_dir*.
     - `hp_pnts` - Dictionary where keys are data file numbers and values are arrays of tuples in the form `(sample_index, cycle_index)` that represent handpicked PSAL points. Generated in `QC_hand_pick_pnts.py`.
	 - `special_areas_polygons.shp` polygons that represent special areas.
   - Outputs:
     - Prints name of data file when beginning QC for that data file, along with tqdm progress bar in console.
     - QC plots - Shows raw vs. QC data for each file. Saves TEMP plots to `QC_TEMP_figs` in *figures_dir* and PSAL plots to `QC_PSAL_figs` in *figures_dir*, formatted as `'QCprofs_TEMP.nn.pdf'` and `'QCprofs_PSAL.nn.pdf'`.
     - Flag array files - Saved to *flag_array_dir* as `'flag_arrays.nn.nc'`.
   - Program Parameters:
     - `flag_QC` - If True, will generate flag arrays for bad points and bad profiles based on method defined in Summary.
     - `plot_QC` - If True, will generate and save plots showing raw vs. QC data.
     - `save_QC` - If True, will save flag array files to *flag_array_dir*.
     - `overwrite_flag_files` - If True, will overwrite existing flag array files in *flag_array_dir*.
   - Notes:
     - If no data files are present in *data_dir*, `FileNotFoundError` will be raised.
     - Will only generate flag array files for data files that don't yet have a corresponding flag array file.

+ `QC_hand_pick_pnts.py` - Generates handpicked points dictionaries:
   - Inputs:
     - Directories - Uses current directory as `main_dir` by default, but path of `main_dir` can be inputted. Other directories are defined based on `main_dir`.
     - Data files - Contained in `data_dir`.
   -*Outputs:
     - `hp_pnts` - Python Dictionary where keys are data file numbers and values are arrays of tuples in the form `(sample_index, cycle_index)` that represent handpicked PSAL points.
   - Functions:
     - `pnt_tuples(filenum, cyc_i, sample_i_arr=None, all_pnts=False, start_sample=0, end_sample=-1)`:
       - Generates an array of tuples in the format `(sample_index, cycle_index)` where the indices correspond to a bad point to be flagged in the main program.
       - Arguments:
         - `filenum` - Number of the data file.
         - `cyc_i` - Index of the cycle to flag.
         - `sample_i_arr` - An array of sample indices can be provided.
         - `start_sample` & `end_sample` - Starting and ending (exclusive) sample indices can be provided.
         - `all_pnts` - True if all points of the file need to be flagged.
       - Returns: NumPy array of tuples in the format `(sample_i, cyc_i)`.


+ `dynamic_height_calculation.py` - Computes dynamic height for each data file at set levels:
   - Inputs:
     - Directories - Uses current directory as `main_dir` by default, but path of `main_dir` can be inputted. Other directories are defined based on `main_dir`.
     - Data files - Contained in `data_dir`.
     - Flag array files - Contained in `flag_array_dir`.
   - Outputs:
     - Dynamic heights - Dynamic height for 29 levels for each `'global_profiles.nn.nc'` file in `data_dir` with a corresponding  `'flag_arrays.nn.nc'` file. Outputs to `dynamic_height_dir`
   - Functions:
     - Uses functions from `QC_functions`.
     - `compute_dh_level(interp_func,upper_pres,lower_pres)`
        - Initialized with Dask. Creates delayed dask item to integrate the specific volume anomaly interpolated function between two levels.
        - Arguments:
          - `interp_func` - the pchip interpolated specific volume anomaly functions, scipy object
          - `upper_pres` - upper bound of integration
          - `lower_pres` - lower bound of integration
        - Returns: dynamic height for the `upper_pres` value as a single number

+ `QC_analysis.ipynb` - Analyzes QC flags and handpicked points, and generates plots and tables:
   - Inputs:
     - Directories - Uses current directory as `main_dir` by default, but path of `main_dir` can be inputted. Other directories are defined based on `main_dir`.
     - Data files - Contained in `data_dir`.
     - `hp_pnts` - Dictionary where keys are data file numbers and values are arrays of tuples in the form `(sample_index, cycle_index)` that represent handpicked PSAL points. Generated in `QC_hand_pick_pnts.py`.
     - `global_dh.nc` files in `dynamic_height_dir`.
   - Outputs:
     - Handpicked points plots - Multipanel visualization of handpicked salinity and temperature points. Saved in `hp_pnts_figs` folder in the `figures_dir`, as `'hp_TEMP.pdf'` and `'hp_PSAL.pdf'`.
     - Data information - Info on extent, size, etc. of data files. Printed to console.
     - Box and whisker plots - Distribution of QC flags. Saved in `figures_dir` as `'prof_flag_bxp.pdf'`, `'special_areas_flag_bxp.pdf'`, and `'pnts_flag_bxp.pdf'`.
     - Tables - Counts and percentages of QC flags. Printed to console in LaTeX table format.

+ `QC_functions.py` - Contains utility functions used in QC_flag_argo.py and dynamic_height_calculation.py::
	- Functions:
		- `find_nc_filenames(path_to_dir, suffix=".nc") `- Finds all files in given directory that end with the suffix .nc. 
			- Arguments: `path_to_dir` -- filepath to directory
			- Returns: list of filepaths for each .nc file in the directory
		- `make_dir(dirname, clean=False)`- Make a directory if it does not exist.
	    	- Use `clean=True` to clobber the existing directory.
		- `load_polygons(fp)` - Loads polygons from a file.
			- Arguments: `fp` -- filepath to .shp file containing polygon geometries
			- Returns: list of `shapely.geometry.Polygon` objects


## Package dependencies
#### name of package, version number
numpy, 1.21.5  
xarray, 2022.3.0  
pandas, 1.4.1  
python-dateutil, 2.8.2  
tqdm, 4.64.0  
matplotlib, 3.5.1  
gsw,  3.4.0
scipy, 1.7.3
dask, 2021.7.2

Python version used: 3.9.7

## Global Argo float data files
80 .nc files named with the format: `global_profiles.nn.nc`, where `nn` ranges from 01 to 80.

Variables (dimensions in parentheses):
- `FLOATNUM` (cycle)
- `CYCLENUM` (cycle)
- `DATENUM` (cycle)
- `LAT` (cycle)
- `LON` (cycle)
- `POS_QC` (cycle)
- `POS_ADJUSTED_QC` (cycle)
- `POS_RECOVERED` (cycle)
- `DATA_MODE` (cycle)
- `MAX_PRES_ADJUSTED_ERROR` (cycle)
- `PI_NAME` (string_length1, cycle)
- `DATA_CENTRE` (string_length2, cycle)
- `PRES` (sample, cycle)
- `PRES_OUTOFRANGE` (sample, cycle)
- `PRES_QC` (sample, cycle)
- `PRES_ADJUSTED_QC` (sample, cycle)
- `TEMP` (sample, cycle)
- `TEMP_OUTOFRANGE` (sample, cycle)
- `TEMP_QC` (sample, cycle)
- `TEMP_ADJUSTED_QC` (sample, cycle)
- `PSAL` (sample, cycle)
- `PSAL_OUTOFRANGE` (sample, cycle)
- `PSAL_QC` (sample, cycle)
- `PSAL_ADJUSTED_QC` (sample, cycle)
- `ProfilePressure` (cycle)
- `ProfilePressureSource` (cycle)
- `ParkingPressure` (cycle)
- `ParkingPressureSource` (cycle)

Further description of these variables is in Robert Drucker's 'Global Profiles.pdf'.

## Flag Array Files

80 .nc files named with the format: `flag_arrays.nn.nc`, where `nn` ranges from 01 to 80, each corresponding to a global profiles data file with the same file number `nn`.

Variables (dimensions in parentheses):
- `bad_pnts` (sample, cycle)
- `bad_prof` (cycle)
- `special_areas` (cycle)

Further description of these variables is in the QC Summary.


<!--		- uses code from QC_correction_argo.py** to prepared data for interpolation --!>
	
<!--**Auxiliary: QC_correction_argo.py - Isolated procedure for correcting profiles flagged 1-3
	inputs:
		directories - Uses current directory as main_dir by default
		data files - Contained in data_dir
	outputs:
		Prints name of data file when beginning QC for that data file along with tqdm progress bar.
		ds - variable in program that contains corrected profiles

	--!>