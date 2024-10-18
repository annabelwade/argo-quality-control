### Last updated April 1, 2024 by Annabel Wade

# External imports
import os, numpy as np, shutil
from os.path import join

### Functions to use for QC methods and plotting ###

def find_nc_filenames(path_to_dir, suffix=".nc"):
	"""
	Finds all files in given directory that end with the suffix '.nc' 
	unless file suffix otherwise specified.

	Input:
		path_to_dir: filepath to directory, string

	Returns:
		list of filepaths for each .nc file in the directory
	"""

	filenames = os.listdir(path_to_dir)
	return [ os.path.join(path_to_dir, filename) for filename in filenames if filename.endswith( suffix ) ]

def make_dir(dirname, clean=False):
    """
    Make a directory if it does not exist.

    Input:
		dirname: directory name string
		clean: boolean, whether to clobber/overwrite the directory if it already exists
    """
    if clean == True:
        shutil.rmtree(dirname, ignore_errors=True)
        os.mkdir(dirname)
    else:
        try:
            os.mkdir(dirname)
        except OSError:
            pass # assume OSError was raised because directory already exists

def load_polygons(fp):
	"""
	Input:
		fp: filepath to .shp file containing polygon geometries
			
	Returns:
		list of shapely.geometry.Polygon objects
	"""
	import pickle
	with open(fp, 'rb') as f:
		loaded_lst = pickle.load(f)
	return loaded_lst
