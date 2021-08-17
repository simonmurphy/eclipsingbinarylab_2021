import numpy as np
from scipy import interpolate
from astropy.table import Table

def BHAC15(iso_ages,mass_range,tracks='BHAC15_tracks+structure.txt'):
	""" Stripped back to just return mass, radius and age """
	tracks = Table.read(tracks,format='ascii',comment='!')
	masses = np.unique(tracks['col1'])
	idx = (masses >= mass_range[0]) & (masses <= mass_range[1])
	masses = masses[idx]
	iso_ages = np.log10(iso_ages) + 6
	mass_interp = []
	age_interp = []
	radius_interp = []
	for m in masses:
		idx = tracks['col1'] == m
		tmp = tracks[idx]
		int_radius = interpolate.interp1d(tmp['col2'],tmp['col6'])
		radius = int_radius(iso_ages)
		mass_interp = np.concatenate((mass_interp,np.ones(len(iso_ages))*m))
		radius_interp = np.concatenate((radius_interp,radius))
		age_interp = np.concatenate((age_interp,10**(iso_ages - 6)))

	return mass_interp,radius_interp,age_interp