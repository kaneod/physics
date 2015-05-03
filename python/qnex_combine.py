#!/usr/bin/env python
################################################################################
#
# qnex_combine.py
#
# Given a radius in angstroms and an atomic index, adds a suffix to all atoms
# within range of that specified atom. Handles periodic boundary conditions.
#
################################################################################
#
# Copyright 2015 Kane O'Donnell
#
#     This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# 
# NOTES
#
# 1. This is a fairly specialised helper code that combines specifically-named
# 	 outputs from the Quantum Espresso code XSpectra in various ways. To get 
#    XSpectra output to work, you need a series of files named xanes.atomX.ABC.dat
#    where X is an atom number (1,2,3 etc) and ABC is each combination of {100, 110,
#    010, 011, 001, 101}, representing XSpectra outputs with E-vectors of each of
#		 those orientations.
#
#	2. The output file format is the same for all files and consists of the headings E,
#    xx, yy, zz, xy, xz, yz where the latter six are the cross section components. If
#    the files total_energies and frontier_levels exist in the directory where the 
#    script is run, the energy scale will properly incorporate chemical shifts and the
#    DFT transition energy. The format of these files is revealed in the code below.
#
# 3. There is a symmetrized output option which looks for a symmetry-generating function
#    in the directory where the script is executed, if not there is a stub function
#    in this file. All the symmetry function has to do is to take a 3-vector (numpy array)
#    and return a list of symmetry-related 3-vectors (numpy arrays). The module
#    must be called qnex_symmetries and the function must be called symmetry_generate.
#
################################################################################

from __future__ import division
from numpy import *
from numpy.linalg import norm
import argparse
import sys
import os.path

Ry2eV = 13.605698066
SMALL = 1.0e-6		# small floating point number for equality comparisons
DEBUG=1

######################
#
# Functions
# 
######################

def calculate_spectrum(evec, components):
	""" Uses the passed electric field vector to combine components into a particular angle
	of NEXAFS spectrum."""
	
	# Check evec is normalized.
	nvec = evec / norm(evec)
	
	# Now generate spectrum. 
	spec = zeroes((components['E'].shape[0], 2))
	for i in range(spec.shape[0]):
		spec[i,0] = components['E'][i]
		spec[i,1] = (nvec[0] ** 2) * components['xx'][i] + (nvec[1] ** 2) * components['yy'][i] + \
								(nvec[2] ** 2) * components['zz'][i] + \
								(nvec[0] * nvec[1]) * components['xy'][i] + \
								(nvec[0] * nvec[2]) * components['xz'][i] + \
								(nvec[1] * nvec[2]) * components['yz'][i]
	
	return spec
	
def symmetry_stub(evec):
	""" This is the symmetry function used if no other is provided - all it does is apply
	the symmetry y -> x, x -> -y. """

	evec2 = evec.copy()
	evec2[0] = evec[1]
	evec2[1] = -1 * evec[0]
	evec2[2] = evec[2]
	return [evec, evec2]
	
def angle_vec(evec, angle):
	""" Given an in-plane vector evec, returns a vector with the same in-plane components
	but with an out-of-plane component set so the final vector is at an angle "angle" w.r.t.
	the original vector. """
	
	# First, if the angle is 90 degrees, there is only one answer (technically violates
	# the hypotheses because x, y components not the same but they have to be zero here anyway!
	if (angle - 90.0) < SMALL:
		return array([0, 0, 1])
	else:
		factor = (cos(angle * pi / 180.0))**(-2.0) - 1.0
		return array([evec[0], evec[1], norm(evec) * factor])

def combine_spectra(x1, y1, x2, y2, npts):
	""" Given two spectra on possibly slightly-different X-axes, combine into one with a
	common X-axis defined by the overlapping region."""
	
	# Generate the new axis
	newmin = max(amin(x1), amin(x2))
	newmax = min(amax(x1), amax(x2))
	newx = linspace(newmin, newmax, npts)
	
	# Interpolate the existing y ranges onto the new grid.
	ny1 = interp(newx, x1, y1)
	ny2 = interp(newx, x2, y2)
	
	# Add the two together and put in an energy column.
	spec = zeros((npts, 2))
	spec[:,0] = newx
	spec[:,1] = ny1[:] + ny2[:]
	
	return spec
	
def combine_sites(sites):
	""" Uses combine_spectra to combine the cross section components from multiple atomic
	sites. """
	
	

parser = argparse.ArgumentParser(description="Take XSpectra output from multiple atoms and combine it.")

parser.add_argument('num_atoms', type=int, help="Number of atomic absorber sites.")
parser.add_argument('--sites', '-s', dest='sites', action='store_true', default=False, help="Output combined spectra for each individual atomic site.")
parser.add_argument('--angles', '-a', dest='angles', action='store_true', default=False, help="Generate 20, 35, 55, 75, 90 degree spectra instead of components.")
parser.add_argument('--symmetry', '-y', dest='symmetry', action='store_true', default=False, help="Look for symmetries function and apply symmetries.")
args = parser.parse_args()

# First task is to load all the spectra. Structure is a list of atoms where each atom is a dictionary keyed by '100' etc.
atoms = []

for i in range(args.num_atoms):
	this_atom = {}
	for v in ['100', '110', '010', '011', '001', '101']:
		tmp = loadtxt("xanes.atom" + str(i+1) + "." + v + ".dat")
		# If the calculation was spin polarized this file will have four columns and we only
		# want the "up" spin (that's where the excited electron goes).
		if tmp.shape[1] == 4:
			this_atom[v] = tmp[:,0::2].copy() ## AHA! I have defeated you, you sly View beast!
		else:
			this_atom[v] = tmp.copy()
	atoms.append(this_atom)

# Next task is to use the 100, etc outputs to generate components xx, yy, etc.
# A single 7-column array might be "neater" here but the dictionary makes the code easier
# to read/think about.

atom_components = []

for i in range(args.num_atoms):
	this_atom = {}
	this_atom['E'] = atoms[i]['100'][:,0]
	this_atom['xx'] = atoms[i]['100'][:,1]
	this_atom['yy'] = atoms[i]['010'][:,1]
	this_atom['zz'] = atoms[i]['001'][:,1]
	this_atom['xy'] = 2.0 * atoms[i]['110'][:,1] - atoms[i]['100'][:,1] - atoms[i]['010'][:,1]
	this_atom['xz'] = 2.0 * atoms[i]['101'][:,1] - atoms[i]['100'][:,1] - atoms[i]['001'][:,1]
	this_atom['yz'] = 2.0 * atoms[i]['011'][:,1] - atoms[i]['010'][:,1] - atoms[i]['001'][:,1]
	atom_components.append(this_atom)

# Regardless of whether we want output as components or angles, we have to have everything
# on a sensible energy scale, so we do this now. Note that we exit with an error message
# if no energy scale is available and combined spectra are requested.

if os.path.isfile("total_energies"):
	f = open("total_energies", 'r')
	lines = f.readlines()
	f.close()
	total_energies = {}
	for l in lines:
		total_energies[l.split()[0]] = float(l.split()[1]) * Ry2eV
	# We need the ground state to get a decent transition energies, otherwise just use the average
	# of the energies and the difference as the shift.
	transitions = []
	if "ground" in total_energies.keys():
		ground = total_energies['ground']
		for i in range(args.num_atoms):
			transitions.append(total_energies[str(i+1)] - ground)
	else:
		# Dodgy hack for average of a dictionary.
		ground = 0
		for k in total_energies.keys():
			ground += total_energies[k]
		ground /= len(total_energies.keys())
		for i in range(args.num_atoms):
			transitions.append(total_energies[str(i+1)] - ground)
	# Now we need the difference between the HOMO and LUMO and subtract this off the
	# existing transitions.
	if os.path.isfile("frontier_levels"):
		tmp = loadtxt("frontier_levels")
		for i in range(args.num_atoms):
			transitions[i] -= tmp[i,2] - tmp[i,1]
	# Now set up the energy scale for all the atomic sites.
	for i in range(args.num_atoms):
		atom_components[i]['E'] += transitions[i]

if DEBUG:
	# Print some output at this point to check we're not crazy.
	print transitions
	
# Now we enter the output phase. Four cases: we are either outputting angles/components, combined/separate.

# CASES 1 and 2: component output, separate or combined.
# Note in the output, the xy, xz and yz components are:
# XY = <f|del_x|i><f|del_y|i>* + <f|del_y|i><f|del_x|i>* etc
# This is consistent with XX being |<f|del_x|i>|^2 etc. So, we can only get the "square"
# of matrix elements and not the actual complex matrix elements themselves.
if not args.angles:
	if args.sites:
		for i in range(args.num_atoms):
			f = open("qnex.atom" + str(i+1) + ".dat", 'w')
			f.write("# NEXAFS components output by qnex_combine.py\n")
			f.write("# Atom: "+ str(i+1) + "\n")
			f.write("# E    xx    yy    zz    xy    xz    yz\n")
			tmp = atom_components[i]
			for j in range(atom_components[i]['E'].shape[0]):
				f.write("%f    %f    %f    %f    %f    %f    %f\n" % (tmp['E'][j], tmp['xx'][j], tmp['yy'][j], tmp['zz'][j], tmp['xy'][j], tmp['xz'][j], tmp['yz'][j]))
		f.close()
	else:
		combined_components = combine_sites(atom_components)
		f = open("qnex.combined_sites.dat", 'w')
		f.write("# NEXAFS components output by qnex_combine.py\n")
		f.write("# All atomic sites combined.\n")
		f.write("# E    xx    yy    zz    xy    xz    yz\n")
		for j in range(combined_components['E'].shape[0]):
			f.write("%f    %f    %f    %f    %f    %f    %f\n" % (combined_components['E'][j], combined_components['xx'][j], combined_components['yy'][j], combined_components['zz'][j], combined_components['xy'][j], combined_components['xz'][j], combined_components['yz'][j]))
		f.close()

# CASES 3 and 4: we want ANGLE output. Angles are taken with respect to the Z axis. We now also
# need to consider whether symmetries are to be used. By default, the only symmetry applied is that
# both X and Y contributions are equal, that is, the angle vectors are all x + y + fz where f is 
# a fraction used to get the angle with respect to the sample correct. 
# To generate your own symmetries, you simply need a function called symmetry_generate in a module
# called qnex_symmetries in the calling directory. The function has to take a single vector and 
# return a list of vectors that are equivalent. Just to be precise, these vectors are made up
# of the electric field components ex, ey, ez, and these generate the coefficients ex^2, ey^2,
# ez^2, ex*ey, ex*ez and ey*ez. This code will always normalize all vectors prior to generating
# a spectrum.

