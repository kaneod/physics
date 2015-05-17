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
#    outputs from the Quantum Espresso code XSpectra in various ways. To get 
#    XSpectra output to work, you need a series of files named xanes.atomX.ABC.dat
#    where X is an atom number (1,2,3 etc) and ABC is each combination of {100, 110,
#    010, 011, 001, 101}, representing XSpectra outputs with E-vectors of each of
#    those orientations.
#
# 2. The output file format is the same for all files and consists of the headings E,
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
import imp

Ry2eV = 13.605698066
SMALL = 1.0e-6    			# small floating point number for equality comparisons
DEBUG = 1								# Print debug statements
FIX_OCC_CUT_DEFECT = 1	# XSpectra's method of cutting occupied states leaves a defect,
												# where the first point below zero eV is too large. Fix by
												# interpolating from the surrounding points.
NUM_POINTS = 2000   		# Number of interpolation points for combining spectra.

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
  spec = zeros((components['E'].shape[0], 2))
  for i in range(spec.shape[0]):
    spec[i,0] = components['E'][i]
    spec[i,1] = (nvec[0] ** 2) * components['xx'][i] + (nvec[1] ** 2) * components['yy'][i] + \
                (nvec[2] ** 2) * components['zz'][i] + \
                (nvec[0] * nvec[1]) * components['xy'][i] + \
                (nvec[0] * nvec[2]) * components['xz'][i] + \
                (nvec[1] * nvec[2]) * components['yz'][i]
  
  return spec
  
def spherical_to_cartesian(theta, phi):
  """ Converts spherical coordinates (theta, phi) to cartesian (x,y,z). Theta and phi in
  degrees, not radians. Note theta is azimuth from x in x-y plane, phi is angle from z-axis."""
  
  tr = pi * theta * 1.0 / 180 # Careful to avoid int division
  pr = pi * phi * 1.0 / 180
  
  x = cos(tr) * sin(pr)
  y = sin(tr) * sin(pr)
  z = cos(pr)
  
  return array([x,y,z])
  
def cartesian_to_spherical(vec):
  """ Converts a cartesian vector to spherical coords IGNORING THE RADIUS. That is, only theta 
  and phi are returned. Careful! """
  
  r = norm(vec)
  pr = arccos(vec[2] / r)
  tr = arctan2(vec[1], vec[0])
  
  theta = 180.0 * tr / pi
  phi = 180.0 * pr / pi
  
  return theta, phi
  
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
    
  # Check that it really is an in-plane vector. If not, project to plane and warn.
  if abs(dot(evec, array([0.0, 0.0, 1.0]))) > SMALL:
    print "Warning - in-plane vector isn't actually in plane!"
    vec = array([evec[0], evec[1], 0])
  else:
    vec = evec
    
  # If the angle is 90 degrees, there is only one answer (technically violates
  # the hypotheses because x, y components not the same but they have to be zero here anyway!
  if abs(angle - 90.0) < SMALL:
    return array([0, 0, 1.0])
  else:
    z = norm(vec) * tan(angle * pi / 180.0) 
    return array([vec[0], vec[1], z])
    
def prune_duplicate_vectors(vecs):
  """ Returns a list of unique vectors, where two vectors are unique if they are not
  parallel *or* anti-parallel."""
  
  pruned = [vecs[0]]
  
  for v in vecs[1:]:
    dupe = False
    for p in pruned:
      if abs(dot(v, p) - norm(p) * norm(v)) < SMALL:
        dupe = True
        break
      elif abs(dot(v, p) + norm(p) * norm(v)) < SMALL:
        dupe = True
        break
    if not dupe:
      pruned.append(v)
  
  if DEBUG:
    print "Inside prune_duplicate_vectors, input list was:"
    print vecs
    print "Output list is:"
    print pruned
    
  return pruned
        
def fix_occupancy_cutting_defect(esigma):
	""" Given the raw XSpectra E, sigma input, finds the sigma value closest to but below
	zero eV and replaces it with a linear interpolation from the nearest points. This is to
	correct the intensity "spike" that XSpectra puts there as a result of a somewhat faulty
	occupancy cutting algorithm. """
	
	# First find index of smallest non-negative energy.
	smallest = -100000
	index = 0
	while esigma[index,0] < 0:
		index += 1
	
	x1 = esigma[index-2,0]
	x2 = esigma[index,0]
	xi = esigma[index-1,0]
	for j in range(1,esigma.shape[1]):
		y1 = esigma[index-2,j]
		y2 = esigma[index,j]
		esigma[index-1,j] = y1 + ((y2 - y1)/(x2 - x1)) * (xi - x1)
	
	return esigma
	
def combine_two_spectra(s1, s2, npts):
  """ Given two spectra on possibly slightly-different X-axes, combine into one with a
  common X-axis defined by the overlapping region."""
  
  # Generate the new axis
  newmin = max(amin(s1[:,0]), amin(s2[:,0]))
  newmax = min(amax(s1[:,0]), amax(s2[:,0]))
  newx = linspace(newmin, newmax, npts)
  
  # Interpolate the existing y ranges onto the new grid.
  ny1 = interp(newx, s1[:,0], s1[:,1])
  ny2 = interp(newx, s2[:,0], s2[:,1])
  
  # Add the two together and put in an energy column.
  spec = zeros((npts, 2))
  spec[:,0] = newx
  spec[:,1] = ny1[:] + ny2[:]
  
  return spec
  
def combine_spectra(spectra, npts):
  """ Use combine_two_spectra to combine a list of spectra arrays. """
  
  if len(spectra) == 0:
    return []
  else:
    combined = spectra.pop()
    while len(spectra) > 0:
      combined = combine_two_spectra(combined, spectra.pop(), npts)
    return combined 
  
  
def combine_two_sites(s1, s2, npts):
  """ Rather than re-use combine_spectra, we rewrite specifically for a components
  structure. That is, s1 and s2 are dictionaries with keys 'E', 'xx' etc. """
  
  combined = {}
  # Generate the new axis
  newmin = max(amin(s1['E']), amin(s2['E']))
  newmax = min(amax(s1['E']), amax(s2['E']))
  combined['E'] = linspace(newmin, newmax, npts)
  
  
  for key in ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']:
    # Interpolate components onto new grid.
    newy1 = interp(combined['E'], s1['E'], s1[key])
    newy2 = interp(combined['E'], s2['E'], s2[key])
    combined[key] = newy1[:] + newy2[:]
  
  return combined
  
def combine_sites(sites, npts):
  """ Uses combine_two_sites to interpolate all the atomic site components together. """
  
  if len(sites) == 0:
    if DEBUG:
      print "Warning (combine_sites): tried to combine an empty list of sites."
    return []
  else:
    combined = sites.pop(sites.keys()[-1])
    while len(sites) > 0:
      combined = combine_two_sites(combined, sites.pop(sites.keys()[-1]), npts)
    return combined

parser = argparse.ArgumentParser(description="Take XSpectra output from multiple atoms and combine it.")

parser.add_argument('atoms', type=int, nargs="+", help="Indices of atomic absorber sites.")
parser.add_argument('--sites', '-s', dest='sites', action='store_true', default=False, help="Output combined spectra for each individual atomic site.")
parser.add_argument('--angles', '-a', dest='angles', action='store_true', default=False, help="Generate 20, 35, 55, 75, 90 degree spectra instead of components.")
parser.add_argument('--symmetry', '-y', dest='symmetry', action='store_true', default=False, help="Look for symmetries function and apply symmetries.")
args = parser.parse_args()

if DEBUG:
  print "DEBUG input check:"
  print "atoms = ", args.atoms
  print "sites = ", args.sites
  print "angles = ", args.angles
  print "symmetry = ", args.symmetry
  
# Init task: check atoms input is valid/remove duplicates.
indices = set(args.atoms)

# Init task: decide whether we have and/or need an external symmetry function, and
# assign as necessary.
if args.symmetry:
  try:
    imp.find_module('qnex_symmetries')
    sym_found = True
  except ImportError:
    sym_found = False
  if sym_found:
    import qnex_symmetries
    symmetry_function = qnex_symmetries.symmetry_generate
  else:
    symmetry_function = symmetry_stub
else:
  symmetry_function = symmetry_stub


# First real task is to load all the spectra. Structure is a dict of atoms where each atom is a dictionary keyed by '100' etc.
atoms = {}

for i in indices:
  this_atom = {}
  for v in ['100', '110', '010', '011', '001', '101']:
    tmp = loadtxt("xanes.atom" + str(i) + "." + v + ".dat")
    if FIX_OCC_CUT_DEFECT:
    	tmp = fix_occupancy_cutting_defect(tmp)
    # If the calculation was spin polarized this file will have four columns and we only
    # want the "up" spin (that's where the excited electron goes).
    if tmp.shape[1] == 4:
      this_atom[v] = tmp[:,0::2].copy() ## AHA! I have defeated you, you sly View beast!
    else:
      this_atom[v] = tmp.copy()
  
  atoms[i] = this_atom

# Next task is to use the 100, etc outputs to generate components xx, yy, etc.
# A single 7-column array might be "neater" here but the dictionary makes the code easier
# to read/think about.

atom_components = {}

for i in indices:
  this_atom = {}
  this_atom['E'] = atoms[i]['100'][:,0]
  this_atom['xx'] = atoms[i]['100'][:,1]
  this_atom['yy'] = atoms[i]['010'][:,1]
  this_atom['zz'] = atoms[i]['001'][:,1]
  this_atom['xy'] = 2.0 * atoms[i]['110'][:,1] - atoms[i]['100'][:,1] - atoms[i]['010'][:,1]
  this_atom['xz'] = 2.0 * atoms[i]['101'][:,1] - atoms[i]['100'][:,1] - atoms[i]['001'][:,1]
  this_atom['yz'] = 2.0 * atoms[i]['011'][:,1] - atoms[i]['010'][:,1] - atoms[i]['001'][:,1]
  atom_components[i] = this_atom

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
  transitions = {}
  if "ground" in total_energies.keys():
    ground = total_energies['ground']
    for i in indices:
      transitions[i] = total_energies[str(i)] - ground
  else:
    # Dodgy hack for average of a dictionary.
    ground = 0
    for k in total_energies.keys():
      ground += total_energies[k]
    ground /= len(total_energies.keys())
    for i in indices:
      transitions[i] = total_energies[str(i)] - ground
  # Now we need the difference between the HOMO and LUMO and subtract this off the
  # existing transitions.
  if os.path.isfile("frontier_levels"):
    tmp = {}
    f = open("frontier_levels", 'r')
    lines = f.readlines()
    f.close()
    for l in lines:
    	tmp[l.split()[0]] = float(l.split()[2]) - float(l.split()[1]) # Elumo - Ehomo
    for i in indices:
      transitions[i] -= tmp[str(i)]
  # Now set up the energy scale for all the atomic sites.
  for i in indices:
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
    for i in indices:
      f = open("qnex.atom" + str(i) + ".dat", 'w')
      f.write("# NEXAFS components output by qnex_combine.py\n")
      f.write("# Atom: "+ str(i) + "\n")
      f.write("# E    xx    yy    zz    xy    xz    yz\n")
      tmp = atom_components[i]
      for j in range(atom_components[i]['E'].shape[0]):
        f.write("%f    %f    %f    %f    %f    %f    %f\n" % (tmp['E'][j], tmp['xx'][j], tmp['yy'][j], tmp['zz'][j], tmp['xy'][j], tmp['xz'][j], tmp['yz'][j]))
    f.close()
  else:
    combined_components = combine_sites(atom_components, NUM_POINTS)
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

# Regardless of our output choice (combined or separate) we need angle spectra generated for all
# angles for all atomic sites, so we do that first.

if args.angles:
  base_vecs = [spherical_to_cartesian(45,20), spherical_to_cartesian(45, 35), spherical_to_cartesian(45,55), spherical_to_cartesian(45,75), spherical_to_cartesian(45,90)]
  base_names = ["20", "35", "55", "75", "90"]
  
  # Treat each angle separately even if combining atomic sites.
  for v,b in zip(base_vecs, base_names):
    sym_vecs = symmetry_function(v)
    # Some of these might be duplicates, so we have to prune the list.
    sym_vecs = prune_duplicate_vectors(sym_vecs)
    # For each symmetry vector, need a spectrum for each site. Since we always add the 
    # symmetry-related spectra, there is no need to store separately - just add them 
    # together as they are calculated.
    atomic_spectra = []
    for i in indices:
      tmp_specs = []
      for v in sym_vecs:
        tmp_specs.append(calculate_spectrum(v, atom_components[i]))
      # There might be a different number of symmetry vectors for difference angles,
      # so normalize by the number of vectors in the set.
      combined = combine_spectra(tmp_specs, NUM_POINTS)
      combined[:,1] /= len(sym_vecs)
      atomic_spectra.append(combined)
    
    # If we want site-specific output, do this now.
    if args.sites:
      for i in indices:
        f = open("qnex.atom" + str(i) + ".angle" + b + ".dat", 'w')
        f.write("# NEXAFS angle output by qnex_combine.py\n")
        f.write("# Atom: " + str(i) + " Angle: " + b + " degrees.\n")
        f.write("# E    sigma\n")
        for j in range(atomic_spectra[i].shape[0]):
          f.write("%f    %f\n" % (atomic_spectra[i][j,0], atomic_spectra[i][j,1]))
        f.close()
    else:
      # Combine spectra, and output.
      combined = combine_spectra(atomic_spectra, NUM_POINTS)
      f = open("qnex.combined.angle" + b + ".dat", 'w')
      f.write("# NEXAFS angle output by qnex_combine.py\n")
      f.write("# All atomic sites combined, angle " + b + " degrees.\n")
      f.write("# E    sigma\n")
      for j in range(combined.shape[0]):
        f.write("%f    %f\n" % (combined[j,0], combined[j,1]))
      f.close()
      
print "Finished qnex_combine!"