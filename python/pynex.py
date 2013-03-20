#!/usr/bin/env python

################################################################################
#
# pynex.py
#
# Usage: pynex.py [--trans_file=file] SEED list atom
#
# Constructs names "SEED_X" where X is a number from the file "list", then 
# combines all the NEXAFS spectra for atom "atom" together and writes to a file.
#
# Uses libpytep to do the heavy lifting. Can specify a transitions file that
# adds the transition energies to each spectrum.
#
################################################################################
#
# Created October2012 by Kane O'Donnell (Australian Synchrotron).
#
################################################################################
#
# Copyright 2012 Kane O'Donnell
#
#     This script is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This script is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this script.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# 
# NOTES
#
# 1. 
#
#
################################################################################

import argparse
import libpytep as lp
from numpy import zeros, savetxt, amin, amax, arange, array, linspace
from scipy.interpolate import interp1d

lpe = lp.core_level_spectra

parser = argparse.ArgumentParser(description="Calculate individual NEXAFS spectra from different core-hole sites.")

parser.add_argument('seed', help="CASTEP base seedname")
parser.add_argument('atlist', help="List of atomic suffixes.")
parser.add_argument('atom', help="Castep orbital number. This is hard to guess so just run pynex once and then use the number equal\nto the number of core projectors as usually the special core hole atom is last in the species.")
parser.add_argument('-t', '--transition_file', default=None, help="Transition energies file.")

args = parser.parse_args()

prefix = args.seed

# Open the suffixes file.
af = open(args.atlist, 'r')
atlist = [int(x) for x in " ".join(af.readlines()).split()]
af.close()

# If present get the transition energies
if args.transition_file is not None:
  tf = open(args.transition_file, 'r')
  transitions = [float(x) for x in " ".join(tf.readlines()).split()]
  tf.close()
else:
  transitions = None

species = int(args.atom)
spectra = []

for i, a in enumerate(atlist):
  print "Processing seed: ", prefix+ "_" + str(a)
  seed = prefix + "_" + str(a)
  lpe.init(seed)
  # Set spectral preferences here!
  lpe.lorentzian_width = 0.05
  lpe.gaussian_broadening = 0.05
  lpe.linear_broadening = 0.1
  lpe.generate_single_spectrum(species)
  if transitions is not None:
    # If you use the convention that an energy difference is excited - neutral
    # in Gao's transition energy calculation, the transition energy is POSITIVE,
    # so add it here. If the transition energy is negative, subtract here instead.
    spectra.append([lpe.w.copy() + transitions[i], lpe.lastspec.copy()])
  else:
    spectra.append([lpe.w.copy(), lpe.lastspec.copy()])

# Output the single atomic spectra
for i, xy in enumerate(spectra):
  x = xy[0]
  y = xy[1]
  tmp = zeros((y.shape[0], y.shape[1]+1))
  tmp[:,0] = x
  tmp[:,1:] = y
  savetxt(prefix+ "_" + str(atlist[i]) + ".nexafs", tmp)

# Combine the individual spectra. This is non-trivial in general because the
# different transition energies make the x-axis different. So we divide into 
# two cases:
#
# 1. We aren't using transition energies. In this case, just add the spectra 
#    together pointwise because they are all on the same x-scale. Note that
#    this is totally invalid in a core-hole calculation with different core-hole
#    sites.
#
# 2. We're using transition energies. Now what we do is figure out the step size
#    of the spectra (will all be the same) and the global minimum and maximum
#    energy, then generate a linspace over that range with the same step size. 
#    In this way the point density of the original spectra stays intact, but it
#    means the number of points in the global spectrum is not the same as the
#    individual spectra.

if transitions is None:
  print "Generating combined output with no transition energies."
  tmp = zeros((spectra[0][1].shape[0], spectra[0][1].shape[1] + 1))
  for xy in spectra:
    x = xy[0]
    y = xy[1]
    tmp[:,1:] += y
  tmp[:,0] = spectra[0][0]
  savetxt(prefix + "_combined.nexafs", tmp)
else:
  # Get the step size
  step = spectra[0][0][1] - spectra[0][0][0]
  # Find global min and max energy.
  minx = None
  maxx = None
  for xy in spectra:
    if minx is None:
      minx = amin(xy[0])
    elif amin(xy[0]) < minx:
      minx = amin(xy[0])
    if maxx is None:
      maxx = amax(xy[0])
    elif amax(xy[0]) > maxx:
      maxx = amax(xy[0])
  # Generate a spectrum between minx and maxx with the given step. We want to
  # include the endpoint here so we have to add a step.
  x = arange(minx, maxx+float(step)/2, float(step)/2)
  #x = linspace(minx, maxx, 2000)

  print "Generating combined output on xaxis starting at ", minx, ", stopping at ", maxx, " and with step size ", step/2, "."
  print "This gives ", len(x), " points in the spectrum."
  ## Now generate interpolators and interpolate 
  tmp = zeros((len(x), 7))
  for xy in spectra:
    xt = xy[0]
    yt = xy[1]
    for i in range(6):
      ixy = interp1d(xt, yt[:,i], bounds_error=False, fill_value=0.0)
      tmp[:,i+1] += array([ixy(E) for E in x])
  tmp[:,0] = x
  savetxt(prefix + "_combined.nexafs", tmp)
    
