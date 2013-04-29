#!/usr/bin/env python

################################################################################
#
# pynex.py
#
# Usage: pynex.py preferences.in
# 
# Contents of prefs file (see more detail below).
# seed - the base seedname of the CASTEP calculation.
# list - file containing the list of atom identifiers.
# transitions_file - contains transition energies for each core hole excitation.
# spins_file - file with total spin for each core hole calculation.
# nelectrons_file - file with number of electrons in total for each core hole calculation.
# neutral_file - file containing info about the neutral system.
# atom - atom for which to calculate the spectrum.
# lorentizian width (eV)
# gaussian broadening (eV)
# linear broadening (eV)
#
# Constructs names "SEED_X" where X is a number from the file "list", then 
# combines all the NEXAFS spectra for atom "atom" together and writes to a file.
#
# Uses libpytep to do the heavy lifting. Must specify a transitions file that
# adds the transition energies to each spectrum.
#
# The transitions, spins, nelectrons and neutral files are generated by transitions.sh,
# from Kane's GitHub repo under bash.
#
# UPDATE March2012: Removed the non-core hole capability. You're going to have to
# fake it by putting "neutral" into the atoms.list file and running transitions.sh.
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

DEBUG = 1

parser = argparse.ArgumentParser(description="Calculate individual NEXAFS spectra from different core-hole sites.")

#parser.add_argument('seed', help="CASTEP base seedname")
#parser.add_argument('atlist', help="List of atomic suffixes.")
#parser.add_argument('atom', help="Castep orbital number. This is hard to guess so just run pynex once and then use the number equal\nto the number of core projectors as usually the special core hole atom is last in the species.")
#parser.add_argument('transition_file', help="Transition energies file (generated by transitions.sh).")
#parser.add_argument('spins_file', help="File containing the spin of each core level calculation (generated by transitions.sh).")
#parser.add_argument('nelec_file', help="File containing the number of electrons in each core level calculation (generated by transitions.sh).")
#parser.add_argument('neutral_file', help="File containing the neutral information (generated by transitions.sh).")

# Update: Just need an input file now.
parser.add_argument('prefs', help="Input file, see pynex.py code for formatting details.")
args = parser.parse_args()

pf = open(args.prefs, 'r')
prefs = [x for x in " ".join(pf.readlines()).split()]
pf.close()

#prefix = args.seed
prefix = prefs[0]

# Open the suffixes file.
#af = open(args.atlist, 'r')
af = open(prefs[1], 'r')
atlist = [x for x in " ".join(af.readlines()).split()]
af.close()

# Get the transition energies
#tf = open(args.transition_file, 'r')
tf = open(prefs[2], 'r')
transitions = [float(x) for x in " ".join(tf.readlines()).split()]
tf.close()

# Get the spins
#sf = open(args.spins_file, 'r')
sf = open(prefs[3], 'r')
spins = [float(x) for x in " ".join(sf.readlines()).split()]
sf.close()

# Get the number of electrons
#ef = open(args.nelec_file, 'r')
ef = open(prefs[4], 'r')
nelecs = [float(x) for x in " ".join(ef.readlines()).split()]
ef.close()

# Get the neutral information
#nf = open(args.neutral_file, 'r')
nf = open(prefs[5], 'r')
neutrals = [float(x) for x in " ".join(nf.readlines()).split()]
nf.close()

#species = int(args.atom)
species = int(prefs[6])
spectra = []
rawspectra = []

for i, a in enumerate(atlist):
  print "Processing seed: ", prefix+ "_" + str(a)
  seed = prefix + "_" + str(a)
  
  # Prior to init, we need to figure out the hole charges on the spin up/spin down
  # electrons.
  # First figure out populations for the neutral system.
  nnu = 0.5 * (neutrals[2] + neutrals[1])
  nnd = nnu - neutrals[1]
  if DEBUG:
    print "For the neutral calculation, nup and ndown are: ", nnu, nnd
    print nelecs[i], spins[i]
  # Now figure out populations for the core hole system.
  cnu = 0.5 * (nelecs[i] + spins[i])
  cnd = cnu - spins[i]
  # Set these in libpytep now.
  lpe.nelectrons = array([cnu, cnd])
  if DEBUG:
    print "For atom %d, populations are: " % i, cnu, cnd
  # Core hole charges are then the difference between each.
  chu = cnu - nnu
  chd = cnd - nnd
  lpe.hole_charge = array([chu, chd])
  if DEBUG:
    print "For atom %d, the core hole charges are: " % i, chu, chd
  lpe.init(seed)
  # Set spectral preferences here!
  #lpe.lorentzian_width = 0.05
  #lpe.gaussian_broadening = 0.05
  #lpe.linear_broadening = 0.1
  #lpe.lorentzian_width = 0.05
  #lpe.gaussian_broadening = 0.05
  #lpe.linear_broadening = 0.05
  lpe.lorentzian_width = float(prefs[7])
  lpe.gaussian_broadening = float(prefs[8])
  lpe.linear_broadening = float(prefs[9])
  lpe.generate_single_spectrum(species)

  # If you use the convention that an energy difference is excited - neutral
  # in Gao's transition energy calculation, the transition energy is POSITIVE,
  # so add it here. If the transition energy is negative, subtract here instead.
  spectra.append([lpe.w.copy() + transitions[i], lpe.lastspec.copy()])
  rawspectra.append([lpe.w.copy() + transitions[i], lpe.rawspec.copy()])

# Output the single atomic spectra
for i, xy in enumerate(spectra):
  x = xy[0]
  y = xy[1]
  tmp = zeros((y.shape[0], y.shape[1]+1))
  tmp[:,0] = x
  tmp[:,1:] = y
  savetxt(prefix+ "_" + str(atlist[i]) + ".nexafs", tmp)

# Output the unbroadened atomic spectra.  
for i, xy in enumerate(rawspectra):
  x = xy[0]
  y = xy[1]
  tmp = zeros((y.shape[0], y.shape[1]+1))
  tmp[:,0] = x
  tmp[:,1:] = y
  savetxt(prefix+ "_" + str(atlist[i]) + ".raw.nexafs", tmp)

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
# Now generate interpolators and interpolate 
tmp = zeros((len(x), 7))
for xy in spectra:
  xt = xy[0]
  yt = xy[1]
  for i in range(6):
    ixy = interp1d(xt, yt[:,i], bounds_error=False, fill_value=0.0)
    tmp[:,i+1] += array([ixy(E) for E in x])
tmp[:,0] = x
savetxt(prefix + "_combined.nexafs", tmp)
    
# Do same for raw spectra
tmp = zeros((len(x), 7))
for xy in rawspectra:
  xt = xy[0]
  yt = xy[1]
  for i in range(6):
    ixy = interp1d(xt, yt[:,i], bounds_error=False, fill_value=0.0)
    tmp[:,i+1] += array([ixy(E) for E in x])
tmp[:,0] = x
savetxt(prefix + "_combined.raw.nexafs", tmp)