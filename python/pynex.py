#!/usr/bin/env python

################################################################################
#
# pynex.py
#
# Usage: pynex.py [--trans_file=file --output=[single|combined]] SEED list atom
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
from numpy import zeros, savetxt

lpe = lp.core_level_spectra

parser = argparse.ArgumentParser(description="Combine individual NEXAFS spectra from different core-hole sites.")

parser.add_argument('seed', help="CASTEP base seedname")
parser.add_argument('atlist', help="List of atomic suffixes.")
parser.add_argument('atom', help="Atom number of core hole site.")
parser.add_argument('-t', '--transition_file', default=None, help="Transition energies file.")
parser.add_argument('-o', '--output', default="single", help="Output single files or combined.")

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
  seed = prefix + "_" + str(a)
  lpe.init(seed)
  # Set spectral preferences here!
  lpe.lorentzian_width = 0.15
  lpe.gaussian_broadening = 0.3
  lpe.linear_broadening = 0.1
  lpe.generate_single_spectrum(species)
  if transitions is not None:
    spectra.append([lpe.w.copy() + transitions[i], lpe.lastspec.copy()])
  else:
    spectra.append([lpe.w.copy(), lpe.lastspec.copy()])

if args.output == "single":
  for i, xy in enumerate(spectra):
    x = xy[0]
    y = xy[1]
    tmp = zeros((y.shape[0], y.shape[1]+1))
    tmp[:,0] = x
    tmp[:,1:] = y
    savetxt(prefix+ "_" + atlist[i] + ".nexafs", tmp)

  
