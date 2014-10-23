#!/usr/bin/env python
################################################################################
#
# qnex2nexafs.py
#
# Version 1.0
#
# Usage qnex2nexafs.py SEED
#
# Converts the series of xspectra outputs (in the form of .qnex files from Kane's
# BASH script) to the .nexafs format created by nexspec, for use with the Igor PRO
# nexspec loader.
#
################################################################################
#
# Copyright 2014 Kane O'Donnell
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
# 1. 
#
################################################################################

from __future__ import division
from numpy import *
import argparse
import sys

# Parse the command line arguments
parser = argparse.ArgumentParser(description="Collect XSpectra calculations into .nexafs format.")
parser.add_argument('SEED', help="Input files have the format xanes_SEED_vv.qnex")
args = parser.parse_args()

cmpts = {1:"xx", 2:"yy", 3:"zz", 4:"xy", 5:"xz", 6:"yz"}
prefix = "xanes"
extension = ".qnex"

f = open(args.SEED+".nexafs", 'w')
f.write("# NEXAFS core-level spectrum calculated by xspectra with Quantum Espresso input.\n")
f.write("# Spin channel:            1\n")
f.write("# Omega (eV) Mxx Myy Mzz Mxy Mxz Myz\n")

M = None
for i in cmpts.keys():
  print "Trying filename: " + prefix + "_" + args.SEED + "_" + cmpts[i] + extension
  tmp = loadtxt(prefix + "_" + args.SEED + "_" + cmpts[i] + extension)
  if tmp is None:
    print "ERROR: couldn't open the " + cmpts[i] + " file. Exiting."
    sys.exit(0)
  if M is None:
    M = zeros((tmp.shape[0], 7))
    M[:,0] = tmp[:,0]
  M[:,i] = tmp[:,2]
  
# Use the angle terms (xy) to deduce the matrix cross terms.
# (MxMby + MbxMy) = 2*spectrum(xy) - spectrum(xx) - spectrum(yy), etc.
M[:,4] = 2 * M[:,4] - M[:,1] - M[:,2]
M[:,5] = 2 * M[:,5] - M[:,1] - M[:,3]
M[:,6] = 2 * M[:,6] - M[:,2] - M[:,3]

for i in range(M.shape[0]):
  f.write("%g %g %g %g %g %g %g\n" % (M[i,0], M[i,1], M[i,2], M[i,3], M[i,4], M[i,5], M[i,6]))

f.close()
  
  