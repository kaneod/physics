#!/usr/bin/env python
################################################################################
#
# mix_bulk_bands.py
#
# Takes a set of SEED_X.bands_xy files and reorders so that SEED_b0.bands_xy contains
# all the first bands, SEED_b1.bands_xy contains all the second bands and so on.
#
# Usage: mix_bulk_bands.py SEED num_files
#
################################################################################
#
# Copyright 2013 Kane O'Donnell
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
# 1. Assumes the existence of files continuously from 0 to num_files-1. Error will
# occur if one of the files is missing.
#
################################################################################

from numpy import savetxt, array, loadtxt, matrix
import argparse

parser = argparse.ArgumentParser(description="Reorders a set of SEED_X.bands_xy files into individual bands.")
parser.add_argument('seed', help="The .bands_xy file seed to convert.")
parser.add_argument('num_files', type=int, help="Number of files to convert.")
args = parser.parse_args()  

banddata = []
for i in range(args.num_files):
  banddata.append(loadtxt("%s_%d.bands_xy" % (args.seed, i)))

# Reorder the bands into [band0], [band1], etc.
bandsets = []
for i in range(banddata[0].shape[1]):
  tmpbands = []
  for bd in banddata:
    tmpbands.append(bd[:,i])
  tmpbands = matrix(tmpbands).T
  bandsets.append(tmpbands)

for i, bs in enumerate(bandsets):
  savetxt("%s_b%d.bands_xy" % (args.seed, i), bs)


  

