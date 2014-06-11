#!/usr/bin/env python
################################################################################
#
# castep_constrain.py
#
# Usage: castep_constraints.py some.cell threshold --axis=(a, b or c)
#
# Generates constraints on all atoms below "threshold" (units fractional by default) along
# the specified axis (c axis is default). 
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
# 1.
#
################################################################################

from __future__ import division
import argparse
import os.path
import os
import sys
import esc_lib as el

el.DEBUG = 1

parser = argparse.ArgumentParser(description="Generate CASTEP constraints for a cell file.")

parser.add_argument('inputfile', help="Input .cell file.")
parser.add_argument('-t', '--threshold', help="Constrain all atoms less than this threshold (if not given, will ask).")
parser.add_argument('-a', '--axis', default="c", help="Threshold axis direction (a, b or c).")
args = parser.parse_args()

cell = el.Atoms(args.inputfile, "castep,cell")

if args.axis == "a":
  axis = 0
elif args.axis == "b":
  axis = 1
elif args.axis == "c":
  axis = 2
else:
  print "Unknown value for threshold axis - use a, b or c. Exiting..."
  sys.exit(0)

nat = len(cell.positions[0])
print "Found %d atoms in total, here given in atomic coordinates:" % (nat)
print ""
for i in range(nat):
  print "%s\t%.3g\t%.3g\t%.3g" % (el.elements[cell.species[0][i]], cell.positions[0][i][0], cell.positions[0][i][1], cell.positions[0][i][2])

# Ask user for threshold if not present
if args.threshold is None:
  threshold = float(raw_input("Enter threshold value: "))
else:
  threshold = float(args.threshold)
  
# Get all atom indices less than the appropriate threshold.
atoms = []
for i in range(nat):
  if abs(cell.positions[0][i][axis]) < threshold:
    atoms.append(i)

print ""
print "Found %d atoms to be constrained." % (len(atoms))

# Now use esc_lib's built in constraints function!

string_constraints = cell.generateCASTEPConstraints(atoms)

# Write by appending to the original file - note this means running repeatedly on the
# same cell file will result in multiple constraints blocks.
with open(args.inputfile, 'a') as f:
  f.write("\n")
  f.write(string_constraints)

