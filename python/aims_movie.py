#!/usr/bin/env python
################################################################################
#
# aims_movie.py
#
# Makes an XYZ movie from FHI-aims output file.
#
################################################################################
#
# Copyright 2015 Kane O'Donnell
#
#    This library is free software: you can redistribute it and/or modify
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
# 1. Yes, create_xyz_movie.pl exists already in the aims utilities. However,
# it doesn't seem to work for periodic geometries, and it's easier for me to
# rewrite in Python than correct the Perl. 
#
# 2. Output goes to terminal just like create_xyz_movie.pl - redirect to a file
# as necessary.
#
################################################################################

import argparse
import sys

parser = argparse.ArgumentParser(description="Create XYZ movie from FHI-aims output.")

parser.add_argument('inputfile', help="Input FHI-aims output file.")
args = parser.parse_args()

try:
  f = open(args.inputfile, 'r')
  lines = f.readlines()
  f.close()
except IOError:
  print "ERROR: Failed to open file %s - check file name and try again." % (args.inputfile)
  sys.exit(0)
  
num_atoms = -1
# First get the number of atoms
for l in lines:
  if "| Number of atoms" in l:
    num_atoms = int(l.split()[5])
    break

if num_atoms == -1:
  print "ERROR: Number of atoms could not be read, are you sure you specified a FHI-aims output?"
  sys.exit(0)
else:
  print num_atoms

# Now get the input geometry 
for i, l in enumerate(lines):
  if "Input geometry:" in l:
    start_block = i
    break

if "Unit cell:" in lines[start_block + 1]:
  start_block += 7
  has_unit_cell = True
else:
  start_block += 4
  has_unit_cell = False

print 0
for l in lines[start_block:start_block + num_atoms]:
  print l.split()[3], l.split()[4], l.split()[5], l.split()[6]

# Next pull the geometry updates.
updates = []
for i, l in enumerate(lines):
  if "Updated atomic structure:" in l:
    updates.append(i)

step = 1
for i in updates:
  print num_atoms
  print step
  if has_unit_cell:
    start_block = i + 6
  else:
    start_block = i + 2
  for l in lines[start_block:start_block + num_atoms]:
    print l.split()[4], l.split()[1], l.split()[2], l.split()[3]
  step += 1

# Finally get the final atomic structure (if it's there).
start_block = -1
for i, l in enumerate(lines):
  if "Final atomic structure:" in l:
    start_block = i
    break

if start_block > -1:

  if has_unit_cell:
    start_block += 6
  else:
    start_block += 2

  print num_atoms
  print step
  for l in lines[start_block:start_block + num_atoms]:
    print l.split()[4], l.split()[1], l.split()[2], l.split()[3]
