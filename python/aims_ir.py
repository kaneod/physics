#!/usr/bin/env python
################################################################################
#
# aims_ir.py
#
# Converts the SEED.xyz to a molden file for vibration
# viewing.
#
# Usage: aims_ir.py SEED
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
# 1. Yes, this stuff could have been written into esc_lib, but the XYZ reading
# routine is too specific here because FHI-aims uses the comment field for actual
# data. Grr. So it's been kept out.
# 
# 2. Previously needed both the SEED.xyz and the geometry.in.SEED files, now
# only need the xyz. The old code is left commented out here just in case this
# turns out to be a mistake!
#
################################################################################

import argparse
from esc_lib import reduced2cart, ang2bohr, getElementZ, elements, remove_comments
from numpy import array

DEBUG=0

parser = argparse.ArgumentParser(description="Construct a Molden input file from a FHI-aims calculation.")

parser.add_argument('seed', help="Need the SEED.xyz file for reconstruction.")
args = parser.parse_args()

xyz = open(args.seed+".xyz", 'r')
#geo = open("geometry.in."+args.seed, 'r')

#data = geo.readlines()
#geo.close()
#data = remove_comments(data, "#")
# First read the initial geometry.

#pos = []
#spec = []
#lvec = []
      
#for line in data:
#  if line.split()[0] == "atom":
#    pos.append(array([float(x) for x in line.split()[1:4]]))
#    spec.append(line.split()[4])
#  elif line.split()[0] == "lattice_vector":
#    lvec.append(array([float(x) for x in line.split()[1:4]]))
      
#if len(lvec) == 3:
#  pos = reduced2cart(pos, lvec)

# Alright, pos is now our base geometry. Now we need to read the XYZ - the comments
# have the IR frequencies and intensities.

data = xyz.readlines()
xyz.close()

ncpos = []
ncspec = []
ncdx = []
irfreq = []
irintens = []

i = 0
while i < len(data):
  if DEBUG:
    print "XYZ read: line %d" % i
  numat = int(data[i].split()[0])
  commentbits = data[i+1].split()
  irfreq.append(float(commentbits[3]))
  irintens.append(float(commentbits[8]))
  curpos = []
  curspec = []
  curdx = []
  for j in range(i+2,i+2+numat):
    linebits = data[j].split()
    curpos.append([float(x) for x in linebits[1:4]])
    curdx.append([float(x) for x in linebits[4:7]])
    curspec.append(linebits[0])
  ncpos.append(curpos)
  ncdx.append(curdx)
  ncspec.append(curspec)
  i += 2+numat

# Now we can write everything. Must convert everything to atomic units first!

f = open(args.seed+".molden", 'w')
f.write("[Molden Format]\n")
f.write("[Atoms] AU\n")
for i, p in enumerate(ncpos[0]):
  p = ang2bohr(p)
  f.write("%s\t%d\t%d\t%10.5f\t%10.5f\t%10.5f\n" % (ncspec[0][i], i+1, getElementZ(ncspec[0][i]), p[0], p[1], p[2]))

f.write("[FREQ]\n")
for frq in irfreq:
  f.write("%10.5g\n" % frq)

f.write("[FR-COORD]\n")
for s, p in zip(ncspec[0], ncpos[0]):
  p = ang2bohr(p)
  f.write("%s\t%10.5f\t%10.5f\t%10.5f\n" % (s, p[0], p[1], p[2]))

f.write("[FR-NORM-COORD]\n")
for i, alldx in enumerate(ncdx):
  f.write("vibration %d\n" % (i+1))
  for dx in alldx:
    dx = ang2bohr(dx)
    f.write("%10.5f\t%10.5f\t%10.5f\n" % (dx[0], dx[1], dx[2]))

f.write("[INT]\n")
for i in irintens:
  f.write("%10.5g\n" % i)
  
f.close()