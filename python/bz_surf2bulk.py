#!/usr/bin/env python
################################################################################
#
# bz_surf2bulk.py
#
# Converts surface BZ coordinates back to bulk BZ lattice coordinates.
#
# Usage: load as a module or bz_surf2bulk.py num_slices slice
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
# 1. Configured for diamond unit cell.
#
################################################################################

from __future__ import division
from numpy import *
from numpy.linalg import *
import argparse

# Set up a general coordinate frame
X = array([1.0, 0.0, 0.0])
Y = array([0.0, 1.0, 0.0])
Z = array([0.0, 0.0, 1.0])

# Crystal properties of bulk *primitive* cell specified here.
a = 3.568 # angstroms, lattice constant
a1 = (a/2) * (X + Y)
a2 = (a/2) * (X + Z)
a3 = (a/2) * (Y + Z)

# Conventional cell properties specified here.
C1 = a1 + a2 - a3
C2 = a1 - a2 + a3
C3 = -a1 + a2 + a3

# Orient and size the slab here (in terms of C's or a's).
# This: 2x2 (001) unit cell. Ignore Z (slab thickness).
s1 = C1 + C2
s2 = C1 - C2
s3 = C3

# Recip space vectors (a1 -> b1, C1 -> R1, s1 -> u1 etc)
denom = dot(a1, cross(a2, a3))
b1 = 2 * pi * cross(a2, a3) / denom
b2 = 2 * pi * cross(a3, a1) / denom
b3 = 2 * pi * cross(a1, a2) / denom

denom = dot(C1, cross(C2, C3))
R1 = 2 * pi * cross(C2, C3) / denom
R2 = 2 * pi * cross(C3, C1) / denom
R3 = 2 * pi * cross(C1, C2) / denom

denom = dot(s1, cross(s2, s3))
u1 = 2 * pi * cross(s2, s3) / denom
u2 = 2 * pi * cross(s3, s1) / denom
u3 = 2 * pi * cross(s1, s2) / denom

# Points in SURFACE recip space to be converted to BULK.
pts = [ array([0.5, 0.0, 0.0]), \
        array([0.0, 0.0, 0.0]), \
        array([0.0, 0.5, 0.0])  ]

# Recip space vector (SURFACE recip space) perpendicular
# to BS path, for bulk projection.
# Note: because we left s3 = C3, the length of D here is exactly
# right for the BZ projection to the BZ edge. If we hadn't,
# we would need to scale D here. Sneaky!
D = u3

# Theory: multiplying pts by the uj expresses pts in the cartesian basis.
# That is, the matrix [u1, u2, u3] (uj as rows) converts from u-basis to cart-basis.
# Thus the inverse goes the other way. So we multiple by uj to get to the cartesian
# basis then by the inverse of the b-basis matrix to get to the b-basis.

# Construct the b matrix
bmat = matrix([b1, b2, b3])

# Need our D vector in bulk primitive recip coords.
D = asarray(dot(bmat.I, D.T)).reshape(-1)

bulk_pts = []

for p in pts:
  p_cart = p[0] * u1 + p[1] * u2 + p[2] * u3
  bulk_pts.append(asarray(dot(bmat.I, p_cart.T)).reshape(-1))

# The following output routines are only used if we run as a script.
def output_kpt_path(num_slices, slice):
  dstep = (slice / num_slices) * D
  f = open("kpt_path.txt", 'w')
  for p in bulk_pts:
    pd = p + dstep
    f.write("  %10.5f %10.5f %10.5f\n" % (pd[0], pd[1], pd[2]))
  f.close()
  
def diamond_bz_path():
  """ Returns a path that can be plotted as the diamond primitive BZ. """
    
  path = [[0.5, 0.0, 1.0],
          [0.0, -0.5, 1.0],
          [-0.5, 0.0, 1.0],
          [0.0, 0.5, 1.0],
          [0.5, 0.0, 1.0],
          [1.0, 0.0, 0.5],
          [1.0, 0.5, 0.0],
          [1.0, 0.0, -0.5],
          [1.0, -0.5, 0.0],
          [1.0, 0.0, 0.5],
          [1.0, 0.5, 0.0],
          [0.5, 1.0, 0.0],
          [0.0, 1.0, 0.5],
          [-0.5, 1.0, 0.0],
          [0.0, 1.0, -0.5],
          [0.5, 1.0, 0.0],
          [0.0, 1.0, 0.5],
          [0.0, 0.5, 1.0],
          [-0.5, 0.0, 1.0],
          [-1.0, 0.0, 0.5],
          [-1.0, -0.5, 0.0],
          [-1.0, 0.0, -0.5],
          [-1.0, 0.5, 0.0],
          [-1.0, 0.0, 0.5],
          [-1.0, 0.5, 0.0],
          [-0.5, 1.0, 0.0],
          [0.0, 1.0, -0.5],
          [0.0, 0.5, -1.0],
          [-0.5, 0.0, -1.0],
          [0.0, -0.5, -1.0],
          [0.5, 0.0, -1.0],
          [1.0, 0.0, -0.5],
          [1.0, -0.5, 0.0],
          [0.5, -1.0, 0.0],
          [0.0, -1.0, -0.5],
          [-0.5, -1.0, 0.0],
          [-1.0, -0.5, 0.0],
          [-1.0, 0.0, -0.5],
          [-0.5, 0.0, -1.0],
          [0.0, 0.5, -1.0],
          [0.5, 0.0, -1.0],
          [0.0, -0.5, -1.0],
          [0.0, -1.0, -0.5],
          [0.5, -1.0, 0.0],
          [0.0, -1.0, 0.5],
          [0.0, -0.5, 1.0],
          [0.0, -1.0, 0.5],
          [-0.5, -1.0, 0.0]]
            
  path = array(path)
    
  return path

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Generate bulk bandstructure paths from surface bandstructure path.")
  parser.add_argument('num_slices', type=int, help="Number of slices along the vector perpendicular to the surface BS path.")
  parser.add_argument('slice', type=int, help="Which particular slice you want output into the kpt_path.txt file. 0 - no displacement.")
  args = parser.parse_args()  
  ns = args.num_slices
  s = args.slice
  output_kpt_path(ns, s)
  