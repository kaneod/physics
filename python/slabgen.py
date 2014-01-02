#!/usr/bin/env python
################################################################################
#
# slabgen.py
#
# Generates a diamond slab with specified surface orientation and number of layers.
# Can optionally hydrogen-terminate atoms to make a single-sided slab or a cluster.
#
# Usage: slabgen.py [-h -c -d] -n ABC -s X Y Z
#
# Options: -h hydrogen-terminates everything that isn't The Surface or periodic wall.
# -c specifies a cluster. -d specifies double-sided.
#
# Essentials: -n is the surface orientation and -s is the slab size in multiples of the
# unit cell parameter. 
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
from numpy import *
from numpy.linalg import *
from random import choice

import argparse
import sys

# Constants
TOL = 1.0e-6  # Two vectors are the same if the norm of their difference is less than this.

# Euclidean Vectors
Ex = array([1.0, 0.0, 0.0])
Ey = array([0.0, 1.0, 0.0])
Ez = array([0.0, 0.0, 1.0])

# Lattice vectors
# (in conventional unit cell!)
a = 3.568
Lx = a * Ex
Ly = a * Ey
Lz = a * Ez

# (primitive cell)
Px = 0.5 * Lx + 0.5 * Ly
Py = 0.5 * Ly + 0.5 * Lz
Pz = 0.5 * Lx + 0.5 * Lz

# Because the lattice vectors are not necessarily orthogonal we need an overlap matrix
# (more specifically, the inverse)

Lov = matrix([[dot(Lx, Lx), dot(Lx, Ly), dot(Lx, Lz)],
              [dot(Ly, Lx), dot(Ly, Ly), dot(Ly, Lz)],
              [dot(Lz, Lx), dot(Lz, Ly), dot(Lz, Lz)]])
              
Pov = matrix([[dot(Px, Px), dot(Px, Py), dot(Px, Pz)],
              [dot(Py, Px), dot(Py, Py), dot(Py, Pz)],
              [dot(Pz, Px), dot(Pz, Py), dot(Pz, Pz)]])

# Atomic basis (8 atoms)
basis = array([[0.0,0.0,0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5], \
               [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], \
               [0.25, 0.75, 0.75]])
species = ["C"] * 8
valence = [4] * 8

# Functions
def radialDistance(vectors, index):
  """ radialDistance(vectors, index)
  
  """
  
  norms = []
  for i, v in enumerate(vectors):
    if i != index:
      norms.append(norm(v - vectors[index]))
  
  return norms
  
def uniqueFloats(floats, precision):
  """ returns a dict of the unique floats where two floats are the same if they are the 
  same after rounding to the given number of decimal places.
  """
  unique_floats = {}
  rounded = around(floats, decimals=precision)
  for f in rounded:
    if f in unique_floats.keys():
      unique_floats[f] += 1
    else: 
      unique_floats[f] = 1
  
  return unique_floats
  
def normsort(vectors):
  """ normsort(vectors)
  
  Sorts a list of vectors by norm in descending order. Uses quicksort algorithm.
  
  """
  
  if len(vectors) <= 1:
    return vectors
  else:
    less = []
    more = []
    pivot = choice(vectors)
    same = [pivot]
    np = norm(pivot)
    for v in vectors:
      nv = norm(v)
      if nv < np:
        less.append(v)
      if nv > np:
        more.append(v)
      if nv == np:
        same.append(v)
    less = normsort(less)
    more = normsort(more)
    return more + same + less
  
# Main program
if __name__ == '__main__':

  print "slabgen version 0.1"
  print ""
  print "Written by Kane O'Donnell, December 2013"
  print ""
  
  # Parse the command line arguments
  parser = argparse.ArgumentParser(description="Generate crystal slab.")
  parser.add_argument("--surface", "-s", type=int, nargs=3, default=[0, 0, 1], help="Lattice vector multipliers (3 ints) giving the surface normal direction.")
  parser.add_argument("--inplane", "-p", type=float, nargs=3, default=[0.5, 0.5, 0.0], help="Lattice vector multipliers (3 floats) giving an in-plane surface vector.")
  parser.add_argument("--dimensions", "-d", type=float, nargs=2, default=[10.0, 10.0], help="Requested cell x,y dimensions in angstrom. The actual size will be as close to this as possible without breaking lattice periodicity, unless cluster mode is activated in which case the dimensions are exact.")
  parser.add_argument("--offset", "-o", type=float, nargs=3, default=[0.0, 0.0, 0.0], help="Add an offset to the atomic positions. Make this small relative to the unit cell dimensions or all hell will break loose.")
  parser.add_argument("--thickness", "-z", type=float, default=0.0, help="Thickness of slab in angstrom. If not set, the length of the surface orientation vector is used.")
  parser.add_argument("--cluster", "-c", action="store_true", help="If set, all unsaturated atoms are hydrogen-terminated and no periodicity checking is done.")
  #parser.add_argument('Z', help="El")
  #parser.add_argument('n', type=int, help="Make the nth atom of element Z special. 1-based.")
  args = parser.parse_args()

  n = array(args.surface)
  p = array(args.inplane)
  cell_size = array(args.dimensions)
  is_cluster = args.cluster
  Voff = array(args.offset)
    
  
  # Generate the surface normal vector and the in-plane vector. Check the in-plane 
  # vector is actually in-plane!
  
  Vn = n[0] * Lx + n[1] * Ly + n[2] * Lz
  Vp1 = p[0] * Lx + p[1] * Ly + p[2] * Lz
  
  if dot(Vn, Vp1) != 0.0:
    print "Error: the in-plane vector is not perpendicular to the surface normal vector!"
    sys.exit(0)
    
  # Take the cross product of the two vectors to get another in-plane vector. Thus Vp1,
  # Vp2 and Vn are the lattice vectors for the new unit cell.
    
  Vp2 = cross(Vn, Vp1)
  
  print "Original lattice:", Lx, Ly, Lz
  print ""
  print "Requested lattice:", Vp1, Vp2, Vn
  print ""
  
  # Generate our first guess at the cell based solely on the user input. 
  
  Sx = cell_size[0] * Vp1 / norm(Vp1)
  Sy = cell_size[1] * Vp2 / norm(Vp2)
  if args.thickness > 0.0:
    Sz = args.thickness * Vn / norm(Vn)
  else:
    Sz = Vn
  
  # Now, we need to check that Sx and Sy are integer
  # linear combinations of (primitive) lattice vectors. Because the lattice vectors don't 
  # need to be orthogonal, we need to use the Pxyz matrix inversion on each vector. Note 
  # that we already know Sz is an integer multiple as we force it that way in the input.

  print "Checking requested lattice vectors match periodicity..."
  print ""
  newS = []
  for s in [Sx, Sy]:
    x = array([dot(Px, s), dot(Py, s), dot(Pz, s)])
    lmn = dot(Pov.I, x)
    print "For s", s, "lmn is", lmn
    lmn = array(lmn).flatten()
    lmn = around(lmn)
    print "Requested vector ", s, "changed to ", lmn[0] * Px + lmn[1] * Py + lmn[2] * Pz
    newS.append(lmn[0] * Px + lmn[1] * Py + lmn[2] * Pz)
  Sx = newS[0]
  Sy = newS[1]
  
  # Also need unit vectors...
  Ux = Sx / norm(Sx)
  Uy = Sy / norm(Sy)
  Uz = Sz / norm(Sz)
  
  # If we are making a cluster, trim the cell to the right XY dimensions.
  if is_cluster:
    Sx = cell_size[0] * Ux
    Sy = cell_size[1] * Uy
    
  # Get vector norms...
  nx = norm(Sx)
  ny = norm(Sy)
  nz = norm(Sz)
  
  # Cell central point
  Sc = 0.5 * Sx + 0.5 * Sy + 0.5 * Sz
  
  # Project each of the vectors onto each of the lattice vectors to get the min and max
  # integer displacements of the unit cell required.
  
  minx = floor(min(dot(Sx, Ex)/norm(Lx), dot(Sy, Ex)/norm(Lx), dot(Sz, Ex)/norm(Lx))) - 1
  maxx = ceil(max(dot(Sx, Ex)/norm(Lx), dot(Sy, Ex)/norm(Lx), dot(Sz, Ex)/norm(Lx))) + 1
  miny = floor(min(dot(Sx, Ey)/norm(Ly), dot(Sy, Ey)/norm(Ly), dot(Sz, Ey)/norm(Ly))) - 1
  maxy = ceil(max(dot(Sx, Ey)/norm(Ly), dot(Sy, Ey)/norm(Ly), dot(Sz, Ey)/norm(Ly))) + 1
  minz = floor(min(dot(Sx, Ez)/norm(Lz), dot(Sy, Ez)/norm(Lz), dot(Sz, Ez)/norm(Lz))) - 1
  maxz = ceil(max(dot(Sx, Ez)/norm(Lz), dot(Sy, Ez)/norm(Lz), dot(Sz, Ez)/norm(Lz))) + 1
  
  # Give some output
  print "Actual lattice: ", Sx, Sy, Sz
  print ""
  print "Min and max i: ", minx, maxx
  print "Min and max j: ", miny, maxy
  print "Min and max k: ", minz, maxz
  
  positions = []
  for i in arange(minx, maxx+1):
    for j in arange(miny, maxy+1):
      for k in arange(minz, maxz+1):
        # Displacement vector for these indices
        Vd = i * Lx + j * Ly + k * Lz
        # For each atomic coordinate A in the basis, see if Vd + A is inside the 
        # rhombohedron defined by the cell basis.
        for b in basis:
          tmp = b[0] * Lx + b[1] * Ly + b[2] * Lz + + Voff + Vd
          ax = abs(dot(Sc - tmp, Ux))
          ay = abs(dot(Sc - tmp, Uy))
          az = abs(dot(Sc - tmp, Uz))
          # Decide if tmp is inside the cube. Don't worry about periodic duplicates, we
          # sort these out later.
          if ax - 0.5 * nx <= TOL and ay - 0.5 * ny <= TOL and az - 0.5 * nz <= TOL:
              positions.append(tmp)
  
  # We sort the positions by distance from the origin. This is to make the uniqify routine
  # (which happens next) always drop the farthest lattice duplicates rather than the
  # closer ones.
  
  positions = normsort(positions)
  
  # Now let's go through each pair of positions and check for any that differ by a lattice
  # vector (eliminate the duplicates). 
  
  unique_positions = []
  is_unique = False
  for i in range(len(positions)):
    is_unique = True
    for j in range(i+1, len(positions)):
      dv = positions[i] - positions[j]
      if norm(dv) < TOL:
        #print "Found identical atoms:"
        #print "Details: atom %d is the same as atom %d" % (i, j)
        is_unique = False
        break
      elif norm(dv - Sx) < TOL or norm(dv - Sy) < TOL or norm(dv - Sz) < TOL:
        #print "Found periodicity-connected pair!"
        #print "Details: ", i, j, positions[i], positions[j]
        is_unique = False
        break
      elif norm(dv + Sx) < TOL or norm(dv + Sy) < TOL or norm(dv + Sz) < TOL:
        #print "Found periodicity-connected pair!"
        #print "Details: ", i, j, positions[i], positions[j]
        is_unique = False
        break
    if is_unique:
      unique_positions.append(positions[i])
 
  # If we want to terminate, do it here.
  # Not done yet!
  
  floats = radialDistance(unique_positions, 35)
  ufloats = uniqueFloats(floats, 4)
  for k in sorted(ufloats.keys()):
    print "%.4g %d" % (k, ufloats[k])
    
  #for p in unique_positions:
  #  print norm(p)
  
  # Write a quick output file. If cluster, XYZ, if not XSF.
  if is_cluster:
    f = open("test.xyz", 'w')
    f.write(str(len(unique_positions)) + "\n")
    f.write("TEST OUTPUT\n")
    for p in unique_positions:
      f.write("C    %.6g    %.6g    %.6g\n" % (p[0], p[1], p[2]))
    f.close()
  else:
    f = open("test.xsf", 'w')
    f.write("CRYSTAL\n")
    f.write("PRIMVEC\n")
    for v in [Sx, Sy, Sz]:
      f.write("%.6g    %.6g    %.6g\n" % (v[0], v[1], v[2]))
    f.write("CONVVEC\n")
    for v in [Lx, Ly, Lz]:
      f.write("%.6g    %.6g    %.6g\n" % (v[0], v[1], v[2]))
    f.write("PRIMCOORD\n")
    f.write("%d    1\n" % (len(unique_positions)))  
    for p in unique_positions:
      f.write("6    %.6g    %.6g    %.6g\n" % (p[0], p[1], p[2]))
    f.close()