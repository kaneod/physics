#!/usr/bin/env python
################################################################################
#
# esc_radius_relabel.py
#
# Given a radius in angstroms and an atomic index, adds a suffix to all atoms
# within range of that specified atom. Handles periodic boundary conditions.
#
################################################################################
#
# Copyright 2015 Kane O'Donnell
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
from numpy.linalg import norm
import argparse
import sys
import esc_tools as et

DEBUG=1

parser = argparse.ArgumentParser(description="Add a suffix to all atoms within range of a specified atom.")

parser.add_argument('atom_index', type=int, help="1-based atomic index of atom focus.")
parser.add_argument('radius', type=float, default=2.0, help="Radius to relabel (angstroms, default 2.0).")
parser.add_argument('input_format', help="Type of input (aims, castep, qe or pdb).")
parser.add_argument('input_file', help="Input filename.")
args = parser.parse_args()

if args.input_format == "aims":
  geom = et.Atoms(args.input_file, "aims,geometry")
elif args.input_format == "castep":
  geom = et.Atoms(args.input_file, "castep,cell")
elif args.input_format == "qe":
  geom = et.Atoms(args.input_file, "qe,input")
elif args.input_format == "pdb":
  geom = et.Atoms(args.input_file, "pdb")
else:
  print "Error: can't recognize input file format specification. Must be castep, qe, pdb or aims. Exiting..."
  sys.exit(0)
  
# Case 1: no unit cell. Doesn't work for castep, of course. 

if geom.lattice is None:
  special_p = geom.positions[args.atom_index - 1]
  index_list = []
  for i, p in enumerate(geom.positions):
    if (i != args.atom_index - 1) and norm(special_p - p) < args.radius:
      index_list.append(i)
  if DEBUG:
    print "No PBCs, index atom is:", args.atom_index - 1
    print "List is: ", index_list  
# Case 2: unit cell present, periodic boundary conditions need to be accounted for.
else:
  special_p = geom.positions[args.atom_index - 1]
  index_list = []
  g_list = [-1, 0, 1]
  g1 = geom.lattice[0]
  g2 = geom.lattice[1]
  g3 = geom.lattice[2]
  for i, p in enumerate(geom.positions):
    if (i != args.atom_index - 1):
      for a in g_list:
        for b in g_list:
          for c in g_list:
            ppg = special_p + a * g1 + b * g2 + c * g3
            if norm(ppg - p) < args.radius:
              index_list.append(i)
  if DEBUG:
    print "PBCs present, index atom is:", args.atom_index - 1
    print "List is: ", index_list  

# Write out in the same format as read.

ch_opt = {"corehole" : args.atom_index - 1, "suffix atoms" : index_list}

if args.input_format == "aims":
  geom.writeAims(args.input_file, opt=ch_opt)
elif args.input_format == "castep":
  geom.writeCASTEP(args.input_file, opt=ch_opt)
elif args.input_format == "qe":
  geom.writeQEInput(args.input_file, opt=ch_opt)
elif args.input_format == "pdb":
  geom.writePDB(args.input_file, opt=ch_opt)
else:
  print "Error: Output file format unrecognized - this should be impossible at this stage of the program! Someone has messed with the code. Exiting..."
  sys.exit(0)
      