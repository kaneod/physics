#!/usr/bin/env python
################################################################################
#
# esc_switch.py
#
# Convenience tool to switch files between QE, FHI-aims, CASTEP and pdb format.
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
# 1. Uses esc_tools.
#
################################################################################

from __future__ import division
import argparse
import sys
import esc_tools as et

DEBUG=0

parser = argparse.ArgumentParser(description="Switch between electronic structure file formats.")

parser.add_argument('-i', '--input_format', default='pdb', help="File format of input - aims, castep, qe or pdb.")
parser.add_argument('-o', '--output_format', default='castep', help="File format for output - aims, castep, qe or pdb.")
parser.add_argument('-c', '--constrain_hydrogens', default=False, dest='constrain_hydrogens', action='store_true', help="Constrain hydrogen atoms.")
parser.add_argument('-t', '--constrain_below_threshold', type=float, default=-1.0, help="Constrain all atoms below the given z value.")
parser.add_argument('-x', '--charge', type=float, default=0.0, help="Put a charge on the resulting cell (QE only).")
parser.add_argument('-p', '--prefix', default="esc", help="Labels the calculation with a name (QE only).")
parser.add_argument('-b', '--bands', type=int, default=0, help="States the number of bands for the calculation (QE only).")
parser.add_argument('inputfile', help="Input file.")
parser.add_argument('outputfile', help="Output file.")
args = parser.parse_args()

if args.input_format == "aims":
  geom = et.Atoms(args.inputfile, "aims,geometry")
elif args.input_format == "castep":
  geom = et.Atoms(args.inputfile, "castep,cell")
elif args.input_format == "qe":
  geom = et.Atoms(args.inputfile, "qe,input")
elif args.input_format == "pdb":
  geom = et.Atoms(args.inputfile, "pdb")
else:
  print "Error: can't recognize input file format specification. Must be castep, qe, pdb or aims. Exiting..."
  sys.exit(0)

opts = {}
# Construct options dictionary
if args.constrain_hydrogens is True:
	opts["constrain species"] = ["H"]

if args.constrain_below_threshold > 0:
  # Search for all atoms below the threshold.
  atoms = []
  for i in range(len(geom.positions)): 
    if geom.positions[i][2] < args.constrain_below_threshold:
      atoms.append(i)
  opts["constrain atoms"] = atoms

if args.charge != 0.0:
  opts["charge"] = args.charge

opts["prefix"] = args.prefix

if args.bands != 0:
  opts["bands"] = args.bands
	
if args.output_format == "aims":
  geom.writeAims(args.outputfile, opt=opts)
elif args.output_format == "castep":
  geom.writeCASTEP(args.outputfile, opt=opts)
elif args.output_format == "qe":
  geom.writeQEInput(args.outputfile, opt=opts)
elif args.output_format == "pdb":
  geom.writePDB(args.outputfile)
else:
  print "Error: can't recognize output file format specification. Must be castep, qe, pdb or aims. Exiting..."
  sys.exit(0)
