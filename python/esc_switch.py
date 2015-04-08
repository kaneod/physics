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

if args.output_format == "aims":
  geom.writeAims(args.outputfile)
elif args.output_format == "castep":
  geom.writeCASTEP(args.outputfile)
elif args.output_format == "qe":
  geom.writeQEInput(args.outputfile)
elif args.output_format == "pdb":
  geom.writePDB(args.outputfile)
else:
  print "Error: can't recognize output file format specification. Must be castep, qe, pdb or aims. Exiting..."
  sys.exit(0)
