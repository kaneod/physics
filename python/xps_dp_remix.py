#!/usr/bin/env python
################################################################################
#
# xps_dp_remix.py
#
# Usage: xps_dp_remix.py casa_output.txt
#
# Turns an elemental depth profile into a materials-based one based on some
# user input.
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
# 1. Urgh, user interface. Look, it's the easy way, ok.
#
################################################################################

from __future__ import division
import argparse
import os.path
import os
import sys
from numpy import *
from numpy.linalg import *
from numpy import __version__ as npyversion

parser = argparse.ArgumentParser(description="Turn an elemental depth profile into a materials one.")

parser.add_argument('inputfile', help="Input .txt report from CasaXPS.")
args = parser.parse_args()

# First open the inputfile and grab the column headers on line 3.
f = open(args.inputfile, 'r')
headers = f.readlines()[2]
headers = headers.split('\t')[:-1]

# Now display a list of the headers to the user.
print ""
print "Welcome to XPS Depth Profile Remix!"
print ""
print "These are the headers found in the file you specified:"

for i,bit in enumerate(headers):
	print i, ":", bit

print ""

success = 0
print "At the prompt, enter the columns you want to process separated by spaces. You don't need to list the first column (etch time)."
while not success:
	columns = raw_input("Columns -> ")
	print "You've selected the following columns:"
	print columns
	user = raw_input("Is that correct? [y/n] ->")
	if user.lower() == "n" or user.lower() == "no":
		success = 0
	else:
		success = 1
	# Check for sensible input
	for column in columns.split():
		try:
			int(column)
		except ValueError:
			print "Column %s is not an integer! Try again." % (column)
			success = 0
			break
		if int(column) < 1 or int(column) > len(headers)-1:
			print "Column %s is not valid, try again." % (column)
			success = 0
			break

columns = [int(x) for x in columns.split()]
species = [headers[x].split()[0] for x in columns]
print ""
print "Ok, this is where it gets messy. Enter a list of space-separated labels you want to use for the materials you'll break the depth profile into."
success = 0
while not success:
	compounds = raw_input("Compounds -> ")
	print "You've input the following compounds:"
	print compounds
	user = raw_input("Is that correct? [y/n] ->")
	if user.lower() == "n" or user.lower() == "no":
		success = 0
	else:
		success = 1

# We now know the dimensions of our SxC matrix. Now we need to populate the matrix.
print ""
print "Enter the number of atoms in each 'unit' of the compounds, space separated again. E.g. 5 for Al2O3. I'm going to stop asking you if it's correct from now on, so watch it! Make sure you're only counting atoms you've got quantification data for in the columns you specified above."
print "Compounds ->", compounds
natoms = raw_input("Num atoms -> ")
natoms = [int(x) for x in natoms.split()]

print "Now I need to populate the compound concentration matrix with your help. For each species below, enter the number of atoms of that species per unit of the compound."
print "Compounds ->", compounds
sc_matrix = zeros((len(species), len(compounds.split())))
for i, s in enumerate(species):
	temp = raw_input("Num atoms %s -> " % (s))
	temp = [int(x) for x in temp.split()]
	for j in range(sc_matrix.shape[1]):
		sc_matrix[i,j] = temp[j] / natoms[j]

# Now we just have to invert the matrix and then apply the inverted matrix selectively to the 
# appropriate columns of the actual depth profile data.

elemental = loadtxt(args.inputfile, skiprows=3)

etch_time = elemental[:,0]
elemental_matrix = zeros((elemental.shape[0], len(columns)))
for j in range(elemental_matrix.shape[1]):
	elemental_matrix[:,j] = elemental[:,columns[j]]

sc_matrix = matrix(sc_matrix)

compound_matrix = zeros((elemental.shape[0], len(compounds.split())))
for i in range(elemental.shape[0]):
	tempS = matrix(elemental_matrix[i,:]).T
	tempC = sc_matrix.I * tempS
	compound_matrix[i,:] = array(tempC).reshape(-1,)
	# Normalize to 100% by row
	#compound_matrix[i,:] = 100 * compound_matrix[i,:] / sum(compound_matrix[i,:])
	# Note: We don't normalize here because small values of the sum make it look ridiculous.

# Done! Write to file.
output_matrix = zeros((compound_matrix.shape[0], compound_matrix.shape[1]+1))
output_matrix[:,0] = etch_time
output_matrix[:,1:] = compound_matrix
header_string = "Materials decomposition of %s done using xps_dp_remix.py (Kane O'Donnell 2014)\n\n"
header_string += "\t".join(["Etch_time"] + compounds.split())

try:
	savetxt("xps_dp_output.txt", output_matrix, header=header_string)
except TypeError:
  # Numpy version not recent enough for the header option.
  print ("""Your numpy version (%s) is not recent enough to support the\nheader option in savetxt: upgrade to a version >= 1.7 to get\nthis functionality.""" % (npyversion))
  savetxt("xps_dp_output.txt", output_matrix)

print ""
print "Done! Output written to xps_dp_output.txt. Goodbye!"
