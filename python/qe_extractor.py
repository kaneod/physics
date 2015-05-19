#!/usr/bin/env python
################################################################################
#
# qe_extractor.py
#
# Pulls all sorts of information from a QE output file and writes to standard
# output, e.g. the command "qe_extractor.py INPUTFILE homo" uses the number of 
# electrons and the output KS eigenvalues to print the KS homo.
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
# 1. A list of allowed commands appears early in the code below. Multiple commands
# leads to multiple output lines in the same order.
#
# 2. Output is simply "printed" (to stdout), so redirect to a file or a variable
# if you need the value.
#
# 3. The output isn't "safe", e.g. the code will give you a homo if you ask for 
# one even if the input is a metal. That's by design, because the target use
# for my own work considers cases where smearing has been used for an insulating
# system to help convergence (small bandgap) but the homo and lumo are needed.
# Quantum Espresso is pretty silly about this case (again, by design) and only
# reports the (non-physical) fermi level.
#
# 4. Energy outputs are in eV, because Ry and Hartree are insane.
#
# 5. The output isn't (yet) clever to geometry steps - you might get an output for
# every single SCF cycle, or you might not, depending on the command.
#
################################################################################

from __future__ import division
import argparse
import sys
import os.path
from math import floor

Ry2eV = 13.605698066
SMALL = 1.0e-6    			# small floating point number for equality comparisons
DEBUG = 1	

valid_commands = ["homo", \
                  "lumo", \
                  "num_atoms", \
                  "num_electrons", \
                  "num_bands", \
                  "num_kpoints", \
                  "total_energy", \
                  "total_ae_energy", \
                  "efermi" ]
                  
def get_eigs_from_string(str):
  """ This is string.split() tweaked to address a bug in Quantum Espresso's formatted fortran
  output where sometimes it prints two negative floats without a space e.g. -130.3940-120.6023.
  In these cases, split() doesn't work directly. """
  
  eigs = []
  bits = str.split()
  for b in bits:
    if b is not '':
      try:
        tmpf = float(b)
        eigs.append(tmpf)
      except ValueError: 
        negs = b.split('-')
        eigs += [-1 * float(c) for c in negs if c is not ''] # This is a bit obscure I know...
  
  return eigs

parser = argparse.ArgumentParser(description="Extract information from QE(PWSCF) output file and print to stdout.")

parser.add_argument('inputfile', help="Quantum Espresso pw.x output file for input.")
parser.add_argument('commands', nargs="+", help="Parameters to be extracted.")

args = parser.parse_args()		

# Check we have valid commands

for c in args.commands:
  if c not in valid_commands:
    print "ERROR: command %s is not valid, see source file for a list of possible commands." % (c)

# Some of the commands are easy, others require more complex parsing. Deal with all the easy ones
# first.

f = open(args.inputfile, 'r')
lines = f.readlines()
f.close()

found_fermi = False
found_homo = False
found_lumo = False
found_ae_energy = False

output_text = {}

for l in lines:
  if "number of atoms/cell      =" in l:
    if "num_atoms" in args.commands:
      output_text["num_atoms"] = l.split()[4]
  if "number of electrons       =" in l:
    if "num_electrons" in args.commands:
      output_text["num_electrons"] = l.split()[4]
    nelec = float(l.split()[4])
  if "number of Kohn-Sham states=" in l:
    if "num_bands" in args.commands:
      output_text["num_bands"] = l.split()[4]
    nband = int(l.split()[4])
  if "number of k points=" in l:
    if "num_kpoints" in args.commands:
      output_text["num_kpoints"] = l.split()[4]
    nkpt = int(l.split()[4])
  if "!    total energy              =" in l:
    if "total_energy" in args.commands:
      output_text["total_energy"] = float(l.split()[4]) * Ry2eV
  if "total all-electron energy =" in l:
    if "total_ae_energy" in args.commands:
      output_text["total_ae_energy"] = float(l.split()[4]) * Ry2eV
      found_ae_energy = True
  if "highest occupied, lowest unoccupied level (ev):" in l:
    # This might not be present - opportunistic!
    homo = float(l.split()[6])
    found_homo = True
    lumo = float(l.split()[7])
    found_lumo = True
    if "homo" in args.commands:
      output_text["homo"] = homo
    if "lumo" in args.commands:
      output_text["lumo"] = lumo
  if "the Fermi energy is" in l:
    efermi = float(l.split()[4])
    found_fermi = True
    if "efermi" in args.commands:
      output_text["efermi"] = efermi
  if "highest occupied level (ev):" in l:
    homo = float(l.split()[4])
    found_homo = True
    if "homo" in args.commands:
      output_text["homo"] = homo


# If PAW potentials aren't used, an all-electron energy won't be reported so if the user
# asked for one, give an error.
if "total_ae_energy" in args.commands and not found_ae_energy:
  print "ERROR - All-electron energy not found. Check calculation used PAW and that it finished correctly."

# Ok, now for the slightly trickier ones - homo, lumo and fermi_level. First, QE might actually
# give us values, which we picked up earlier. If not, we need to do a bit more work.

if ("homo" in args.commands and found_homo is False) or \
  ("lumo" in args.commands and found_lumo is False) or \
  ("efermi" in args.commands and found_fermi is False):
    
  # Lots of things to worry about here and we have to loop a lot. For performance, find
  # the important section of the file.
  for i,l in enumerate(lines):
    if "End of self-consistent calculation" in l:
      istart = i
    if "convergence has been achieved in" in l:
      iend = i
  
  has_spin = False
  for i,l in enumerate(lines[istart:iend]):
    if "SPIN UP" in l:
      has_spin = True
      istart = i
    if "SPIN DOWN" in l:
      iend = i
  
  # Find the k-point block locations
  ks = []
  for i,l in enumerate(lines[istart:iend]):
    if "k =" in l:
      if DEBUG:
        print l
      ks.append(i+istart + 1)
  
  if DEBUG:
    print "K-point indices are:"
    print ks
  
  # Add the iend value to act as an endpoint for the 
  # eigenvalue search.
  
  ks.append(iend)
  
  # For each k, look for eigenvalues until we have enough.
  
  eigsk = []
  
  for i in range(len(ks)-1):
    eigs = []
    for l in lines[ks[i]+1:ks[i+1]]:
      # Now - pay attention! There is a bug in the output of espresso that means split() might
      # not work here. This means we have to play a silly game here assuming eigenvalues are
      # output in increasing order.
      # This is done with a recursive function defined at the top of the file.
      eigs += get_eigs_from_string(l)
      if DEBUG:
        print "Current length of eigs is %d, num_bands is %d." %(len(eigs), nband)
      if len(eigs) == nband:
        eigsk.append(eigs)
        break
  
  # Now, use number of electrons to figure out where the homo is.
  max_occ = -1e6
  min_unocc = 1e6
  idx_homo = int(floor(nelec / 2)) - 1 # The -1 is because we have 0-based indices in python.

  for ek in eigsk:
    if ek[idx_homo] > max_occ:
      max_occ = ek[idx_homo]
    if ek[idx_homo + 1] < min_unocc:
      min_unocc = ek[idx_homo + 1]
  
  if not found_homo:
    homo = max_occ
    output_text["homo"] = homo
  if not found_lumo:
    lumo = min_unocc
    output_text["lumo"] = lumo
  if not found_fermi:
    efermi = (homo + lumo) / 2
    output_text["efermi"] = efermi
  
# Print output in the order requested.
for c in args.commands:
  print output_text[c]    
  