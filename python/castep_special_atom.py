#!/usr/bin/env python
################################################################################
#
# castep_special_atom.py
#
# Version 2.0
#
# Usage castep_special_atom.py SEED Z n
#
# Takes a .cell file SEED.cell, reads it in, changes the "nth" atom of type Z
# to Zspec, writes the file to SEED_special.cell with otherwise no changes.
#
# Designed for calculations where you have to repeat a calculation over all atoms
# of type Z with minor modifications to the pseudopotential for one special atom.
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
# 1. WHY DID I NOT USE ESC_LIB? Because esc_lib.py has scipy dependencies that
# make it hard to get running on a cluster where the default python is from the
# dark ages and installing scipy is a pain anyway.
#
# 2. WHY DID I NOT JUST USE THE SAME DAMN ALGORITHMS AS ESC_LIB? Because it's 
# a couple of years later and I have better ways of doing things now.
#
# 3. I did reuse SOME bits of esc_lib, ok? :o)
#
################################################################################

# This will fail on old python versions but is fine for 2.6 and 2.7
from __future__ import division
import argparse

# Debugging flag - set to 1 to see debug messages.
DEBUG=1

# SPECIAL text - what do you want added after the atomic Z to show it as special?
SPECIAL=":exc"

# Element dictionary - note that element 0 is given the UKN moniker.

elements = { 1 : "H", 2 : "He", 3 : "Li", 4 : "Be", 5 : "B", 6 : "C", 7 : "N", \
            8 : "O", 9 : "F", 10 : "Ne", 11 : "Na", 12 : "Mg", 13 : "Al", 14 : "Si", \
            15 : "P", 16 : "S", 17 : "Cl", 18 : "Ar", 19 : "K", 20 : "Ca", 21 : "Sc", \
            22 : "Ti", 23 : "V", 24 : "Cr", 25 : "Mn", 26 : "Fe", 27 : "Co", 28 : "Ni", \
            29 : "Cu", 30 : "Zn", 31 : "Ga", 32 : "Ge", 33 : "As", 34 : "Se", \
            35 : "Br", 36 : "Kr", 37 : "Rb", 38 : "Sr", 39 : "Y", 40 : "Zr", 41 : "Nb", \
            42 : "Mo", 43 : "Tc", 44 : "Ru", 45 : "Rh", 46 : "Pd", 47 : "Ag", \
            48 : "Cd", 49 : "In", 50 : "Sn", 51 : "Sb", 52 : "Te", 53 : "I", 54 : "Xe", \
            55 : "Cs", 56 : "Ba", 57 : "La", 58 : "Ce", 59 : "Pr", 60 : "Nd", \
            61 : "Pm", 62 : "Sm", 63 : "Eu", 64 : "Gd", 65 : "Tb", 66 : "Dy", \
            67 : "Ho", 68 : "Er", 69 : "Tm", 70 : "Yb", 71 : "Lu", \
            72 : "Hf", 73 : "Ta", 74 : "W", 75 : "Re", 76 : "Os", 77 : "Ir", 78 : "Pt", \
            79 : "Au", 80 : "Hg", 81 : "Tl", 82 : "Pb", 83 : "Bi", 84 : "Po", \
            85 : "At", 86 : "Rn", 87 : "Fr", 88 : "Ra", 89 : "Ac", 90 : "Th", \
            91 : "Pa", 92 : "U", 93 : "Np", 94 : "Pu", 95 : "Am", 96 : "Cm", 97 : "Bk", \
            98 : "Cf", 99 : "Es", 100 : "Fm", 101 : "Md", 102 : "No", 103 : "Lr", \
            104 : "Rf", 105 : "Db", 106 : "Sg", 107 : "Bh", 108 : "Hs", 109 : "Ds", \
            110 : "Ds", 111 : "Rg", 112 : "Uub", 113 : "Uut", 114 : "Uuq", 115 : "Uup", \
            116 : "Uuh", 117 : "Uus", 118 : "Uuo", 0 : "UKN"}
            

def getElementZ(elstr):
    """ Z = getElementZ(elstr)
    
    Given a string that contains either a Z number OR an element
    abbreviation like Cu, MG, whatever, generates and returns the
    appropriate integer Z.
    """
    
    # Is it an integer?
    try:
        Z = int(elstr)
        return Z
    except ValueError:
        # Not an integer.
        if elstr.title() not in elements.values():
            print "(libesc.getElementZ) Warning: Element %s is not in the elements dictionary. Returning 0 for element Z." % elstr
            return 0
        else:
            for key, value in elements.items():
                if elstr.title() == value:
                    return key

def remove_comments(lines, comment_delim="#",just_blanks=False):
    """ stripped = remove_comments(lines, comment_delim="#", just_blanks=False)
    
    Takes a sequence of lines (presumably from a data file)
    and strips all comments, including ones at the end of
    lines that are otherwise not comments. Note that we can
    only deal with one kind of comment at a time - just apply
    multiple times to strip multiple comment types.
    
    Note that we also remove *blank* lines in here, just in
    case we're going to do line-by-line processing subsequently
    rather than join/splitting (eg. like for XSF files).
    
    If just_blanks is specified, we only eliminate the blanks.
    
    """
    
    stripped = []
    for line in lines:
      if just_blanks:
        if line.strip() != "":
          stripped.append(line.strip())
      else:
        if not (line.strip().startswith(comment_delim) or line.strip() == ""):
          stripped.append(line.partition(comment_delim)[0].strip())
    
    return stripped
    
def uniqify(sequence, trans=None):
    """ unique = uniqify(sequence, trans)
    
    Produces an order-preserved list of unique elements in the passed
    sequence. Supports a transform function so that passed elements
    can be transformed before comparison if necessary.
    
    """
    
    if trans is None:
        def trans(x): return x
    seen = {}
    unique = []
    for item in sequence:
        marker = trans(item)
        if marker in seen: continue
        seen[marker] = 1
        unique.append(item)
    return unique

def identify_block(line):
  """ blockstr = identify_block(line)
  
  Identifies the actual block title of a castep block like "%block positions_frac".
  Can happily cope with comments after the block and all the usual perturbations
  such as POSITIONS_FRAC, positionsFrac, positions frac, etc.
  
  """
  
  # Our working string is wstr. First remove any comments
  wstr = line.split("!")[0]
  
  # Remove the actual block 
  wstr = wstr.split()[1:]
  
  # Combine everything that remains, then split by _ and rejoin.
  wstr = "".join(wstr)
  wstr = wstr.split("_")
  wstr = "".join(wstr)
  
  # Convert to lowercase and return
  return wstr.lower()

def parse_cell(lines):
  """ positions, species, properties = parse_cell(lines)
  
  Reads a .cell file and returns the positions, the species and a properties
  dictionary that contains whether positions_frac or positions_abs was used,
  the extra functional lines from the file and anything else useful.
  
  Note we remove all comments.
  
  """
  
  # Strip all the blank lines and comments.
  lines = remove_comments(lines, "!")
  
  # We don't use the for loop here because we may want to sub-loop over blocks.
  i = 0
  species = []
  positions = []
  extralines = []
  properties = {}
  while i < len(lines):
    cline = lines[i]
    if cline.split()[0].lower() == "%block":
      btitle = identify_block(cline)
      if DEBUG:
        print "Found a block", btitle
      if btitle == "positionsfrac" or btitle == "positionsabs":
        if btitle == "positionsfrac":
          properties["fractional"] = True
        else:
          properties["fractional"] = False
        i += 1
        cline = lines[i]
        # Might be a units line
        if len(cline.split()) == 1:
          properties["has units"] = True
          properties["units"] = cline.split()[0].strip()
          i += 1
          cline = lines[i]
        else:
          properties["has units"] = False
        while cline.split()[0].lower() != "%endblock":
          species.append(getElementZ(cline.split()[0]))
          positions.append([float(x) for x in cline.split()[1:4]])
          i += 1
          cline = lines[i]
      else:
        # We aren't interested in the other blocks: just add them to the extra lines.
        extralines.append(cline)
    else:
      # If we aren't reading a block and it isn't a line we're interested in, add
      # to the extra lines.
      extralines.append(cline)
    i += 1
    
  properties["extra lines"] = extralines
  
  return positions, species, properties

def write_new_cell(seed, pos, spec, props, special_index):
  """ write_new_cell(seed, pos, spec, props, special_index)
  
  Writes the new .cell file to seed+"_special.cell" and with the special atom listed
  as "Z:exc" or whatever the custom SPECIAL text is set at the top of the file.
  
  """
  
  f = open(seed+"_special.cell", 'w')
  
  # Write all the extra lines first
  for line in props["extra lines"]:
    f.write(line+"\n")
    if "%endblock" in line.lower():
      f.write("\n")
  
  # Now write the positions block.
  f.write("\n")
  if props["fractional"] == True:
    f.write("%block positions_frac\n")
  else:
    f.write("%block positions_abs\n")
    if props["has units"]:
      f.write(props["units"]+"\n")
  for i in range(len(pos)):
    if i == special_index:
      pstr = "%s%s %g %g %g\n" % (elements[spec[i]], SPECIAL, pos[i][0], pos[i][1], pos[i][2])
    else:
      pstr = "%s %g %g %g\n" % (elements[spec[i]], pos[i][0], pos[i][1], pos[i][2])
    
    f.write(pstr)
  if props["fractional"] == True:
    f.write("%endblock positions_frac\n")
  else:
    f.write("%endblock positions_abs\n")
  
  f.close()
  
  
def index_of_species_index(species, Z, n):
  """ i = index_of_species_index(Z, n)
  
  Returns the absolute index in the species list of the nth species of element Z.
  
  """
  
  si = -1
  for i, s in enumerate(species):
    if DEBUG:
      print "Testing atom %d, element %d to see if it matches %d." % (i, s, n)
    if s == Z:
      si += 1
      if si == n:
        return i
  
  # Didn't find it
  return -1

# Main program
if __name__ == '__main__':

  print "castep_special_atom version 2.0"
  print ""
  print "Written by Kane O'Donnell, September 2013"
  print ""
  
  # Parse the command line arguments
  parser = argparse.ArgumentParser(description="Make a named atom special in a .cell file")
  parser.add_argument('SEED', help="Input file is SEED.cell")
  parser.add_argument('Z', help="Element (e.g. C or 6) to make special.")
  parser.add_argument('n', type=int, help="Make the nth atom of element Z special. 1-based.")
  args = parser.parse_args()
  
  # Convert from 1-based n and possibly-string Z to 0-based n and integer Z.
  n = args.n - 1
  z = getElementZ(args.Z)
  
  # Sanity check on inputs
  if n < 1:
    print "ERROR: It looks like you've specified a negative number for the atomic index. Try again. Exiting..."
    exit(0)
  if z == 0:
    print "ERROR: You've specified an unknown element - it has to be one of the ones in our universe, not imaginary ones! Try again. Exiting..."
    exit(0)
  
  # Read the cell and check for any species "0" - i.e. if we already have special atoms.
  f = open(args.SEED+".cell", 'r')
  lines = f.readlines()
  f.close()
  p, s, props = parse_cell(lines)
  
  if 0 in s:
    print "ERROR: There are unknown species in this file already - you already have at least one special atom. For safety reasons, cannot continue. Goodbye!"
    exit(0)
  
  si = index_of_species_index(s, z, n)
  
  if si == -1:
    print "ERROR: Didn't find atom %d of species %s. Exiting..." % (n, args.Z)
    exit(0)
  
  write_new_cell(args.SEED, p, s, props, si)
  
  print "Goodbye!"
  exit(0)