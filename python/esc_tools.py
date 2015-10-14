################################################################################
#
# esc_tools.py
#
# Library of electronic-structure related code - updated version of esc_lib.
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
# 1. Internally we *always* use angstrom/eV units (unlike the original esc_lib).
#
# 2. Biggest difference here is that esc_tools is geared towards single-structure
# files (no animations). Animations will be handled as a special case when necessary.
#
################################################################################


from __future__ import division
from numpy import *
from numpy.linalg import *
import os

# Debugging flag - set to 1 to see debug messages.
DEBUG=0

# A small number, for float equality comparison
SMALL = 1.0e-6

# Periodic table

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
            116 : "Uuh", 117 : "Uus", 118 : "Uuo", 0 : "UKN" }
            

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
            print "(libesc.getElementZ) Warning: Element %s is not in the elements dictionary. Returning 0 for element %s." % (elstr, elstr)
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

def abc2cell(lengths, angles):
  """ Convert a lattice from a,b,c,alpha,beta,gamma to three cartesian vectors."""
  
  lvec = []
  angles *= pi/180.0
  lvec.append(array([lengths[0], 0.0, 0.0])) # "a"
  lvec.append(array([lengths[1] * cos(angles[2]), lengths[1] * sin(angles[2]), 0.0]))
  ca = lengths[2] * cos(angles[1])
  cb = lengths[2] * (cos(angles[0]) - cos(angles[1]) * cos(angles[2])) / sin(angles[2])
  cc = sqrt(lengths[2] ** 2 - ca ** 2 - cb ** 2)
  lvec.append(array([ca, cb, cc]))
  
  return lvec
  
def cell2abc(lvec):
  """ Convert cartesian vectors to length/angle format. 
  
  Note that there is ambiguity in the length/angle format and you might need to re-express
  any atomic positions to avoid them being shifted outside the cell.
  
  """
  
  a = norm(lvec[0])
  b = norm(lvec[1])
  c = norm(lvec[2])
  gamma = arccos(dot(lvec[0], lvec[1])/(a * b)) * 180 / pi
  alpha = arccos(dot(lvec[1], lvec[2])/(b * c)) * 180 / pi
  beta = arccos(dot(lvec[2], lvec[0])/(c * a)) * 180 / pi
  
  return a, b, c, alpha, beta, gamma
  
def cart2reduced(position, lattice):
    """ reduced = cart2reduced(position, lattice)
    
    Converts a cartesian coordinate to a reduced coordinate with respect to
    the lattice. Works recursively on lists.
    
    """
    
    if type(position) == type([]):
      return [cart2reduced(x, lattice) for x in position]
    else:
      return array(mat(position) * mat(lattice).I).reshape((3,))
      
def reduced2cart(position, lattice):
    """ cart = reduced2cart(position, lattice)
    
    Converts a reduced coordinate to a cartesian coordinate with respect to the
    supplied lattice. Works recursively on lists.
    
    """
    
    if type(position) == type([]):
      return [reduced2cart(x, lattice) for x in position]
    else:
      return array(mat(position) * mat(lattice)).reshape((3,))
            
def bohr2ang(bohr):
    """ ang = bohr2ang(bohr)
    
    Converts bohr units to angstrom, conversion factor 1 bohr = 0.52917721092 ang
    
    Woohoo recursive function!
    
    """
    
    if type(bohr) == type([]):
        # Call on each element and return
        return [bohr2ang(x) for x in bohr]
    else:
        return bohr * 0.52917721092
    
def ang2bohr(ang):
    """ bohr = ang2bohr(ang)
    
    Converts angstroms to bohr, conversion factor 1 ang = 1.0 / 0.52917721092 bohr
    
    """
    
    if type(ang) == type([]):
        return [ang2bohr(x) for x in ang]
    else:
        return ang / 0.52917721092

def eV2hartree(eV):
    """ hartree = eV2hartree(eV)
    
    Converts eV to hartree, conversion factor 1 Ht = 27.21138505 eV
    
    """
    
    if type(eV) == type([]):
        return [eV2hartree(x) for x in eV]
    else:
        return eV / 27.21138505
    
def hartree2eV(hartree):
    """ eV = hartree2eV(hartree)
    
    Converts hartree to eV, conversion factor 1 eV = 1.0 / 27.21138505 Ht
    
    """
    
    if type(hartree) == type([]):
        return [hartree2eV(x) for x in hartree]
    else:
        return hartree * 27.21138505
        
def eV2rydberg(eV):
    """ rydberg = eV2rydberg(eV)
    
    Converts eV to rydberg, conversion factor 1 Ht = 13.6 eV
    
    """
    
    if type(eV) == type([]):
        return [eV2rydberg(x) for x in eV]
    else:
        return eV / 13.605698066
    
def rydberg2eV(rydberg):
    """ eV = rydberg2eV(rydberg)
    
    Converts rydberg to eV, conversion factor 1 eV = 1.0 / 13.6 Ry
    
    """
    
    if type(rydberg) == type([]):
        return [rydberg2eV(x) for x in rydberg]
    else:
        return rydberg * 13.605698066
        
def theta_z(p, q):
  """ Returns the angle between the z axis (0,0,1) and the vector
      formed between the two passed points p and q. Note the angle
      is always the acute angle - this means the returned angle
      is independent of the order of p and q. """
      
  v = p - q
  theta = arccos(v[2] / norm(v)) * 180.0 / pi
  if theta > 90.0:
    theta -= 90.0
  
  return theta
  
def plane_z(p, q, r):
  """ Returns the angle between the plane formed by the points p, q
  and r, and the z axis (0,0,1). Checks for colinearity of the three
  points and returns an error if they are collinear."""
  
  # Collinearity check: form two vectors from the three points,
  # check scalar product is equal to abs dot product (e.g. parallel).
  pq = p - q
  pr = p - r
  if abs(dot(pq, pr)) - norm(pq) * norm(pr) > SMALL:
    print "Error (plane_z): Points p, q and r are collinear! Cannot form plane."
    return -1.0
  
  # Generate the unit normal to the plane formed by p, q and r
  n = cross(pq, pr)
  n /= norm(n)
  
  theta = arccos(n[2]) * 180.0 / pi
  if theta > 90.0:
    theta -= 90.0
  
  return theta
        
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
    
def indexLine(text, substring, returnAll=False):
  """ index = indexLine(text, substring, returnAll=False)
  
  Returns the index of the line where a substring first occurs
  in a list of strings. If returnAll=True, returns a list of indices
  of all occurrences. 
  
  Note we return None if the string is not found at all.
  
  This is a function that combines the properties of substringInList and 
  substringPositionsInList (the other two are deprecated).
  
  """
    
  indices = [i for i,line in enumerate(text) if substring in line]
  if returnAll:
    if len(indices) == 0:
      return None
    else:
      return indices
  else:
    # The substring might not occur at all, in which case we return None.
    if len(indices) == 0:
      return None
    else:
      return indices[0]
                    
class Atoms:
  """ Atoms class: a collection of atoms (possibly with a crystal structure)
    
  Create using Atoms(xsf_file), Atoms(xyz_file), Atoms(abinit_input), etc.
    
  """
    
  nsteps = 1                  # > 1 if object has time-series data
                                # (ie is animated)
  lattice = None                # array of 3 vectors
  positions = None              # array of vectors
  forces = None                 # Same for forces...
  species = None                # list of atomic numbers
  parameters = {}             # dictionary. Contains non-structural information.
  parameters_type = ""
  animation_species = []
  animation_positions = []
  kpoints_grid = [1,1,1]
  kpoints_offset = [0,0,0]
  

    
  def clean(self):
    """ Atoms.clean()
    
    Internal - wipes out data in the instance object, DON'T USE!
    
    """
      
    self.lattice = None
    self.positions = None
    self.forces = None
    self.species = None
    self.parameters = {}
    self.parameters_type = ""
    self.animated_species = []
    self.animated_positions = []
    self.kpoints_grid = [1,1,1]
    self.kpoints_offset = [0,0,0]
    
  def __init__(self, filename, filetype="pdb", options=None):
      """ atoms = Atoms(filename, filetype="pdb")
      
      Creates an Atoms object from some input file. You have to
      specify the type of file here. Can be:
      
      "XSF" : XCrysden Structure Format, can also be animated axsf.
      
      "abinit,input" : Abinit input file. !!! NOT CODED
      
      "elk,input" : Elk input file. !!! NOT CODED
      
      "castep,cell" : CASTEP .cell file.
      
      "aims,geometry" : geometry.in file from FHI-aims
      
      "aims,output" : Output file from FHI-aims. !!! NOT CODED
      
      "pdb" : Protein Data Bank format.
        
      "qe,input" : Quantum Espresso input file.
        
      """
        
      if filetype == "XSF":
        self.loadFromXSF(filename)
      #elif filetype == "abinit,input":
        #    self.loadFromAbinit(filename)
      elif filetype == "elk,input":
        self.loadFromElk(filename)
      elif filetype == "castep,cell":
        self.loadFromCASTEP(filename)
      elif filetype == "aims,geometry":
        self.loadFromAimsGeometry(filename)
      #elif filetype == "aims,output":
      #  self.loadFromAimsOutput(filename)
      elif filetype == "pdb":
        self.loadFromPDB(filename)
      elif filetype == "qe,input":
        self.loadFromQEInput(filename)
      elif filetype == "test":
        self.loadFromQEInput(filename)
      else:
        print "(esc_lib.Atoms.__init__) ERROR: File type %s not handled at present." % filetype
        return None
          
  def loadFromQEInput(self, filename):
    """ Loads a Quantum Espresso input file into Atoms object.
    
    Note: we always reject input files that aren't specified by CELL_PARAMETERS - unlike QE itself,
    CELL_PARAMETERS is the mandatory cell method for this code.
    
    """
    
    self.clean()
    f = open(filename, 'r')
    data = f.readlines()
    data = remove_comments(data, "!")
    data = remove_comments(data, "#")
    f.close()
    
    self.parameters_type = "QE"
        
    if DEBUG:
      print "Removed comments from file."
    
    # Need this to read the &namelists because they have a different format to the cards.
    def process_namelist(start):
      line = data[start].strip().lower()
      if DEBUG:
        print "process_namelist on for line '%s'" % (line)
      namelist_name = line.split()[0]
      namelist_lines = []
      j = start + 1
      while data[j].strip().lower()[0] != '/':
        namelist_lines.append(data[j].strip())
        j += 1
      params = {}
      namelist_contents = "".join(namelist_lines).split(',')
      namelist_contents = [x for x in namelist_contents if x != '']
      if DEBUG:
        print "namelist_contents is '%s'" % (namelist_contents)
      for param in namelist_contents:
        param_name = param.split("=")[0].strip().lower()
        param_value = param.split("=")[1].strip()
        params[param_name] = param_value
      self.parameters[namelist_name] = params
      return
    
    # Check all mandatory keywords are present.
    mandatory_keywords = ["&control", "&system", "&electrons", "atomic_species", "atomic_positions",
      "k_points", "cell_parameters"]
    
    dataline = "".join(data).lower()
    for m in mandatory_keywords:
      if m not in dataline:
        print "Error: Mandatory block %s missing from QE input file. Exiting..." % m
        return None
    
    # Now locate all blocks. Only the official blocks from the pw.x input are recognized.
    blocks = {}
    keywords = mandatory_keywords + ["&ions", "&cell", "constraints", "occupations", "atomic_forces"]
    
    for i, line in enumerate(data):
      if line.split()[0].strip().lower() in keywords:
        # Need a check here because there's an option "occupations" and a block "occupations"
        # and both will be found by this search.
        if "=" not in line:
          blocks[line.split()[0].strip().lower()] = i
        
    # Need to read the &system namelist first, then can read all the others.
    process_namelist(blocks["&system"])
    
    for k in blocks.keys():
      if k[0] == "&" and k != "&system":
        process_namelist(blocks[k])
      elif k == "cell_parameters":
        i = blocks[k]
        units = data[i].split()[1].lower()
        v1 = array([float(x) for x in data[i+1].strip().split()])
        v2 = array([float(x) for x in data[i+2].strip().split()])
        v3 = array([float(x) for x in data[i+3].strip().split()])
        self.lattice = [v1, v2, v3]
        # Unit conversions
        if units == "bohr":
          self.lattice = bohr2ang(self.lattice)
      elif k == "atomic_positions":
        i = blocks[k]
        units = data[i].strip().split()[1].lower()
        pos = []
        spec = []
        for j in range(int(self.parameters["&system"]["nat"])):
          bits = data[i+j+1].strip().upper().split()
          pos.append(array([float(x) for x in bits[1:]]))
          spec.append(getElementZ(bits[0]))
        self.positions = pos
        self.species = spec
        if units == "bohr":
          self.positions = bohr2ang(self.positions)
        elif units == "alat" or units == "crystal":
          self.positions = reduced2cart(self.positions, self.lattice)
      elif k == "atomic_species":
        i = blocks[k]
        tmp = []
        for j in range(int(self.parameters["&system"]["ntyp"])):
          tmp.append(data[i+j+1].strip())
        self.parameters["atomic_species"] = tmp
      elif k == "k_points":
        i = blocks[k]
        kpts_type = data[i].strip().split()[1].lower()
        # Can only deal with automatic kpoints at the moment.
        if kpts_type == "automatic":
          kpts_line = data[i+1].strip().split()
          self.kpoints_grid = array([int(x) for x in kpts_line[0:3]])
          self.kpoints_offset = array([int(x) for x in kpts_line[3:]])
      elif k == "constraints":
        i = blocks[k]
        nconstr = int(data[i+1].split()[0])
        # Don't parse this card, just put all the lines in a list.
        tmp = [data[i+1]]
        for j in range(nconstr):
          tmp.append(data[i+j+1].strip())
        self.parameters["constraints"] = tmp
      elif k == "occupations":
        # Man I hate this card. Basically, we test each line. If each element
        # of each line after "occupations" can be converted to a float, we take
        # that line, otherwise we end.
        i = blocks[k] + 1
        tmp = []
        while True:
          bits = data[i].split()
          try:
            tmpf = [float(x) for x in bits]
          except ValueError:
            break
          tmp.append(data[i].strip)
          i += 1
        self.parameters["occupations"] = tmp
      elif k == "atomic_forces":
        i = blocks[k]
        tmp = []
        for j in range(int(self.parameters["&system"]["nat"])):
          tmp.append(data[i+j+1].strip())
        self.parameters["atomic_force"] = tmp


  def writeQEInput(self, filename, xtype="ang", opt={}):
    """ Writes a QE input file based on the Atoms object. If the parameters_type
        is set to QE, we use the parameters dictionary to fill out the input file,
        otherwise we list the mandatory blocks and parameters with default settings. """
        
    
    print xtype
    f = open(filename, 'w')
    if not f:
      print "Error: could not write to file %s." % filename
      return None
      
    if xtype == "bohr":
      lv = ang2bohr(self.lattice)
      pos = ang2bohr(self.positions)
    elif xtype == "ang":
      lv = self.lattice
      pos = self.positions
    elif xtype == "alat":
      lv = self.lattice
      pos = cart2reduced(self.positions, self.lattice)
      
    if "corehole" in opt.keys():
      add_corehole = True
      corehole_index = opt["corehole"]
    else:
      add_corehole = False
      
    if "suffix atoms" in opt.keys():
      suffix_some_atoms = True
      suffix = "x"
      atoms_to_suffix = opt["suffix atoms"]
    else:
      suffix_some_atoms = False
    
    if "charge" in opt.keys():
      # Note that charged is used as a template for a core-hole cell, so we add an atom type as well.
      charged_cell = True
      charge = opt["charge"]
    else:
      charged_cell = False
      
    if "prefix" in opt.keys():
      prefix = opt["prefix"]
    else:
      prefux = "pwscf"
    
    if "bands" in opt.keys():
      bands = "nbnd = " + str(opt["bands"])
    else:
      bands = "!nbnd = BANDS"
    
    # This is a bit tricky. If we have QE parameters, write all the namelists first,
    # then write any cards (except atomic_species) AFTER the important cards below.
    def write_namelist(key):
      f.write("%s\n" % key)
      for k in self.parameters[key].keys():
        f.write("  %s = %s,\n" % (k, self.parameters[key][k]))
      f.write("/\n")
      
    if self.parameters_type == "QE":
      write_namelist("&control")
      write_namelist("&system")
      write_namelist("&electrons")
      
      if "&ions" in self.parameters.keys():
        write_namelist("&ions")
      if "&cell" in self.parameters.keys():
        write_namelist("&cell")

    else:
      f.write("&control\n")
      f.write("  calculation = 'scf',\n  title = '',\n  outdir = './',\n  prefix = '%s',\n  pseudo_dir = './',\n  wf_collect = .true.\n/\n" % (prefix))
      f.write("&system\n")
      if charged_cell:
        f.write("  ibrav = 0,\n  nat = %d,\n  ntyp = %d,\n  %s,\n  ecutwfc = 65,\n  ecutrho = 280,\n  tot_charge=%g,\n  occupations = 'smearing',\n  smearing = 'mv',\n  degauss = 0.0073,\n  nspin = 1\n/\n" % (len(self.species), len(uniqify(self.species))+1, bands, charge))
      else:
        f.write("  ibrav = 0,\n  nat = %d,\n  ntyp = %d,\n  %s,\n  ecutwfc = 65,\n  !ecutrho = 280,\n  !tot_charge=+1.0,\n  !occupations = 'fixed',\n  smearing = 'mv',\n  !degauss = 0.02,\n  nspin = 1\n/\n" % (len(self.species), len(uniqify(self.species)), bands))
      f.write("&electrons\n")
      f.write("  conv_thr = 1.0D-7\n/\n")
      f.write("&ions\n/\n")
    
    f.write("ATOMIC_SPECIES\n")
    if charged_cell:
      f.write("  Xh 1.0 PSEUDO_GOES_HERE\n")
    if "atomic_species" in self.parameters.keys():
      for b in self.parameters["atomic_species"]:
        f.write("  %s\n" % (b))
    else:
      for s in uniqify(self.species):
        f.write("  %s 1.0 PSEUDO_GOES_HERE\n" % (elements[s]))
    
    f.write("CELL_PARAMETERS")
    if xtype == "bohr":
      f.write(" bohr\n")
    elif (xtype == "ang") or (xtype == "alat"):
      f.write(" angstrom\n")
    for v in lv:
      f.write("  %g %g %g\n" % (v[0], v[1], v[2]))
    f.write("K_POINTS automatic\n")
    if self.kpoints_grid is not None:
      f.write("  %d %d %d %d %d %d\n" % (self.kpoints_grid[0], self.kpoints_grid[1], self.kpoints_grid[2], self.kpoints_offset[0], self.kpoints_offset[1], self.kpoints_offset[2]))
    else:
      f.write("  1 1 1 0 0 0\n")
    if (xtype == "alat") or (xtype == "bohr"):
      f.write("ATOMIC_POSITIONS %s\n" % xtype)
    elif xtype == "ang":
      f.write("ATOMIC_POSITIONS angstrom\n")
    for i in range(len(pos)):
      v = pos[i]
      s = elements[self.species[i]]
      if add_corehole and (i == corehole_index):
        s += "h"
      if suffix_some_atoms and i in atoms_to_suffix:
        s += suffix
      f.write("%s    %g    %g    %g\n" % (s, v[0], v[1], v[2]))
    
    # Ok now write any cards from self.parameters (except atomic_species).
    def write_card(key):
      f.write("%s\n" % key.upper())
      for b in self.parameters[key]:
        f.write("  %s\n" % b)
        
    if "constraints" in self.parameters.keys():
      write_card("constraints")
    if "occupations" in self.parameters.keys():
      write_card("occupations")
    if "atomic_forces" in self.parameters.keys():
      write_card("atomic_forces")
    f.close()
    return True
    
  def loadFromPDB(self, filename):
    """ atoms = Atoms.loadFromPDB(filename)
    
    Internal: grabs a molecule or crystal from a .pdb file. Not very clever but good
    enough. Anything that comes up as ERROR will return None.
  
    """
  
    self.clean()
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
  
    # We are only interested in the following lines: anything that starts with "ATOM",
    # "CRYST1", "HETATM", "NUMMDL", "MODEL" and "ENDMDL".
    crys = indexLine(data, "CRYST1")
    atoms = indexLine(data, "ATOM", returnAll=True)
    hetatms = indexLine(data, "HETATM", returnAll=True)
    nummdl_idx = indexLine(data, "NUMMDL")
    models = indexLine(data, "MODEL", returnAll=True)
    endmdls = indexLine(data, "ENDMDL", returnAll=True) 
    num_models = None # will set this shortly if necessary
  
    # Basic sanity checks.
    # 1. Check we actually have atoms.
    if atoms is None and hetatms is None:
      print "(Atoms.loadFromPDB) ERROR: PDB file contains no ATOM or HETATM statements - there are no atoms specified!"
      return None
    # 2. If NUMMDL is greater than 1 - need to make sure there
    # are the correct number of blocks and ATOM/HETATM statements.
    if nummdl_idx is not None:
      num_models = int(data[nummdl_idx].split()[1])
      if num_models > 1:
        if len(models) != len(endmdls):
          print "(Atoms.loadFromPDB) ERROR: PDB file is incomplete - number of MODEL and ENDMDL statements do not match."
          return None
        elif len(models) != num_models:
          print "(Atoms.loadFromPDB) ERROR: PDB file is incomplete - number of MODEL blocks is not equal to NUMMDL."
          return None
        # Check that the number of ATOM/HETATM lines is a multiple of NUMMDL.
        combined_atoms = len(atoms) + len(hetatms)
        if combined_atoms % num_models != 0:
          print "(Atoms.loadFromPDB) ERROR: Number of atoms is not a multiple of NUMMDL - there are missing atoms in the model sequence."
          return None
        self.nsteps = num_models
        
    # If we have a crystal lattice specified, convert it to vectors. Note that since
    # in a PDB the crystal is specified in ABC format (three lengths, three angles)
    # the conversion is not unique and needs to be done according to a standard. Here,
    # the CASTEP standard is used: "a" goes along x, b is in the xy plane and the whole
    # thing must have a positive cell volume with a.(bxc) [this fixes c].
    if crys is not None:
      lvec = []
      lengths = array([float(x) for x in data[crys].split()[1:4]])
      angles = array([float(x) for x in data[crys].split()[4:7]])
      lvec = abc2cell(lengths, angles)
    
      self.lattice = lvec
  
    # Read every ATOM/HETATM for positions and species. Note that we combine the two
    # index lists and then sort, so that it doesn't matter whether there is exclusively
    # ATOM, exclusively HETATM or mixed in any combination. The result is then chopped
    # by the number of models in order to get the frames of animation.
    if atoms is not None and hetatms is not None:
      all_atoms = atoms + hetatms
    elif hetatms is not None:
      all_atoms = hetatms
    elif atoms is not None:
      all_atoms = atoms
    # The else case is dealt with earlier.
  
    all_atoms.sort()
    pos = []
    spec = []
    for i in all_atoms:
      at = data[i]
      pos.append(array([float(x) for x in at[30:54].split()]))
      spec.append(getElementZ(at[76:78].strip()))
  
    # The way we deal with animated data in esc_tools is to put all the positions
    # in self.animated_positions and self.animated_species and the last position
    # in self.positions and self.species.
    if num_models is not None:
      chunksize = int(float(len(pos)) / num_models)
      for j in range(num_models):
        self.animated_positions.append(pos[j*chunksize:(j+1)*chunksize])
        self.animated_species.append(spec[j*chunksize:(j+1)*chunksize])
      self.positions = self.animated_positions[-1]
      self.species = self.animated_species[-1]
    else:
      self.positions = pos
      self.species = spec
      
  def writePDB(self, filename, opt={}):
    """ Write a PDB file.
    
    Some notes about compliance:
    
    1. Uses "ATOM" keyword for all atoms (never HETATM).
    2. If a lattice is passed, this is written in ABC format to CRYST1 but the space
    group is not calculated and is therefore void (by default P 1).
    3. Does not include any bonding information.
    4. There is a rotational ambiguity going from cartesian to abc lattice format. Due to this, positions have to be re-calculated here because the PDB format only accepts absolute positions not fractional relative to the abc basis.
    
    """
  
    pos = self.positions
    spec = self.species
    
    if "corehole" in opt.keys():
      add_corehole = True
      corehole_index = opt["corehole"]
    else:
      add_corehole = False
      
    if "suffix atoms" in opt.keys():
      suffix_some_atoms = True
      suffix = "x"
      atoms_to_suffix = opt["suffix atoms"]
    else:
      suffix_some_atoms = False
  
    if self.lattice is not None and self.lattice != []:
      avec = self.lattice
      is_crystal = True
    else:
      is_crystal = False
    
    f = open(filename, 'w')
    f.write("REMARK written by esc_lib.py (Kane O'Donnell, kaneod on GitHub)\n")
  
    if is_crystal:
      # Get the lengths of the vectors and the angles between them.
      a, b, c, alpha, beta, gamma = cell2abc(avec)
    
      # Convert these back to a new set of vectors so we can re-express the positions.
      lvec = abc2cell(array([a,b,c]), array([alpha, beta, gamma]))
    
      # First convert to fractional in the old basis
      pos = cart2reduced(pos, avec)
      # Now convert to cartesian in the new basis.
      pos = reduced2cart(pos, lvec)
      
      if DEBUG:
        print "Found lengths and angles:"
        print "%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f" % (a, b, c, alpha, beta, gamma)
    
      f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1            1\n" % (a, b, c, alpha, beta, gamma))
    
    for i in range(len(pos)):
      p = pos[i]
      s = elements[spec[i]]
      if add_corehole and (i == corehole_index):
        s += "h"
      if suffix_some_atoms and i in atoms_to_suffix:
        s += suffix      
      specname = s.rjust(2)
      f.write("ATOM  %5d %s     X     1    %8.3f%8.3f%8.3f                      %s  \n" % \
                    (i+1, specname.upper(), p[0],p[1],p[2],specname))
  
    f.write("END\n")
    f.close()
    return True

  def loadFromAimsOutput(self, filename):
    """ atoms = Atoms.loadFromAimsOutput(filename)
    
    Internal: loads geometry from a FHI-aims output file. The idea here is that
    we can get the entire geometry sequence, for example, or phonon displacements or
    whatever.
    
    """
    
    self.clean()
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    # Reading the output is harder than the geometry: we have to look for markers
    # in the output to tell us where we are in the file. Also, because the user
    # can cancel the verbatim writing of the control.in and geometry.in files, we can't
    # rely on those to start with: our first marker must be one that is always there.
    # So, the initial geometry read is DIFFERENT from the subsequent ones.
    
    pos = []
    spec = []
    lvec = []
    
    # First thing we look for is "Reading geometry description geometry.in.".
    idx = substringInList("Reading geometry description geometry.in.", data)
    
    # Next, we read how many atoms there are.
    numatidx = substringInList("| Number of atoms", data)
    numat = int(data[numatidx].split()[5])
    
    if idx:
      # Chop data, locate if we have unit cell or not.
      data = data[idx:]
      idx = substringInList("No unit cell requested.", data)
      if idx:
        # Molecule, not crystal. Get the initial positions and species.
        for line in data[(idx+3):(idx+3+numat)]:
          spec.append(getElementZ(line.split()[3]))
          pos.append(array([float(x) for x in line.split()[4:7]]))
        self.positions = pos
        self.species = spec
        self.is_crystal = False
      else:
        # Molecule is a crystal. Get the lattice vectors then the initial pos and spec.
        idx = substringInList("| Unit cell:", list)
        for line in data[(idx+1):(idx+4)]:
          lvec.append(array([float(x) for x in line.split()[1:4]]))
        for line in data[(idx+6):(idx+6+numat)]:
          spec.append(getElementZ(line.split()[3]))
          pos.append(array([float(x) for x in line.split()[4:7]]))
        self.positions = pos
        self.species = spec         
        self.lattice = lvec
        self.is_crystal = True
    else:
      # This out file is broken: Didn't even get to initial calculation.
      print "(Atoms.loadFromAimsOutput) ERROR: File is not a complete FHI-aims output. Returning None."
      return None
    
    # Now get all the coordinate updates.
    idxs = substringPositionsInList("Updated atomic structure:", data)
    for i in idxs:
      pos = [self.positions]
      for line in data[(i+2):(i+2+numat)]:
        pos.append(array([float(x) for x in line.split()[1:4]]))
      self.animated_positions.append(array(pos))
      self.animated_species.append(self.species)
      self.positions = self.animated_positions[-1]
               
  def loadFromAimsGeometry(self, filename):
    """ atoms = Atoms.loadFromAimsGeometry(filename)
    
    Internal, inits an Atoms object from a FHI-aims geometry.in file.
    
    """
    
    self.clean()
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    data = remove_comments(data, "#")
    # aims is case-independent so we are too.
    data = [s.lower() for s in data]
    
    # Easy: just have to read through each line and look for "atom" lines
    # or "lattice_vector" lines or "atom_frac" lines.
    pos = []
    spec = []
    lvec = []
    
    for line in data:
      if line.split()[0] == "atom":
        pos.append(array([float(x) for x in line.split()[1:4]]))
        spec.append(getElementZ(line.split()[4]))
        is_atoms = True # sets a flag to help us later.
      elif line.split()[0] == "lattice_vector":
        lvec.append(array([float(x) for x in line.split()[1:4]]))
      elif line.split()[0] == "atom_frac":
        pos.append(array([float(x) for x in line.split()[1:4]]))
        spec.append(getElementZ(line.split()[4]))
        is_atoms = False
    
    # Ok bit of a bind here. For now, let's assume aims requires either
    # three lattice vectors or none. If there is less than three, we assume
    # the coords are actually abs positions. If there's 3 lattice vectors,
    # we use the flag set up in the previous loop to decide whether we have
    # atoms for atom_fracs (assume all the same).
    
    if len(lvec) < 3:
      self.positions = pos
      self.species = spec
      self.is_crystal = False
    else:
      self.lattice = lvec
      if is_atoms:
        self.positions = pos
      else:
        self.positions = reduced2cart(pos, lvec)
      self.species = spec
      self.is_crystal = True
    
  def loadFromCASTEP(self, filename):
    """ atoms = Atoms.loadFromCASTEP(filename)
    
    Internal, inits an Atoms object from a CASTEP .cell file. Note that
    we don't handle a lot of units here: specifically, we check for the
    presence of an explicit units line but only convert for angstrom or bohr.
    
    We try to handle as much of the input file flexibility as possible:
    statements such as positions_frac and positionsfrac and PositionsFRAC and
    so on are all recognized, and we remove comments and blanks before the
    parse so they don't get in the way. We also allow the assignment variants
    kpoint_mp_grid=1 1 1 or kpoint_mp_grid: 1 1 1 or kpoint_mp_grid 1 1 1.
    
    TODO: The lattice_abc block.
    
    """
    
    self.clean()
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    # Strip away all comments and end of line comments.
    data = remove_comments(data, "!")
    
    # CASTEP is case-independent so we should be too.
    data = [s.lower() for s in data]
    
    if DEBUG:
      for i, line in enumerate(data):
        print i, line 
        
    self.parameters = {}
    self.parameters['CASTEP'] = {}
    
    # Line-by-line parse
    i = 0
    postype = "None"
    while i < len(data):
      line = data[i]
      if line.split()[0] == "%block":
        btitle =  line.split()[1]
        btitle = "".join(btitle.split("_"))
        if DEBUG:
          print "Found block: ", btitle
        if btitle == "latticecart":
          # Either three or four lines in this block.
          if data[i+4].split()[0] == "%endblock":
            # No units, assume angstrom.
            vec1 = array([float(x) for x in data[i+1].split()])
            vec2 = array([float(x) for x in data[i+2].split()])
            vec3 = array([float(x) for x in data[i+3].split()])
            self.lattice = [vec1, vec2, vec3]
            i = i + 4
          elif data[i+5].split()[0] == "%endblock":
            units = data[i+1].split()[0].lower()
            vec1 = array([float(x) for x in data[i+2].split()])
            vec2 = array([float(x) for x in data[i+3].split()])
            vec3 = array([float(x) for x in data[i+4].split()])
            if units == "ang" or units == "angstrom":
              self.lattice = [vec1, vec2, vec3]
            elif units == "bohr":
              self.lattice = bohr2ang([vec1, vec2, vec3])
            i = i + 5
        elif btitle == "positionsabs":
          # Loop to the end of this block
          postype = "absolute"
          pos = []
          specs = []
          unit = "ang"
          for j in range(i+1,len(data)):
            if data[j].split()[0] == "%endblock":     
              i = j
              break
            elif len(data[j].split()) == 1:
              unit = data[j].split()[0].lower()
            else:
              specs.append(getElementZ(data[j].split()[0]))
              pos.append(array([float(s) for s in data[j].split()[1:4]]))
          if unit == "ang" or unit == "angstrom":
            self.positions = pos
          elif unit == "bohr":
            self.positions = bohr2ang(pos)
          self.species = specs
        elif btitle == "positionsfrac":
          # Loop to the end of this block
          postype = "fractional"
          pos = []
          specs = []
          for j in range(i+1,len(data)):
            if data[j].split()[0] == "%endblock":            
              i = j
              break
            else:
              specs.append(getElementZ(data[j].split()[0]))
              pos.append(array([float(s) for s in data[j].split()[1:4]]))
          self.species = specs
      else:
        # Line is outside a block
        # Look for "=" or ":"
        if ":" in line:
          option = "".join(line.split(":")[0].split("_"))
          value = line.split(":")[1].strip()
        elif "=" in line:
          option = "".join(line.split("=")[0].split("_"))
          value = line.split("=")[1].strip()
        else:
          option = "".join(line.split()[0].split("_"))
          value = " ".join(line.split()[1:]).strip()
        
        # Add these to the parameters under the CASTEP dictionary.
        self.parameters['CASTEP'][option] = value
        
        # We have a special place for kpoint parameters.
        if option in ["kpointsmpgrid", "kpointmpgrid"]:
          self.kpoints_grid = [int(x) for x in value.split()]
        if option in ["kpointsmpoffset", "kpointmpoffset"]:
          self.kpoints_offset = [int(x) for x in value.split()]
          
        if DEBUG:
          print "Found option: ", option, " with value ", value
      i = i + 1
    
    if postype == "fractional":
      self.positions = reduced2cart(pos, self.lattice)
    
  def writeAims(self, filename="geometry.in", xtype="ang", opt={}):
    """ succeeded = write_aims(filename, positions, species, lattice, xtype="ang", opt=None)

      Writes a FHI-aims geometry.in file. 

      If a lattice is present but xtype is ang we still write the lattice.
    
      opt is a dictionary where the key is the name of the option and the value is some
      additional data. Options implemented so far are:
    
      "constrain atoms": value is a list of atomic indices to constrain.
    
      "constrain species": value is a list of species to constrain.

    """

    pos = self.positions
    spec = self.species
  
    if "constrain atoms" in opt.keys():
      index_constraints = True
    else:
      index_constraints = False
   
    if "constrain species" in opt.keys():
      species_constraints = True
      constraint_list_species = [getElementZ(x) for x in opt["constrain species"]]
    else:
      species_constraints = False 
      
    if "corehole" in opt.keys():
      add_corehole = True
      corehole_index = opt["corehole"]
    else:
      add_corehole = False
      
    if "suffix atoms" in opt.keys():
      suffix_some_atoms = True
      suffix = "x"
      atoms_to_suffix = opt["suffix atoms"]
    else:
      suffix_some_atoms = False

    if xtype == "frac":
      avec = self.lattice
      pos = cart2reduced(pos, avec)
    elif xtype == "ang":
      if self.lattice is not None:
        avec = self.lattice
      else:
        avec = None
      pos = self.positions
    else:
      print "write_aims ERROR: Must specify xtype=ang or frac."
      return False

    f = open(filename, 'w')
    f.write("# geometry.in written by esc_tools.py\n\n")
    if xtype == "ang":  
      if avec is not None:
        for l in avec:
          f.write("lattice_vector %4.8g %4.8g %4.8g\n" % (l[0], l[1], l[2]))
      for i in range(len(pos)):
        p = pos[i]
        s = elements[spec[i]]
        if add_corehole and (i == corehole_index):
          s += "h"
        if suffix_some_atoms and i in atoms_to_suffix:
          s += suffix
        f.write("atom  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], s))
        if add_corehole and (i == corehole_index):
          # Put an initial spin on the atom with the core hole.
          f.write("initial_moment 1.0\n")
        if index_constraints and i in opt["constrain atoms"]:
          f.write("  constrain_relaxation .true.\n")
        elif species_constraints and s in opt["constrain species"]:
          f.write("  constrain_relaxation .true.\n")
    elif xtype == "frac":
      for l in avec:
        f.write("lattice_vector %4.8g %4.8g %4.8g\n" % (l[0], l[1], l[2]))
    
      f.write("\n")
      for i in range(len(pos)):
        p = pos[i]
        s = elements[spec[i]]
        if add_corehole and (i == corehole_index):
          s += "h"
        if suffix_some_atoms and i in atoms_to_suffix:
          s += suffix
        f.write("atom_frac  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], s))
        if add_corehole and (i == corehole_index):
          # Put an initial spin on the atom with the core hole.
          f.write("initial_moment 1.0\n")
        if index_constraints and i in opt["constrain atoms"]:
          f.write("  constrain_relaxation .true.\n")
        elif species_constraints and elements[spec[i]] in constraint_list_species:
          f.write("  constrain_relaxation .true.\n")
    f.close()

    return True
    
  def writeCASTEP(self, filename, xtype="ang", opt={}):
    """ succeeded = writeCASTEP(filename, xtype="ang", opt={})
  
    Writes a CASTEP .cell file using the passed positions. Options for xtype are 
    "ang", "bohr" or "frac".
    We assume, as always, that EVERYTHING internal is in ang or eV, so a 
    conversion is performed if "frac" is specified using the passed lattice
    vectors. Also, if fractional output is specified, the lattice vectors
    are output in angstroms, the CASTEP default length unit.
  
    opt is a dictionary used to pass in custom write modifications.
  
    """
  
    pos = self.positions
    avec = self.lattice
    spec = self.species
  
    # Do conversions if necessary.
    if xtype == "bohr":
      pos = ang2bohr(pos)
      avec = ang2bohr(avec)
    elif xtype == "frac":
      pos = cart2reduced(pos,avec)
      
    if "corehole" in opt.keys():
      add_corehole = True
      corehole_index = opt["corehole"]
    else:
      add_corehole = False
      
    if "suffix atoms" in opt.keys():
      suffix_some_atoms = True
      suffix = "x"
      atoms_to_suffix = opt["suffix atoms"]
    else:
      suffix_some_atoms = False
  
    f = open(filename, 'w')
  
    f.write("%block lattice_cart\n")
    if xtype == "bohr":
      f.write("  bohr\n")
    for v in avec:
      f.write("  %010e %010e %010e\n" % (v[0], v[1], v[2]))
    f.write("%endblock lattice_cart\n")
    f.write("\n")
    if xtype == "frac":
      f.write("%block positions_frac\n")
    else:
      f.write("%block positions_abs\n")
      if xtype == "bohr":
        f.write("  bohr\n")
    for i in range(len(pos)):
        p = pos[i]
        s = elements[spec[i]]
        if add_corehole and (i == corehole_index):
          s += "h"
        if suffix_some_atoms and i in atoms_to_suffix:
          s += suffix        
        f.write("  %s %010e %010e %010e\n" % (s, p[0], p[1], p[2]))
    if xtype == "frac":
      f.write("%endblock positions_frac\n")
    else:
      f.write("%endblock positions_abs\n")
  
    f.close()
    return True
  
  