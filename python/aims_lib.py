################################################################################
#
# aims_lib.py
#
# Library of electronic-structure related code specifically for FHI-aims
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
# 1. Internally we ALWAYS USE eV and angstroms. This is just because these are the
# native units of FHI-aims and yes, I know it's incompatible with esc_lib. It just
# saves useless conversions.
#
#
################################################################################

from __future__ import division
from numpy import *
import h5py

#-------------------------------------------------------------------------------
#
# Data section
#
#-------------------------------------------------------------------------------

# Periodic table dictionary linking symbol to Z number.
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
            116 : "Uuh", 117 : "Uus", 118 : "Uuo" }

#-------------------------------------------------------------------------------
#
# Function Definitions
#
#-------------------------------------------------------------------------------

def indexLine(text, substring, returnAll=False):
  """ index = indexLine(text, substring, returnAll=False)
  
  Returns the index of the line where a substring first occurs
  in a list of strings. If returnAll=True, returns a list of indices
  of all occurrences. 
  
  Note we return None if the string is not found at all.
  
  """
  
  indices = [i for i,line in enumerate(text) if substring in line]
  if returnAll:
    if len(indices) == 0:
      return None
    else:
      return indices
  else:
    # The substring might not occur at all, in which case we return an empty
    # list to be consistent with the returnAll=True case.
    if len(indices) == 0:
      return None
    else:
      return indices[0]

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
            print "(getElementZ) Warning: Element %s is not in the elements dictionary. Returning 0 for element Z." % elstr
            return 0
        else:
            for key, value in elements.items():
                if elstr.title() == value:
                    return key
#-------------------------------------------------------------------------------
#
# Class definitions
#
#-------------------------------------------------------------------------------

class AIMSError(Exception):
  """ Exception raised for errors within the aims_lib code. """
    
  def __init__(self, expr, msg):
    self.expr = expr
    self.msg = msg
    print expr
    print msg

class AimsOutput:
  """ AimsOutput class documentation
  
  This class is used to read information contained in a FHI-aims output file.
  
  Instantiate using the filename: my_aimsoutput = aims_lib.AimsOutput("some_filename.out")
  
  """
  
  text = None
  nspecies = 0
  natoms = 0
  charge = 0.0
  nspins = 1
  nstates = 0
  eigenvalues = None
  occupancies = None
  positions = None
  species = None
  lattice = None
  
  def __init__(self, filename):
    """ Constructor for AimsOutput, pass the filename/path to the FHI-aims output:
    
    my_aimsoutput = aims_lib.AimsOutput("some_filename.out")
    
    """
    
    f = open(filename, 'r')
    
    if f:
      self.text = f.readlines()
      self.readBasicModel()
      self.readGeometry()
      self.readKSEigenvalues()
      f.close()
    else:
      raise AIMSError("AimsOutput.__init__", "File %s could not be opened for reading." % filename)
    
  def readBasicModel(self):
    """ readBasicModel
    
    Grabs things like the number of atoms, species, charge, spin channels, things that
    are always in an aims.out file.
    
    """
    
    try:
      # All of these MUST be present in any Aims output. If they aren't, we have to stop
      # here.
      idx = indexLine(self.text, "  | Number of species")
      self.nspecies = int(self.text[idx].split(":")[-1])
      idx = indexLine(self.text, "  | Number of atoms")
      self.natoms = int(self.text[idx].split(":")[-1])
      idx = indexLine(self.text, "  | Number of spin channels")
      self.nspins = int(self.text[idx].split(":")[-1])    
    except TypeError:
      raise AIMSError("AimsOutput.readBasicModel", "A fundamental parameter (nspecies, natoms or nspins) is missing from the file. It appears this Aims output is broken.")
    idx = indexLine(self.text, "Charge =")
    # Need to be careful now, user might not have specified the charge explicitly.
    if idx is not None:
      self.charge = float(self.text[idx].split(":")[0].split()[2])
    else:
      # Assume 0 charge if not specified.
      self.charge = 0.0
    # Read total number of KS states.
    idx = indexLine(self.text, "| Number of Kohn-Sham states (occupied + empty):")
    self.nstates  = int(self.text[idx].split(":")[-1])
  
  def readGeometry(self):
    """ readGeometry
    
    If there is a final atomic structure listed (ie, there was a geometry optimization)
    then this is read. If not, the initial atomic structure is read.
    
    """
    
    self.species = zeros((self.natoms))
    self.positions = zeros((self.natoms, 3))

    # Have to play games here because there may or may not be a unit cell and this
    # changes both the input and output formats for geometry. Strategy is to find out
    # about the Unit Cell first and set a flag.
    idx = indexLine(self.text, "| No unit cell requested.")
    if idx is None:
      hasUnitCell = True
      self.lattice = zeros((3,3))
    else:
      hasUnitCell = False
    
    # Now look for a final atomic structure and read it if we can.
    idx = indexLine(self.text, "Final atomic structure:")
    if idx is not None:
      if hasUnitCell:
        for i in range(idx + 2, idx + 5):
          bits = self.text[i].split()
          self.lattice[i - idx - 2,:] = array([float(x) for x in bits[1:4]])
        for i in range(idx + 6, idx + 6 + self.natoms):
          bits = self.text[i].split()
          self.species[i - idx - 6] = getElementZ(bits[4])
          self.positions[i - idx - 6,:] = array([float(x) for x in bits[1:4]])          
      else:
        for i in range(idx + 2, idx + 2 + self.natoms):
          bits = self.text[i].split()
          self.species[i - idx - 2] = getElementZ(bits[4])
          self.positions[i - idx - 2,:] = array([float(x) for x in bits[1:4]])
    else:
      # There was no "Final atomic structure" block - look for the inputs instead.
      idx = indexLine(self.text, "| Atomic structure:")
      if idx is not None:
        for i in range(idx + 2, idx + 2 + self.natoms):
          bits = self.text[i].split()
          self.species[i - idx - 2] = getElementZ(bits[3])
          self.positions[i - idx - 2,:] = array([float(x) for x in bits[4:7]])
      else:
        raise AIMSError("AimsOutput.readGeometry", "Apparently there isn't even an input geometry specified in this file. Something is very wrong!")
      # Get unit cell if there is one.
      idx = indexLine(self.text, "| Unit cell:")
      if idx is not None:
        for i in range(idx + 1, idx + 4):
          bits = self.text[i].split()
          self.lattice[i - idx - 1,:] = array([float(x) for x in bits[1:4]])
        
  def readKSEigenvalues(self):
    """ readKSEigenvalues
    
    Reads the last set of eigenvalues listed in the aims.out file. Also populates the 
    occupancies array, both have slots [spin, state].
    
    """
    
    idx = indexLine(self.text, "Writing Kohn-Sham eigenvalues.", returnAll=True)
    if idx is not None:
      # We only want the last place.
      idx = idx[-1]
      self.eigenvalues = zeros((self.nspins, self.nstates))
      self.occupancies = zeros((self.nspins, self.nstates))
      # If we have a unit cell, there is an extra kpoint line in the output
      # that needs to be accounted for.
      if self.lattice is not None:
        offset = 1
      else:
        offset = 0
      if self.nspins == 1:
        for i in range(idx + 3 + offset, idx + 3 + offset + self.nstates):
          bits = self.text[i].split()
          self.occupancies[0,i - idx - 3 - offset] = float(bits[1])
          self.eigenvalues[0,i - idx - 3 - offset] = float(bits[3])
      else:
        # Need to do twice, one for each spin.
        for i in range(idx + 5 + offset, idx + 5 + offset + self.nstates):
          bits = self.text[i].split()
          self.occupancies[0,i - idx - 5 - offset] = float(bits[1])
          self.eigenvalues[0,i - idx - 5 - offset] = float(bits[3])
        idx += 5 + self.nstates
        for i in range(idx + 5 + offset, idx + 5 + offset + self.nstates):
          bits = self.text[i].split()
          self.occupancies[1,i - idx - 5 - offset] = float(bits[1])
          self.eigenvalues[1,i - idx - 5 - offset] = float(bits[3])        
    else:
      raise AIMSError("AimsOutput.readKSEigenvalues", "No eigenvalues present in output file - something is very wrong!")
      
class AimsMomentumMatrixHDF5:
  """ AimsMomentumMatrixHDF5 class documentation
  
  This class is used to read the hdf5-formatted momentum matrix output available from 
  version 061913 of FHI-aims and onwards. 
  
  Instantiate using the filename: my_aimsmommat = aims_lib.AimsMomentumMatrixHDF5("some_filename.h5")
  
  Why do we need this class? Because the momentum matrix elements are stored in a list that runs like:
  
  <1 | del | 2>
  <1 | del | 3>
  ...
  <2 | del | 3>
  ...
  <3 | del | 4>
  ...
  
  so it is like a handshake problem. If we are doing a spectrum calculation and, for
  example, want ALL the transitions <i | del | j>, if i > 1 we need to find the transitions
  <1 | del | i>, for example, and take their complex conjugate. So, this class
  exists mainly to provide energy-ordered lists of momentum matrix elements for a given
  transition.
  
  """
  
  nstates = 0
  nkx = 0
  nky = 0
  nkz = 0
  kpoints = None # Once initialized has shape [nkx, nky, nkz, 4]
  eigenvalues = None # Once initialized has shape [nkx, nky, nkz, nstates]
  matrix_elements = None # Once initialized has shape [nkx, nky, nkz, (nstates+1)*nstates/2, 6] 
  
  
  def __init__(self, filename):
    """ Constructor for AimsMomentumMatrixHDF5, pass the filename/path to the .h5 file containing
    the matrix elements:
    
    my_aimsmommat = aims_lib.AimsMomentumMatrixHDF5("some_filename.h5")
    
    """
    
    f = h5py.File(filename, 'r')
    
    # Copy all the relevant data into our object.
    self.nstates = f['E_bands'].shape[3]
    self.nkx = f['k_points'].shape[0]
    self.nky = f['k_points'].shape[1]
    self.nkz = f['k_points'].shape[2]
    self.kpoints = f['k_points'][:].copy()
    self.eigenvalues = f['E_bands'][:].copy()
    self.matrix_elements = f['Momentummatrix'][:].copy()
    
    f.close()
    
  def getElementsWithInitial(self, i, kx, ky, kz):
    """ mat_elements = getElementsWithInitial(i, kx, ky, kz)
    
    Returns a [nstates - 1, 6]-shape matrix with all the elements <i| del |j> 
    for j != i and the specified k-point.
    
    Note that i is 0-based, not 1-based.
    
    """
    
    # The trick here is that  for all j > i we can just read the matrix elements
    # directly from matrix_elements once we figure out an appropriate set of indices.
    # For j < i, we need to hunt for the bits and take the complex conjugates.
    
    mat_elements = zeros((self.nstates - 1, 6))
    transition_idx = 0
    total_idx = 0
        
    for idx_i in range(0, self.nstates):
      for idx_j in range(idx_i + 1, self.nstates):
        if idx_i == i:
          # Read directly 
          mat_elements[transition_idx,:] = self.matrix_elements[kx, ky, kz, total_idx, :]
          transition_idx += 1
        elif idx_j == i:
          # Complex conjugates.
          mat_elements[transition_idx,0] = self.matrix_elements[kx, ky, kz, total_idx, 0]
          mat_elements[transition_idx,1] = -1.0 * self.matrix_elements[kx, ky, kz, total_idx, 1]
          mat_elements[transition_idx,2] = self.matrix_elements[kx, ky, kz, total_idx, 2]
          mat_elements[transition_idx,3] = -1.0 * self.matrix_elements[kx, ky, kz, total_idx, 3]
          mat_elements[transition_idx,4] = self.matrix_elements[kx, ky, kz, total_idx, 4]
          mat_elements[transition_idx,5] = -1.0 * self.matrix_elements[kx, ky, kz, total_idx, 5]
          transition_idx += 1
        # We only increment transition_idx if we actually are on a transition. 
        # Otherwise just increment the total_idx.
        total_idx += 1
    
    return mat_elements
    
class AimsMomentumMatrixText:
  """ AimsMomentumMatrixText class documentation
  
  This class is used to read the momentum matrix output text files available from 
  version 061913 of FHI-aims and onwards. 
  
  Instantiate using the filename: my_aimsmommat = aims_lib.AimsMomentumMatrixText("some_filename.dat")
  
  Why do we need this class? Because the momentum matrix elements are stored in a list that runs like:
  
  <1 | del | 2>
  <1 | del | 3>
  ...
  <2 | del | 3>
  ...
  <3 | del | 4>
  ...
  
  so it is like a handshake problem. If we are doing a spectrum calculation and, for
  example, want ALL the transitions <i | del | j>, if i > 1 we need to find the transitions
  <1 | del | i>, for example, and take their complex conjugate. So, this class
  exists mainly to provide energy-ordered lists of momentum matrix elements for a given
  transition.
  
  """
  
  nstates = 0
  kpoint = None
  kpoint_idx = 0
  eigenvalues = None 
  matrix_elements = None # Once initialized has shape [(nstates+1)*nstates/2, 6] 
  
  
  def __init__(self, filename):
    """ Constructor for AimsMomentumMatrixText, pass the filename/path to the .h5 file containing
    the matrix elements:
    
    my_aimsmommat = aims_lib.AimsMomentumMatrixText("some_filename.dat")
    
    """
    
    f = open(filename, 'r')
    text = f.readlines()
    f.close()
    
    self.kpoint_idx = int(text[0].split()[2])
    self.kpoint = array([float(x) for x in text[0].split()[4:7]])
    
    text = text[6:]
    
    self.nstates = int((-1 + sqrt(1 + 8.0 * len(text))) / 2)
    self.eigenvalues = zeros((self.nstates))
    self.matrix_elements = zeros((len(text), 6))
    
    transition_idx = 0
    total_idx = 0

    # First pass: read the whole file and just pull the matrix elements.
    for i in range(len(text)):
      bits = text[i].split()
      self.matrix_elements[i,:] = array([float(x) for x in bits[4:10]])
    
    # Second pass: read only the first nstates-1 lines to pull all the eigenvalues.
    for i in range(self.nstates - 1):
      bits = text[i].split()
      if i == 0:
        self.eigenvalues[0] = float(bits[1])
        self.eigenvalues[1] = float(bits[3])
      else:
        self.eigenvalues[i+1] = float(bits[3])
    
  def getElementsWithInitial(self, i):
    """ mat_elements = getElementsWithInitial(i)
    
    Returns a [nstates - 1, 6]-shape matrix with all the elements <i| del |j> 
    for j != i.
    
    Note that i is 0-based, not 1-based.
    
    """
    
    # The trick here is that  for all j > i we can just read the matrix elements
    # directly from matrix_elements once we figure out an appropriate set of indices.
    # For j < i, we need to hunt for the bits and take the complex conjugates.
    
    mat_elements = zeros((self.nstates - 1, 6))
    transition_idx = 0
    total_idx = 0
        
    for idx_i in range(0, self.nstates):
      for idx_j in range(idx_i + 1, self.nstates):
        if idx_i == i:
          # Read directly 
          mat_elements[transition_idx,:] = self.matrix_elements[total_idx, :]
          transition_idx += 1
        elif idx_j == i:
          # Complex conjugates.
          mat_elements[transition_idx,0] = self.matrix_elements[total_idx, 0]
          mat_elements[transition_idx,1] = -1.0 * self.matrix_elements[total_idx, 1]
          mat_elements[transition_idx,2] = self.matrix_elements[total_idx, 2]
          mat_elements[transition_idx,3] = -1.0 * self.matrix_elements[total_idx, 3]
          mat_elements[transition_idx,4] = self.matrix_elements[total_idx, 4]
          mat_elements[transition_idx,5] = -1.0 * self.matrix_elements[total_idx, 5]
          transition_idx += 1
        # We only increment transition_idx if we actually are on a transition. 
        # Otherwise just increment the total_idx.
        total_idx += 1
    
    return mat_elements