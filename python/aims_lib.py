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
from numpy.linalg import *
import h5py
import libaims

#-------------------------------------------------------------------------------
#
# Data section
#
#-------------------------------------------------------------------------------

# Debugging switch - set to True to get debug printouts.
DEBUG = True

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
            
# Conversion factors
hartree2eV = 27.211396132

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
    # The substring might not occur at all, in which case we return None.
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
  nelectrons = 0
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
      self.computeNeutralLUMO()
      self.setArtificialFermiLevel()
      self.computeTransitions()
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
      idx = indexLine(self.text, "The structure contains")
      self.nelectrons = float(self.text[idx].split()[9])
    except TypeError:
      raise AIMSError("AimsOutput.readBasicModel", "A fundamental parameter (nspecies, natoms or nspins) is missing from the file. It appears this Aims output is broken.")
    idx = indexLine(self.text, "Charge =")
    # Need to be careful now, user might not have specified the charge explicitly.
    if idx is not None:
      self.charge = float(self.text[idx].split(":")[0].split()[2])
    else:
      # Assume 0 charge if not specified.
      self.charge = 0.0
    # Read total number of KS states. Need to be careful here, because the number
    # of KS states requested will be bigger than the actual number if it exceeds
    # the number of basis functions. We do a check.
    idx = indexLine(self.text, "| Number of Kohn-Sham states (occupied + empty):")
    self.nstates  = int(self.text[idx].split(":")[-1])
    idx = indexLine(self.text, "| Maximum number of basis functions            :")
    basis_size = int(self.text[idx].split(":")[-1])
    if self.nstates > basis_size:
      self.nstates = basis_size
  
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
    if DEBUG:
      print "Kohn Sham line (1-based): ", idx[-1] + 1
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
        if DEBUG:
          print "First spin read from line (1-based): ", idx + 5 + offset + 1
          print "Occupancy shape: ", self.occupancies.shape
          print "Eigenvalues shape: ", self.eigenvalues.shape
        for i in range(idx + 5 + offset, idx + 5 + offset + self.nstates):
          bits = self.text[i].split()
          self.occupancies[0,i - idx - 5 - offset] = float(bits[1])
          self.eigenvalues[0,i - idx - 5 - offset] = float(bits[3])
        idx += 4 + offset + self.nstates
        if DEBUG:
          print "End of spin read line: ", idx + 1
          print "Second spin read from line (1-based: ", idx + 5 + offset + 1
        for i in range(idx + 5 + offset, idx + 5 + offset + self.nstates):
          bits = self.text[i].split()
          self.occupancies[1,i - idx - 5 - offset] = float(bits[1])
          self.eigenvalues[1,i - idx - 5 - offset] = float(bits[3])     
    else:
      raise AIMSError("AimsOutput.readKSEigenvalues", "No eigenvalues present in output file - something is very wrong!")
  
  def computeTransitions(self):
    """ computeTransitions
    
    Calculates all possible optical (spin-polarized) and direct transitions between KS
    eigenvalues. If nspins=1 the optical and direct transitions are the same. This only
    works for the k-point given in the output file (and hence, generally only useful for
    molecules).
    
    """
    
    self.optical_transitions = zeros((1, self.nspins, self.nstates, self.nstates - 1))
    self.direct_transitions = zeros((1, self.nspins, self.nstates, self.nstates - 1))
    
    if self.nspins == 1:
      for m in range(self.nstates):
        transitions = []
        for j in range(self.eigenvalues.shape[1]):
          if j != m:
            transitions.append(self.eigenvalues[0,j] - self.eigenvalues[0,m])
        self.optical_transitions[0,0,m,:] = array(transitions)
      self.direct_transitions = self.optical_transitions.copy()
    else:
      for s in range(self.nspins):
        for m in range(self.nstates):
          transitions_optical = []
          transitions_direct = []
          for j in range(self.eigenvalues.shape[1]):
            if j != m:
              # Optical transition: take Ej from the opposite spin eigenvalues to Em
              transitions_optical.append(self.eigenvalues[(s+1)%2,j] - self.eigenvalues[s,m])
              # Direct transition: take Ej from the same spin eigenvalues as Em
              transitions_direct.append(self.eigenvalues[s,j] - self.eigenvalues[s,m])
          self.optical_transitions[0,s,m,:] = array(transitions_optical)
          self.direct_transitions[0,s,m,:] = array(transitions_direct)
      
  def computeNeutralLUMO(self):
    """ computeNeutralLUMO
    
    Using the occupancies and the system charge (if any), this function tries to calculate
    the index of the lowest unoccupied KS state and whether it is partially or fully empty
    when the system is NEUTRAL. The reason we do this is because if a core hole has been
    included in the calculation, the reported occupancies reflect the excited system and
    not the neutral system, but we need the neutral LUMO to determine an artificial fermi
    level for molecular NEXAFS calculations. 
    
    Note that this only reads from the occupancies that are output in the main file and
    typically this is just the gamma point. In a system with essentially flat bands, this
    is fine, and in a system with a dense k-point mesh and dispersed bands it probably
    won't matter either because there are plenty of states around the onset of the
    unoccupied states so the spectral weighting at onset will not be substantially 
    affected. YMMV.
    
    """
    
    # First figure out if the system has a core hole by checking for unoccupied states
    # well below the "fermi" level.
    if self.nspins == 1:
      full_occupancy = 2.0
    else:
      full_occupancy = 1.0
    
    first_unoccupied = []
    for i,occ in enumerate(self.occupancies[0,:]):
      if occ < full_occupancy:
        first_unoccupied.append(i)
        break
    if self.nspins == 2:
      for i,occ in enumerate(self.occupancies[1,:]):
        if occ < full_occupancy:
          first_unoccupied.append(i)
          break
    
    print "Found first unoccupied states at i = ", first_unoccupied
    
    # Now decide if these are core holes by looking if there are fully occupied states
    # above.
    self.hole_charge = zeros((self.nspins))
    for i in range(self.nspins):
      if full_occupancy in self.occupancies[i,first_unoccupied[i]+1:]:
        self.hole_charge[i] = full_occupancy - self.occupancies[i,first_unoccupied[i]]
    
    print "Determined we have a hole charge of:", self.hole_charge
    
    # Now find the apparent LUMO (raw occupancies):
    lumo = zeros((self.nspins))
    nelectrons = zeros((self.nspins))
    for i in range(self.nspins):
      nelectrons[i] = sum(self.occupancies[i,:])
      esum = 0.0
      for j in range(self.occupancies.shape[1]):
        if esum == nelectrons[i]:
          lumo[i] = j
          break
        else:
          esum += self.occupancies[i,j]
    
    print "Apparent LUMO for each spin channel is: ", lumo
    
    # If there is no hole charge, or the system charge matches the hole charge, then
    # we are done. 
    if sum(self.hole_charge) == 0.0 or sum(self.hole_charge) == self.charge:
      self.lumo = int(amin(lumo))
      return
    
    # If not, we have an at least partially neutral excitation. Separate cases for 
    # spin polarized and not spin-polarized.
    if self.nspins == 1:
      if full_occupancy - self.occupancies[0,lumo[0]] == self.hole_charge[0]:
        self.lumo = int(amin(lumo))
        return
      elif self.occupancies[0,lumo[0]] == 0.0:
        self.lumo = int(lumo[0] - 1)
        return
      else:
        print "(aims_lib.AimsOutput.computeNeutralLUMO) WARNING: System has a PARTIALLY neutral excitation - we cannot handle these at present. LUMO is probably wrong."
        return
    else:
      if lumo[0] == lumo[1]:
        self.lumo = int(lumo[0] - 1)
        return
      else:
        self.lumo = int(amin(lumo))
        return
  
  def setArtificialFermiLevel(self):
    """ setArtificialFermilevel
    
    Use the neutral LUMO and LUMO - 1 (HOMO) eigenvalues to set a fermi level halfway 
    in between. Note this is meaningless - it simply ensures that transitions from the
    LUMO and up are included in any spectra generated.
    
    """
    
    e_lumo = self.eigenvalues[0,self.lumo]
    e_homo = self.eigenvalues[0,self.lumo-1]
    self.efermi = (e_lumo + e_homo) / 2
    
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
  nkpts = 0
  kpoints = None # Once initialized has shape [nkpts, 4]
  eigenvalues = None # Once initialized has shape [nkpts, nstates]
  matrix_elements = None # Once initialized has shape [nkpts, nstates, nstates - 1, 6]
  
  efermi = None
  
  
  def __init__(self, filename):
    """ Constructor for AimsMomentumMatrixHDF5, pass the filename/path to the .h5 file containing
    the matrix elements:
    
    my_aimsmommat = aims_lib.AimsMomentumMatrixHDF5("some_filename.h5")
    
    """
    
    f = h5py.File(filename, 'r')
    
    # Copy all the relevant data into our object.
    self.nstates = f['E_bands'].shape[3]
    self.ntransitions = self.nstates - 1
    nkx = f['k_points'].shape[0]
    nky = f['k_points'].shape[1]
    nkz = f['k_points'].shape[2]
    self.nkpts = nkx * nky * nkz
    self.efermi = hartree2eV * f['Fermi_energy'][0]   
    
    # The shape of these is a bit silly so we reshape to go by kpt index rather than kpt.
    kpoints = f['k_points'][:]
    eigenvalues = hartree2eV * f['E_bands'][:] # Given in Hartree in the h5 file.
    matrix_elements = f['Momentummatrix'][:]
    
    self.kpoints = zeros((self.nkpts, 3))
    self.matrix_elements = zeros((self.nkpts, matrix_elements.shape[3], 6))
    self.optical_matrix = zeros((self.nkpts, self.nstates, self.ntransitions, 6))
    self.optical_transitions = zeros((self.nkpts, self.nstates, self.ntransitions))
    self.eigenvalues = zeros((self.nkpts, self.nstates))
    
    for i in range(nkx):
      for j in range(nky):
        for k in range(nkz):
          ik = int(kpoints[i,j,k,0])
          self.kpoints[ik,:] = kpoints[i,j,k,1:]
          self.eigenvalues[ik,:] = eigenvalues[i,j,k,:]
          self.matrix_elements[ik,:,:] = matrix_elements[i,j,k,:,:]
          for m in range(self.nstates):
            self.optical_matrix[ik, m, :, :] = self.getElementsWithInitial(m, ik)
            self.optical_transitions[ik, m, :] = self.getTransitionsWithInitial(m, ik)
    f.close()
    
  def getElementsWithInitial(self, i, ik):
    """ mat_elements = getElementsWithInitial(i, ik)
    
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
          mat_elements[transition_idx,:] = self.matrix_elements[ik, total_idx, :]
          transition_idx += 1
        elif idx_j == i:
          # Complex conjugates.
          mat_elements[transition_idx,0] = self.matrix_elements[ik, total_idx, 0]
          mat_elements[transition_idx,1] = -1.0 * self.matrix_elements[ik, total_idx, 1]
          mat_elements[transition_idx,2] = self.matrix_elements[ik, total_idx, 2]
          mat_elements[transition_idx,3] = -1.0 * self.matrix_elements[ik, total_idx, 3]
          mat_elements[transition_idx,4] = self.matrix_elements[ik, total_idx, 4]
          mat_elements[transition_idx,5] = -1.0 * self.matrix_elements[ik, total_idx, 5]
          transition_idx += 1
        # We only increment transition_idx if we actually are on a transition. 
        # Otherwise just increment the total_idx.
        total_idx += 1
    
    return mat_elements
  
  def getTransitionsWithInitial(self, i, ik):
    """ transitions = getTransitionsWithInitial(i, ik)
    
    Returns a [nstates - 1]-shape matrix with the transition energies Ej - Ei for each
    j != i for a single k-point. Note that this uses Ei FROM THE GIVEN K.
    
    Note that i is 0-based, not 1-based.
    
    """ 
    
    transitions = []
    
    for j in range(self.eigenvalues.shape[1]):
      if j != i:
        transitions.append(self.eigenvalues[ik, j] - self.eigenvalues[ik, i])
    
    return array(transitions) 
           

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
  
  The text version (as opposed to the HDF5 version) is slightly different in that it
  only has one kpoint. However, it retains identical array structures to the HDF5 version
  for compatibility.
  
  """
  
  nstates = 0
  nkpts = 1
  #kpoints = None # Once initialized has shape [nkpts, 4]
  eigenvalues = None # Once initialized has shape [1, nstates]
  matrix_elements = None # Once initialized has shape [1, nstates, nstates - 1, 6]
  
  efermi = None
  
  
  def __init__(self, filename):
    """ Constructor for AimsMomentumMatrixText, pass the filename/path to the .h5 file containing
    the matrix elements:
    
    my_aimsmommat = aims_lib.AimsMomentumMatrixText("some_filename.dat")
    
    """
    
    f = open(filename, 'r')
    text = f.readlines()
    f.close()
    
    #self.kpoint_idx = int(text[0].split()[2])
    #self.kpoint = array([float(x) for x in text[0].split()[4:7]])
    
    text = text[6:]
    
    self.nstates = int((-1 + sqrt(1 + 8.0 * len(text))) / 2)
    self.ntrans = self.nstates - 1
    self.eigenvalues = zeros((1, self.nstates))
    self.matrix_elements = zeros((len(text), 6))
    self.optical_matrix = zeros((1, self.nstates, self.ntrans, 6))
    self.optical_transitions = zeros((1, self.nstates, self.ntrans))
    
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
        self.eigenvalues[0,0] = float(bits[1])
        self.eigenvalues[0,1] = float(bits[3])
      else:
        self.eigenvalues[0,i+1] = float(bits[3])
    
    # Third pass: create the optical_matrix and optical_transitions arrays.
    for m in range(self.nstates):
      self.optical_matrix[0, m, :, :] = self.getElementsWithInitial(m)
      self.optical_transitions[0, m, :] = self.getTransitionsWithInitial(m)
    
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

  def getTransitionsWithInitial(self, i):
    """ transitions = getTransitionsWithInitial(i)
    
    Returns a [nstates - 1]-shape matrix with the transition energies Ej - Ei for each
    j != i.
    
    Note that i is 0-based, not 1-based.
    
    """ 
    
    transitions = []
    for j in range(self.eigenvalues.shape[1]):
      if j != i:
        transitions.append(self.eigenvalues[0,j] - self.eigenvalues[0,i])
    
    return array(transitions)
  
class AimsNEXAFS:
  """ AimsNEXAFS class documentation
  
  The purpose of this class is to generate NEXAFS spectra from FHI-aims calculations. It
  can use a variety of inputs that can give increasingly more-useful output. The default 
  constructor:
  
  my_spec = aims_lib.AimsNEXAFS(filename)
  
  assumes that the passed filename points to a HDF5 file and will use a AimsMomentumMatrixHDF5
  object to generate the spectrum. In that case, the system fermi level is used to figure
  out where the empty states start.
  
  If you specify both a .dat file and an output file:
  
  my_spec = aims_lib.AimsNEXAFS(filename, HDF5=False, outputFile=outFilename)
  
  a AimsMomentumMatrixText object will be used along with an AimsOutput object. These
  two are needed together because in this case figuring out where the unoccupied states
  start matters in a more precise way.
  
  """
  
  # These variables are set within the constructors.
  mommatrix = None
  momhdf5 = None
  output = None
  
  # These variables are set from generateSpectrum()
  w =  None # This is the photon energy array itself, shape (num_points)
  spectral_cmpts = None # Individual components of the spectrum
  
  # Set all the following before calling generateSpectrum()
  # If unset, default values are used after printing a warning.
  w_start = None
  w_end = None
  num_points = None
  
  # Smearing factors
  linear_smearing = None
  lorentzian_smearing = None
  gaussian_smearing = None
  
  def __init__(self, filename, HDF5=True, outputFile=None):
    """ my_spec = aims_lib.AimsNEXAFS(filename, HDF5=True, outputFile=None)
    
    Constructor for AimsNEXAFS object.
    
    """
    
    # Instantiate the mommatrix object using the file.
    if HDF5:
      self.mommatrix = AimsMomentumMatrixHDF5(filename)
      self.momhdf5 = True
    else:
      if outputFile:
        self.mommatrix = AimsMomentumMatrixText(filename)
        self.momhdf5 = False
        self.output = AimsOutput(outputFile)
        self.mommatrix.efermi = self.output.efermi
      else:
        raise AimsERROR("aims_lib.AimsNEXAFS.__init__", "If using a non-HDF5 momentum matrix, a FHI-aims output file is also required.")    
   
  def generateSpectrum(self, i_core, override_w=True):
    """ worked = AimsNEXAFS.generateSpectrum(i_core, override_w=True)
    
    Populates the spectral_cmpts matrix using the momentum matrix elements. Note that
    we DO NOT SMEAR here - we set the smearing to half the step size as explained in
    the libpytep code. Smearing has to happen later on.
    
    i_core is the KS index of the state used as the core level.
    
    By default (override_w=True), the w_start and w_end parameters are overwritten
    by new automatically generated values. If override_w=False, the existing ones are
    left in place if they exist at all.
    
    """
    
    # First job is to set up the energy spectrum. If w_start and w_end haven't
    # been specified, we run from -5 eV below the Fermi level up to 5 eV past the highest 
    # transition energy which is amax(eigenvalues) - eigenvalues[icore] + 5.
    
    # We make the (reasonable) assumption here that the core state i_core is not 
    # dispersed with k, so any k choice will do if HDF5 input is used.
    
    if override_w or (not self.w_start):
      self.w_start = self.mommatrix.efermi - self.mommatrix.eigenvalues[0,i_core] - 5.0
    if override_w or (not self.w_end):
      maxe = amax(self.mommatrix.eigenvalues)
      self.w_end = maxe - self.mommatrix.eigenvalues[0,i_core] + 5.0
    
    if not self.num_points:
      self.num_points = 2000
      
    self.w = linspace(self.w_start, self.w_end, self.num_points)
    self.spectral_cmpts = zeros((self.num_points, 6))

    # Zero the spectral_cmpts matrix.
    self.spectral_cmpts = zeros((self.num_points, 6), 'f')
    
    # Due to the arrays being set up the same in both AimsMomentumMatrixHDF5 and
    # AimsMomentumMatrixText, we can use the same function call to libaims here to 
    # generate the spectral components. However, note that the efermi level passed here
    # is ARTIFICIAL for the text case (it is computed by the AimsMomentumMatrixText
    # object and not from FHI-aims itself. 
    self.spectral_cmpts = libaims.spectroscopy.generate_spectrum(self.mommatrix.optical_matrix, \
                                                self.mommatrix.optical_transitions, \
                                                self.mommatrix.eigenvalues, \
                                                self.w, self.mommatrix.efermi, i_core+1)
                                                
  def experimentalSpectrum(self, polarization_vector):
    """ experimentalSpectrum
    
    Constructs a single spectrum (shape [num_points]) for the specified polarization
    vector and using the smearing parameters of the AimsNEXAFS object.
        
    """              
    
    # Set some defaults if they aren't set.
    if not self.linear_smearing:
      self.linear_smearing = 0.001
      
    if not self.lorentzian_smearing:
      self.lorentzian_smearing = 0.1
      
    if not self.gaussian_smearing:
      self.gaussian_smearing = 0.1
      
    # Normalize the polarization vector.
    p = array(polarization_vector) / norm(array(polarization_vector))
    
    # Make sure the spectral components have actually been generated!
    if self.spectral_cmpts is None:
      print "(aims_lib.AimsNEXAFS.experimentalSpectrum) ERROR - you have not yet generated any spectral components. Call generateSpectrum first."
      return
      
    self.spectrum = p[0]**2 * self.spectral_cmpts[:,0] + \
                    p[1]**2 * self.spectral_cmpts[:,1] + \
                    p[2]**2 * self.spectral_cmpts[:,2] + \
                    p[0]*p[1] * self.spectral_cmpts[:,3] + \
                    p[0]*p[2] * self.spectral_cmpts[:,4] + \
                    p[1]*p[2] * self.spectral_cmpts[:,5]
    
    self.spectrum = libaims.utilities.lorentzian_linear_convolute(self.w, \
                                      self.spectrum, self.lorentzian_smearing, \
                                      self.linear_smearing)
    
    self.spectrum = libaims.utilities.gaussian_convolute(self.w, self.spectrum, \
                                      self.gaussian_smearing)