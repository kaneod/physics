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

def indexLine(text, substring, returnAll=False):
  """ indexLine(text, substring, returnAll=False)
  
  Returns the index of the line where a substring first occurs
  in a list of strings. If returnAll=True, returns a list of indices
  of all occurrences.
  
  """
  
  indices = [i for i,line in enumerate(text) if substring in line]
  if returnAll:
    return indices
  else:
    # The substring might not occur at all, in which case we return an empty
    # list to be consistent with the returnAll=True case.
    if len(indices) == 0:
      return []
    else:
      return indices[0]

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
  
  text = []
  nspecies = 0
  natoms = 0
  charge = 0.0
  nspins = 1
  nstates = 0
  eigenvalues = []
  occupancies = []
  positions = []
  species = []
  
  def __init__(self, filename):
    """ Constructor for AimsOutput, pass the filename/path to the FHI-aims output:
    
    my_aimsoutput = aims_lib.AimsOutput("some_filename.out")
    
    """
    
    f = open(filename, 'r')
    
    if f:
      self.text = f.readlines()
      self.readBasicModel()
      self.readKSEigenvalues()
      f.close()
    else:
      raise AIMSError("AimsOutput.__init__", "File %s could not be opened for reading." % filename)
    
  def readBasicModel(self):
    """ readBasicModel
    
    Grabs things like the number of atoms, species, charge, spin channels, things that
    are always in an aims.out file.
    
    """
    
    idx = indexLine(self.text, "  | Number of species")
    self.nspecies = int(self.text[idx].split(":")[-1])
    idx = indexLine(self.text, "  | Number of atoms")
    self.natoms = int(self.text[idx].split(":")[-1])
    idx = indexLine(self.text, "  | Number of spin channels")
    self.nspins = int(self.text[idx].split(":")[-1])    
    idx = indexLine(self.text, "Charge =")
    # Need to be careful now, user might not have specified the charge explicitly.
    if type(idx) is not type([]):
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
    idx = indexLine(self.text, "Final atomic structure:")
    if type(idx) is not type([]):
      for i in range(idx + 1, idx + 1 + self.natoms):
        bits = self.text[i].split()
        self.species[i - idx - 1] = getElementZ(bits[4])
        self.positions[i - idx - 1,:] = array([float(x) for x in bits[1:4]
        
  def readKSEigenvalues(self):
    """ readKSEigenvalues
    
    Reads the last set of eigenvalues listed in the aims.out file. Also populates the 
    occupancies array, both have slots [spin, state].
    
    """
    
    idx = indexLine(self.text, "Writing Kohn-Sham eigenvalues.", returnAll=True)
    if idx != []:
      # We only want the last place.
      idx = idx[-1]
      self.eigenvalues = zeros((self.nspins, self.nstates))
      self.occupancies = zeros((self.nspins, self.nstates))
      if self.nspins == 1:
        for i in range(idx + 3, idx + 3 + self.nstates):
          bits = self.text[i].split()
          self.occupancies[0,i - idx - 3] = float(bits[1])
          self.eigenvalues[0,i - idx - 3] = float(bits[3])
      else:
        # Need to do twice, one for each spin.
        for i in range(idx + 5, idx + 5 + self.nstates):
          bits = self.text[i].split()
          self.occupancies[0,i - idx - 5] = float(bits[1])
          self.eigenvalues[0,i - idx - 5] = float(bits[3])
        idx += 5 + self.nstates
        for i in range(idx + 4, idx + 4 + self.nstates):
          bits = self.text[i].split()
          self.occupancies[1,i - idx - 4] = float(bits[1])
          self.eigenvalues[1,i - idx - 4] = float(bits[3])        
    else:
      raise AIMSError("AimsOutput.readKSEigenvalues", "No eigenvalues present in output file - something is very wrong!")