# Library of electronic-structure related code.
# 
# Note that internally we *always* use atomic units
# (bohr, hartree, etc) and provide converter routines
# to deal with other common units such as ang and eV.

from numpy import *

def bohr2ang(bohr):
    """ ang = bohr2ang(bohr)
    
    Converts bohr units to angstrom, conversion factor 1 bohr = 0.52917721092 ang
    
    """
    
    return bohr * 0.52917721092
    
def ang2bohr(ang):
    """ bohr = ang2bohr(ang)
    
    Converts angstroms to bohr, conversion factor 1 ang = 1.0 / 0.52917721092 bohr
    
    """
    
    return ang / 0.52917721092

def eV2hartree(eV)
    """ hartree = eV2hartree(eV)
    
    Converts eV to hartree, conversion factor 1 Ht = 27.21138505 eV
    
    """
    
    return eV / 27.21138505
    
def hartree2eV(hartree)
    """ eV = hartree2eV(hartree)
    
    Converts hartree to eV, conversion factor 1 eV = 1.0 / 27.21138505 Ht
    
    """
    
    return hartree * 27.21138505
    
def uniqify(sequence)
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
        
class ESCError(Exception):
    """ Exception raised for errors within the esc_lib code. """
    
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg
    
class Structure:
    """ Superclass representing a structure (crystal, molecule, etc).
    
    Remember - all internal variables are in Hartree atomic units.
    
    """
    
    species = []
    positions_cart = []             # Elements expected to be numpy arrays.
    
    def __init__(self, species, positions, units="bohr"):
        """ struct = Structure(species, positions, units="bohr")
        
        Returns a Structure object just by setting species and positions
        lists directly. Only checks the lengths of each list are the same.
        
        """
        
        if len(species) != len(positions):
            raise ESCError("Total atoms in species (%d) does not match \
            the number of positions (%d)" % (len(species), len(positions)))
        else:
            self.species = species
            if units == "ang":
                self.positions = [ang2bohr(x) for x in positions]
            else:
                self.positions = positions
        
    
    def __init__(self, species_spec, positions, units="bohr"):
        """ struct = Structure(species_spec, positions, units="bohr")
        
        Returns a Structure object generated from a species specification
        and a list of positions. Each element of the species_spec list contains
        a list of two elements:
        
        species_spec = [[n1, Z1], [n2, Z2], ...]
        
        The positions are then assigned such that the first n1 positions are
        of species Z1 and so on.
        
        If units is "ang", the positions are converted to bohr before storage.
        
        """
        
        self.set_species_spec(species_spec, positions, units)
                
    
    def set_species_spec(self, species_spec, positions, units="bohr"):
        """ Structure.set_species_spec(species_spec, positions, units="bohr")
        
        Sets a Structure object to a species specification
        and a list of positions. Each element of the species_spec list contains
        a list of two elements:
        
        species_spec = [[n1, Z1], [n2, Z2], ...]
        
        The positions are then assigned such that the first n1 positions are
        of species Z1 and so on.
        
        If units is "ang", the positions are converted to bohr before storage.
        
        """
        
        
        # Check correct sized inputs.
        total_atoms=0
        for spec in species_spec:
            total_atoms += spec[0]
        
        if total_atoms != len(positions):
            raise ESCError("Total atoms in species_spec (%d) does not match \
            the number of positions (%d)" % (total_atoms, len(positions)))
        else:
            self.species = []
            self.positions_cart = []
            for spec in species_spec:
                self.species += spec[0] * [spec[1]]
            if units=="ang":
                self.positions_cart = [bohr2ang(x) for x in positions]  
                  
    
    def __init__(self, position_spec, units="bohr"):
        """ struct = Structure(position_spec, units="bohr")
        
        Returns a Structure object generated from a position specification. 
        Each element of position_spec is a list:
        
        position_spec = [[Z1, array(R1)], ...]
        
        If units is "ang", the positions are converted to bohr before storage.
        
        """
        
        self.set_position_spec(position_spec, units)
                
    def set_position_spec(self, position_spec, units="bohr"):
        """
        
        Sets the structure to the given position spec.
        Each element of position_spec is a list:
        
        position_spec = [[Z1, array(R1)], ...]
        
        If units is "ang", the positions are converted to bohr before storage.
        
        """
        
        self.species = []
        self.positions_cart = []
        
        for pos in position_spec:
            self.species.append(pos[0])
            if units=="ang":
                self.positions_cart.append(ang2bohr(pos[1]))
            else:
                self.positions_cart.append(pos[1])
          
            
    def ang_positions(self):
        return [bohr2ang(x) for x in self.positions_cart]
    
class Crystal(Structure):
    """ A crystal is just a structure along with a lattice.
    
    """
    
    lattice_cart = []
    
    def __init__(self, lattice, units="bohr"):
        """ crystal = Crystal(lattice, units="bohr")
        
        Generates a (empty) Crystal object using the specified lattice.
        The lattice object should be 
    

