# Library of electronic-structure related code.
# 
# Note that internally we *always* use atomic units
# (bohr, hartree, etc) and provide converter routines
# to deal with other common units such as ang and eV.

from numpy import array

# Element dictionaries

elements = { 1 : "H", 6 : "C", 7 : "N" }
xsf_keywords = ["ANIMSTEPS", "CRYSTAL", "ATOMS", "PRIMVEC", "PRIMCOORD"]

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
            raise ESCError("getELementZ", "Element %s is not in the elements dictionary. Returning -1." % elstr)
            return -1
        else:
            for key, value in elements.items():
                if elstr.title() == value:
                    return key

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

def eV2hartree(eV):
    """ hartree = eV2hartree(eV)
    
    Converts eV to hartree, conversion factor 1 Ht = 27.21138505 eV
    
    """
    
    return eV / 27.21138505
    
def hartree2eV(hartree):
    """ eV = hartree2eV(hartree)
    
    Converts hartree to eV, conversion factor 1 eV = 1.0 / 27.21138505 Ht
    
    """
    
    return hartree * 27.21138505
    
def uniqify(sequence):
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
        print expr
        print msg
    
class Atoms:
    """ Atoms class: a collection of atoms (possibly with a crystal structure)
    
    Create using Atoms(xsf_file), Atoms(xyz_file), Atoms(abinit_input), etc.
    
    """
    
    nsteps = 1                  # > 1 if object has time-series data
                                # (ie is animated)
    is_crystal = False          # True if we expect to have lattice vectors
    lattice = []                # List of lists of 3 vectors.
    positions = []              # List of lists of atomic positions.
    forces = []                 # Same for forces...
    species = []                # and species.
                               
    
    def __init__(self, xsf_file):
        """ atoms = Atoms(xsf_file)
        
        Creates an Atoms object from an xsf file. Note that we can deal with
        animated (AXSF) files just fine.
        
        """
        
        # Destroy our current data
        self.is_crystal = False
        self.lattice = []
        self.positions = []
        self.forces = []
        self.species = []
        
        f = open(xsf_file)
        if not f:
            raise ESCError("File %s could not be opened - exiting.")
        
        lines = f.readlines()
        f.close()
        data = []
        # First remove all comment and blank lines
        for line in lines:
            if not (line.strip().startswith('#') or line.strip() == ""):
                data.append(line)
        
        keywords = []
        blocks = []
                
        # Locate all keywords
        for i, line in enumerate(data):
            bits = line.split()
            for kw in xsf_keywords:
                if kw in bits:
                    keywords.append(kw)
                    blocks.append(i)
        
        # Cycle through the keywords and deal with each block.
        for i, (s, kw) in enumerate(zip(blocks, keywords)):
            if kw is "ANIMSTEPS":
                self.nsteps = int(data[s].split()[1])
            if kw is "CRYSTAL":
                self.is_crystal = True
            if kw is "PRIMVEC":
                a = array([float(x) for x in data[s+1].split()[0:3]])
                b = array([float(x) for x in data[s+2].split()[0:3]]) 
                c = array([float(x) for x in data[s+3].split()[0:3]])
                self.lattice.append([a, b, c])
            if kw is "PRIMCOORD":
                nat = int(data[s+1].split()[0])
                positions = []
                forces = []
                species = []
                for j in range(nat):
                    bits = data[s+2+j].split()
                    species.append(getElementZ(bits[0]))
                    positions.append(array([float(x) for x in bits[1:4]]))
                    try:
                        forces.append(array([float(x) for x in bits[4:7]]))
                    except (ValueError, IndexError):
                        forces.append(array([0.0, 0.0, 0.0]))
                self.positions.append(positions)
                self.forces.append(forces)
                self.species.append(species)
            if kw is "ATOMS":
                positions = []
                forces = []
                species = []
                try:
                    s1 = blocks[i+1]
                except (IndexError):
                    s1 = len(data)
                for j in range(s+1, s1):
                    bits = data[j].split()
                    species.append(getElementZ(bits[0]))
                    positions.append(array([float(x) for x in bits[1:4]]))
                    try:
                        forces.append(array([float(x) for x in bits[4:7]]))
                    except (ValueError, IndexError):
                        forces.append(array([0.0, 0.0, 0.0]))
                self.positions.append(positions)
                self.forces.append(forces)
                self.species.append(species)