# Library of electronic-structure related code.
# 
# Note that internally we *always* use atomic units
# (bohr, hartree, etc) and provide converter routines
# to deal with other common units such as ang and eV.

from numpy import array

# Element dictionaries

elements = { 1 : "H", 6 : "C", 7 : "N" }

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
        else:
            lines = f.readlines()
            f.close()
            nocoms = []
            # First remove all comment and blank lines
            for line in lines:
                if not (line.strip().startswith('#') or line.strip() == ""):
                    nocoms.append(line)
            
            # Decide if we are animated and whether we are a crystal
            # or a molecule.
            start = 0
            if nocoms[start].split()[0] == "ANIMSTEPS":
                self.nsteps = int(nocoms[0].split()[1])
                start = 1
            
            if nocoms[start].split()[0] == "CRYSTAL":
                self.is_crystal = True
                start = 2
            else:
                start = 1
                           
            # It is not true that if we're not a crystal we're a molecule,
            # however we will assume it is since the other cases (SLAB, etc
            # ) are not seen very often and we can always change later.
            
            j = start
            i = 1
            
            if self.is_crystal:
                # Loop over lines and pull out each animation step.
                # Note that we don't care about CONVVEC or CONVCOORDS
                # here - if PRIMVEC/PRIMCOORD aren't specified, we
                # aren't playing.
                while j < len(nocoms):
                    bits = nocoms[j].split()
                    if bits[0] == "PRIMVEC":
                        # Check we're in step with the current animation step.
                        if self.nsteps > 1:
                            if int(bits[1]) is not i:
                                raise ESCError("Crystal Init", "It looks like the animation steps are out of order in the XSF file, or there are missing steps. Exiting...")
                        # Read and store the lattice coords.
                        print "Got here!"
                        a = ang2bohr(array([float(x) for x in nocoms[j+1].split()[0:3]]))          
                        b = ang2bohr(array([float(x) for x in nocoms[j+2].split()[0:3]]))
                        c = ang2bohr(array([float(x) for x in nocoms[j+3].split()[0:3]]))
                        self.lattice.append([a, b, c])
                        j += 4
                    elif bits[0] == "PRIMCOORD":
                        # Check we're in sync.
                        if self.nsteps > 1:
                            if int(bits[1]) is not i:
                                raise ESCError("Crystal Init", "It looks like the animation steps are out of order in the XSF file, or there are missing steps. Exiting...")
                        # Read number of atoms and use that to grab the coords
                        # Also we might possibly get forces included here.
                        jstep = int(nocoms[j+1].split()[0])
                        curspecies = []
                        curpos = []
                        curforce = []
                        for k in range(jstep):
                            dline = nocoms[j+2+k].split()
                            z = dline[0]
                            curspecies.append(getElementZ(z))
                            pos = ang2bohr(array([float(x) for x in dline[1:4]]))
                            curpos.append(pos)
                            # Try to read forces
                            try:
                                # Note: XSF forces are in angstrom / Ha (??) so 
                                # must convert only the angstrom part.
                                force = ang2bohr(array([float(x) for x in dline[4:7]]))
                                curforce.append(force)
                            except (ValueError, NameError):
                                curforce.append(array([0.0, 0.0, 0.0]))
                            
                        self.positions.append(curpos)
                        self.forces.append(curforce)
                        self.species.append(curspecies)
                        j += jstep + 1
                        i += 1
                    else:
                        j += 1
            else:
                # We aren't a crystal. There is no atom index, we just have to
                # chomp lines until we catch another ATOMS marker.
                
                curspecies = []
                curpos = []
                curforce = []
                
                while j < len(nocoms):
                    bits = nocoms[j].split()
                    if bits[0] == "ATOMS":
                        # Check animation step
                        if int(bits[1]) is not i:
                            raise ESCError("Molecule Init", "It looks like the animation steps are out of order in the XSF file, or there are missing steps. Exiting...")
                        if curspecies:
                            self.species.append(curspecies)
                            self.positions.append(curpos)
                            self.forces.append(curforce)
                        curspecies = []
                        curpos = []
                        curforce = []
                        i += 1
                        j += 1
                    else:
                        z = bits[0]
                        curspecies.append(getElementZ(z))
                        pos = ang2bohr(array([float(x) for x in bits[1:3]]))
                        curpos.append(pos)
                        try:
                            force = ang2bohr(array([float(x) for x in bits[4:7]]))
                            curforce.append(force)
                        except (ValueError, NameError):
                            curforce.append(array([0.0, 0.0, 0.0]))
                        
                        j += 1
    

