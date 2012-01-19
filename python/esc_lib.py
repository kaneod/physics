# Library of electronic-structure related code.
# 
# Note that internally we *always* use atomic units
# (bohr, hartree, etc) and provide converter routines
# to deal with other common units such as ang and eV.

from numpy import array, zeros
from numpy.linalg import norm

# Element dictionaries

elements = { 1 : "H", 6 : "C", 7 : "N" }
xsf_keywords = ["ANIMSTEPS", "CRYSTAL", "ATOMS", "PRIMVEC", "PRIMCOORD"]
bond_lengths = {"CH" : 2.06, "CC" : 2.91, "NC" : 2.78, "NH" : 1.91, "HH" : 2.27}

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
                    
def getBondLengths(positions, species, cutoff=3.0, give_species=False):
        """ bonds = getBondLengths(positions, species, cutoff=3.0, give_species=False)
        
        Returns a list of bond specs [i, j, length] for all pairwise
        distances less than the cutoff distance (default: 3.0 Bohr) for
        the specified animation step (default: first step, ie 0).
        
        If give_species=True, the bond spec includes the species abbreviation:
        
        [i, Zi, j, Zj, length]
        
        """
        
        bonds = []
        for i in range(len(positions)):
            for j in range(i, len(positions)):
                if i is not j:
                    pair = norm(positions[i] - positions[j])
                    if pair < cutoff:
                        if give_species:
                            bonds.append([i,elements[species[i]], j, elements[species[j]], pair])
                        else:
                            bonds.append([i, j, pair])
        
        return bonds

def write_xsf(filename, positions, species, lattice=None):
    """ succeeded = write_xsf(filename, positions, species, lattice=None)
    
    Writes a XSF file containing the passed information. Can be animated,
    in crystal or molecular format, fixed cell or variable cell.
    
    NOTE: NEED TO ADD OPTIONAL FORCES!
    
    """
    
    # Convert everything back to angstroms for XSF
    apos = bohr2ang(positions)
    if lattice is not None:
        alat = bohr2ang(lattice)
    else:
        alat = None
    
    f = open(filename, 'w')
    
    if len(apos) > 1:
        f.write("ANIMSTEPS %d\n" % len(apos))
    
    if alat is not None:
        f.write("CRYSTAL\n")
    
    if alat is not None and len(alat) is 1:
        f.write("PRIMVEC\n")
        f.write("    %g    %g    %g\n" % (alat[0][0][0], alat[0][0][1], alat[0][0][2]))
        f.write("    %g    %g    %g\n" % (alat[0][1][0], alat[0][1][1], alat[0][1][2]))
        f.write("    %g    %g    %g\n" % (alat[0][2][0], alat[0][2][1], alat[0][2][2]))
    
    for i in range(len(apos)):
        if alat is None:
            f.write("ATOMS %d\n" % (i+1))
        if alat is not None and len(alat) > 1:
            f.write("PRIMVEC %d\n" % (i+1))
            f.write("    %g    %g    %g\n" % (alat[i][0][0], alat[i][0][1], alat[i][0][2]))
            f.write("    %g    %g    %g\n" % (alat[i][1][0], alat[i][1][1], alat[i][1][2]))
            f.write("    %g    %g    %g\n" % (alat[i][2][0], alat[i][2][1], alat[i][2][2]))
            f.write("PRIMCOORD %d\n" % (i+1))
            f.write("%d 1\n" % len(apos[i]))
        else:
            f.write("PRIMCOORD %d\n" % (i+1))
            f.write("%d 1\n" % len(apos[i]))
        for j in range(len(apos[i])):
            f.write("%d    %g    %g    %g\n" % (species[i][j], apos[i][j][0], apos[i][j][1], apos[i][j][2]))
        
    f.close()
    return True

def write_abinit(filename, positions, species=None, xtype="bohr"):
    """ succeeded = write_abinit(filename, positions, species=None, xtype="ang")
    
    Writes the passed positions in a format suitable to be copy and pasted
    into an abinit input file. If species are passed, the natom, ntypat,
    typat and znucl parts are also output. Options for xtype are "ang", for
    xangst output (Default) or "bohr" for xcart output.
    
    Note: NEED TO ADD xred OUTPUT!
    
    """
    
    f = open(filename, 'w')
    
    if species is not None:
        f.write("natom       %d\n" % len(species))
        f.write("ntypat      %d\n" % len(uniqify(species)))
        f.write("znucl       %s\n" % " ".join([str(x) for x in uniqify(species)]))
        # Generate typat string
        spec_dict = {}
        typat = []
        for i,s in enumerate(uniqify(species)):
            spec_dict[s] = i+1
        for s in species:
            typat.append(str(spec_dict[s]))
        f.write("typat       %s\n" % " ".join(typat))
    if xtype is "bohr":
        f.write("xcart\n")
        for p in positions:
            f.write("    %010e %010e %010e\n" % (p[0], p[1], p[2]))
    if xtype is "ang":
        f.write("xangst\n")
        for p in bohr2ang(positions):
            f.write("    %010e %010e %010e\n" % (p[0], p[1], p[2]))
    
    f.close()
    return True        
        
        
def bohr2ang(bohr):
    """ ang = bohr2ang(bohr)
    
    Converts bohr units to angstrom, conversion factor 1 bohr = 0.52917721092 ang
    
    Woohoo recursive function!
    
    """
    
    if type(bohr) is type([]):
        # Call on each element and return
        return [bohr2ang(x) for x in bohr]
    else:
        return bohr * 0.52917721092
    
def ang2bohr(ang):
    """ bohr = ang2bohr(ang)
    
    Converts angstroms to bohr, conversion factor 1 ang = 1.0 / 0.52917721092 bohr
    
    """
    
    if type(ang) is type([]):
        return [ang2bohr(x) for x in ang]
    else:
        return ang / 0.52917721092

def eV2hartree(eV):
    """ hartree = eV2hartree(eV)
    
    Converts eV to hartree, conversion factor 1 Ht = 27.21138505 eV
    
    """
    
    if type(eV) is type([]):
        return [eV2hartree(x) for x in eV]
    else:
        return eV / 27.21138505
    
def hartree2eV(hartree):
    """ eV = hartree2eV(hartree)
    
    Converts hartree to eV, conversion factor 1 eV = 1.0 / 27.21138505 Ht
    
    """
    
    if type(hartree) is type([]):
        return [hartree2eV(x) for x in hartree]
    else:
        return hartree * 27.21138505

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
                a = ang2bohr(array([float(x) for x in data[s+1].split()[0:3]]))
                b = ang2bohr(array([float(x) for x in data[s+2].split()[0:3]])) 
                c = ang2bohr(array([float(x) for x in data[s+3].split()[0:3]]))
                self.lattice.append([a, b, c])
            if kw is "PRIMCOORD":
                nat = int(data[s+1].split()[0])
                positions = []
                forces = []
                species = []
                for j in range(nat):
                    bits = data[s+2+j].split()
                    species.append(getElementZ(bits[0]))
                    positions.append(ang2bohr(array([float(x) for x in bits[1:4]])))
                    try:
                        forces.append(ang2bohr(array([float(x) for x in bits[4:7]])))
                    except (ValueError, IndexError):
                        forces.append(array([0.0, 0.0, 0.0]))
                self.positions.append(positions)
                self.forces.append(forces)
                self.species.append(species)
            if kw is "ATOMS":
                # THIS SECTION IS BUGGY!
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
                    positions.append(ang2bohr(array([float(x) for x in bits[1:4]])))
                    try:
                        forces.append(ang2bohr(array([float(x) for x in bits[4:7]])))
                    except (ValueError, IndexError):
                        forces.append(array([0.0, 0.0, 0.0]))
                self.positions.append(positions)
                self.forces.append(forces)
                self.species.append(species)
                
    def getBondLengths(self, cutoff=3.0, animstep=0, give_species=False):
        """ bonds = Atoms.getBondLengths(cutoff=3.0, animstep=1, give_species=False)
        
        Returns a list of bond specs [i, j, length] for all pairwise
        distances less than the cutoff distance (default: 3.0 Bohr) for
        the specified animation step (default: first step, ie 0).
        
        If give_species=True, the bond spec includes the species abbreviation:
        
        [i, Zi, j, Zj, length]
        
        """
        
        return getBondLengths(self.positions[animstep], self.species[animstep], cutoff, give_species)
        
    def autoUnLars(self, animstep=0, convtol=0.01, maxsteps=50, cutoff=2.2):
        """ new_pos, pos_hist = Atoms.autoUnLars(animstep=0, convtol=0.01, maxsteps=50, cutoff=2.2)
        
        Attempts to fix positions generated by Lars automatically. This
        might not work very well. Converges positions to convtol.
        
        """
        
        cur_pos = self.positions[animstep]
        pos_hist = [cur_pos]
        species = self.species[animstep]
        bonds = self.getBondLengths(cutoff=cutoff, animstep=animstep, give_species=True)
        
        for itercount in range(0, maxsteps):
            
            shifts = len(cur_pos) * [array([0.0, 0.0,0.0])]
            print shifts
            
            for bond in bonds:
                # Figure out what the bond length *should* be:
                s1 = species[bond[0]]
                s2 = species[bond[2]]
                if s1 > s2:
                    btype = elements[s1] + elements[s2]
                else:
                    btype = elements[s2] + elements[s1]
                    
                if btype not in bond_lengths.keys():
                    print "%s bond not accounted for." % btype
                else:
                    bproper = bond_lengths[btype]
                    bvector = cur_pos[bond[0]] - cur_pos[bond[2]]
                    bactual = bond[4]
                    shift_factor = bproper / bactual
                    shifts[bond[0]] = shifts[bond[0]] + 0.5 * (shift_factor - 1.0) * bvector
                    shifts[bond[2]] = shifts[bond[2]] + -0.5 * (shift_factor - 1.0) * bvector
                    print "%s bond. Length should be %g, actual is %g. Shift factor is %g." % (btype, bproper, bactual, shift_factor)
            
            # Move all atoms.
            cur_pos = [c+s for (c,s) in zip(cur_pos, shifts)]
            pos_hist.append(cur_pos)
            
            # Check convergence
            diff = 0
            for shift in shifts:
                diff += norm(shift)
                
            if diff < convtol:
                return cur_pos, pos_hist
                
            # Update our bond lengths
            bonds = getBondLengths(cur_pos, species, cutoff=cutoff, give_species=True)
        
        
        return cur_pos, pos_hist
        
            
