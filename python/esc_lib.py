# Library of electronic-structure related code.
# 
# Note that internally we *always* use atomic units
# (bohr, hartree, etc) and provide converter routines
# to deal with other common units such as ang and eV.
from __future__ import division
from numpy import array, zeros, sqrt, reshape, mat
from numpy.linalg import norm
from libabitools import io

# Element dictionaries

elements = { 1 : "H", 2 : "He", 3 : "Li",  6 : "C", 7 : "N" , 8 : "O"}
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

def remove_comments(lines, comment_delim):
    """ stripped = remove_comments(lines, comment_delim)
    
    Takes a sequence of lines (presumably from a data file)
    and strips all comments, including ones at the end of
    lines that are otherwise not comments. Note that we can
    only deal with one kind of comment at a time - just apply
    multiple times to strip multiple comment types.
    
    Note that we also remove *blank* lines in here, just in
    case we're going to do line-by-line processing subsequently
    rather than join/splitting (eg. like for XSF files).
    
    """
    
    stripped = []
    for line in lines:
        if not (line.strip().startswith(comment_delim) or line.strip() == ""):
            stripped.append(line.partition(comment_delim)[0].strip())
    
    return stripped
    
def abinit_parse(istr):
    """ result = abinit_parse(istr)
    
    Tries to parse a string according to abinit input file
    rules. There are only three symbols (*, / and sqrt).
    Returns a result list, possible results are "single",
    for a single data value contained in result[1], "multiple",
    for a list of data values contained in result[1],
    or "fill", for a single data value in result[1] that
    is meant to fill a whole array.
    
    """
    
    # Note the use of the {'sqrt' : sqrt} in the eval
    # statements: this is to restrict the namespace
    # available in the eval function so we don't get
    # a namespace overlap which would happily eval
    # a named variable, causing a logic error.
    
    if "*" in istr.lower():
        factor = istr.lower().partition("*")[0].strip()
        try:
            operand = float(eval(istr.lower().partition("*")[2], {'sqrt' : sqrt}))
        except NameError:
            return "nonsense", istr.lower().partition("*")[2]
        if factor is "":
            return "fill", operand
        else:
            return "multiple", int(factor) * [operand]
    else:
        try: 
            return "single", float(eval(istr.lower(), {'sqrt' : sqrt}))
        except NameError:
            return "non-value", istr.lower()
                
def abinit_value(data, keyword, nvalues):
    """ values = abinit_value(data, keyword, nvalues)
    
    Using the raw abinit data stream, parses the value of the
    given keyword. Because abinit allows rudimentary mathematical
    expressions in the input file, we have to actually parse data
    in most cases rather than just read and convert it. Darn. 
    
    Because we have to parse, we *must* know how many values
    we are looking for. Note that we are also merciless about
    whitespace here - if abinit would misread the input file,
    so will we.
    
    Note that we always read numerical values as floats - 
    
    """
    
    print "Reading in variable %s, looking for %d values." % (keyword, nvalues)
    values = []
    # The action starts at the index of the keyword.
    try:
        start = data.index(keyword)
    except ValueError:
        print "abinit_value WARNING: keyword %s is not in the specified data stream." % keyword
        return None
    
    # Since we don't know how many bits will unpack into
    # the required nvalues items, we need to possible loop
    # over all of them...
    for i, d in enumerate(data[start+1:]):
        try:
            val = float(d)
            values.append(val)
        except ValueError:
            # Need to parse this one.
            result = abinit_parse(d)
            if result[0] == "single":
                values.append(result[1])
            elif result[0] == "multiple":
                values = values + result[1]
            elif result[0] == "fill":
                remaining = nvalues - len(values)
                values = values + remaining * [result[1]]
                return values
        # Note "is" does NOT work here - need ==. Why?
        if len(values) == nvalues:
            return values
            
    # If we get to here, we must not have enough values. Better raise an error.
    

def abinit_unit(data, keyword):
    """ unit = abinit_unit(data, keyword)
    
    Returns the next non-value item after the appearance
    of keyword. At the moment recognizes:
    
    "ry", "angstr", "angstrom", "k", "ev", "ha", "bohr"
    "hartree"
    
    and is case insensitive. 
    
    If a unit is recognized it is returned in the form used
    in this library (ang, eV, etc), otherwise None is returned.
    
    """
    
    try:
        start = data.index(keyword)
    except ValueError:
        print "abinit_unit WARNING: keyword %s is not in the specified data stream." % keyword
        return None
    
    for i, d in enumerate(data[start+1:]):
        print "Checking to see if %s is a unit..." % d
        result = abinit_parse(d)
        if result[0] == "non-value":
            print "%s might be a unit!" % d
            if result[1] in ["angstrom", "angstr"]:
                return "ang"
            elif result[1] == "ev":
                return "eV"
            elif result[1] == "bohr":
                return "bohr"
            elif result[1] in ["ha", "hartree"]:
                return "hartree"
            elif result[1] == "k":
                return "K"
            elif result[1] == "ry":
                return "Ry"
            else:
                return None
    return None

def abinit_int(data, keyword, nvalues):
    """ values = abinit_int(data, keyword, nvalues)
    
    This is just a convenience function that converts a list
    of floats returned by abinit_value into ints. If the list
    only has one member, it returns the integer rather than 
    the list.
    
    """
    
    vals = abinit_value(data, keyword, nvalues)
    
    if len(vals) == 1:
        return int(vals[0])
    else:
        return [int(x) for x in vals]
        
def abinit_array(data, keyword, nvalues, newshape=None):
    """ values = abinit_array(data, keyword, nvalues, newshape=None)
    
    Convenience function that wraps abinit_value but
    returns a numpy array shaped according to the newshape input.
    
    If newshape is left blank, the array is shaped to be n x 3
    (rows x columns) where n = nvalues / 3
    """
    
    vals = abinit_value(data, keyword, nvalues)
    
    # SNEAKY! array(None) is NOT NONE, so have to test
    # prior to conversion. Bloody Numpy!
    if vals is None:
        return None
    else:
        vals = array(vals)
    
    if newshape is None:
        return reshape(vals, (-1, 3))
    else:
        return reshape(vals, newshape)
                
            
def chop128(in_string):
    """ out_string = chop128(in_string)
    
    The abinit input file format requires that each individual
    line be less than 132 characters long otherwise it ignores
    the rest of the input.
    
    NOTE: Not clear if we need to code this yet, waiting for
    confirmation from abinit devs. For now doesn't do anything.
    
    """
    
    return in_string
    
def write_cube_density(filename, positions, species, lattice, densities, timestep=0):
    """ succeeded = write_cube_density(filename, positions, species, lattice, density, timestep=0)
    
    Writes a CUBE file containing the passed information. Note that we
    can't deal with non-crystals so the lattice variable must be passed.
    
    The CUBE file format requires 3D volumetric data, in this case the density.
    
    Since we can't animate, we have to specify a timestep to use. By default this is the
    first timestep.
    
    """
    
    pos = positions[timestep]
    den = densities[timestep]
    spec = species[timestep]
    lat = lattice[timestep]
    
    f = open(filename, 'w')
    # First two lines of a CUBE file are comments.
    f.write("CUBE\n")
    f.write("Output created by esc_lib.\n")
    # Number of atoms then the origin of the density coordinate system.
    f.write("%d 0.000000 0.000000 0.000000\n" % len(spec))
    # Each of the next three lines specifies the number of grid points in an vector
    # direction, then gives the step vector itself.
    nx = den.shape[0]
    ny = den.shape[1]
    nz = den.shape[2]
    f.write("%d %f %f %f\n" % (nx, lat[0][0] / nx, lat[0][1] / nx, lat[0][2] / nx))
    f.write("%d %f %f %f\n" % (ny, lat[1][0] / ny, lat[1][1] / ny, lat[1][2] / ny))
    f.write("%d %f %f %f\n" % (nz, lat[2][0] / nz, lat[2][1] / nz, lat[2][2] / nz))
    # Now list atomic number, charge, position (absolute) for each atom. We don't
    # actually deal with charge here, so set to 0.
    for p, s in zip(pos, spec):
        f.write("%d 0.000000 %f %f %f\n" % (s, p[0], p[1], p[2]))
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                f.write("%f " % den[i,j,k])
                # Throw in a newline every now and then to keep the output readable.
                if k % 6 == 5:
                    f.write("\n")
            f.write("\n")
    
    f.close()
    return True
    
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
    
    if alat is not None and len(alat) == 1:
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
    if xtype == "bohr":
        f.write("xcart\n")
        for p in positions:
            f.write("    %010e %010e %010e\n" % (p[0], p[1], p[2]))
    if xtype == "ang":
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
    densities = []              # List of 3D arrays
    
    test = []                   # Test output
    
    def __init__(self, filename, filetype="XSF"):
        """ atoms = Atoms(filename, filetype="XSF")
        
        Creates an Atoms object from some input file. You have to
        specify the type of file here. Can be:
        
        "XSF" : XCrysden Structure Format, can also be animated axsf.
        
        "abinit" : Abinit input file.
        
        """
        
        if filetype == "XSF":
            self.loadFromXSF(filename)
        elif filetype == "abinit":
            self.loadFromAbinit(filename)
        elif filetype == "abi_density":
            self.loadFromAbinitDensity(filename)
    
    def loadFromAbinitDensity(dens_file):
        """ atoms= Atoms.loadFromAbinitDensity(dens_file)
        
        Internal, inits an Atoms object from a _DEN file written
        by abinit.
        
        """
        
        self.is_crystal = True
        self.lattice = []
        self.positions = []
        self.forces = []
        self.species = []
        self.densities = []
        
        # The libabi2py fortran library does the hard work
        # here, we just need to run the density routine and
        # read out the values into the appropriate python
        # variables.
        
        io.density(dens_file)
        
        nx = io.ngfft[0]
        ny = io.ngfft[1]
        nz = io.ngfft[2]
        
        dens = zeros((nx, ny, nz))
        
        for k in range(nz):
          for j in range(ny):
            for i in range(nx):
              dens[i,j,k] = io.rhor[i + nx * j + nx * ny * k]
        
        dens2 = io.rhor.reshape((nx, ny, nz), order='F')
        
        self.lattice = [[array(x) for x in io.rprimd.tolist()]] 
        self.positions = [[array(x) for x in io.xred.T.tolist()]]
        self.species = [[int(io.znucltypat[x-1]) for x in io.typat]]
        self.densities = [dens]
                                            
    def loadFromAbinit(self, abinit_input):
        """ atoms = Atoms.loadFromAbinit(abinit_input)
        
        Internal, inits an Atoms object from an abinit input file. We
        only allow default values for rprim - input files must contain:
        
        natom, ntypat, typat, znucl, acell, [xangst, xred, xcart]
        
        """
        
        # Destroy current data
        self.is_crystal = True      # Always true for abinit
        self.lattice = []
        self.positions = []
        self.forces = []
        self.species = []
        
        f = open(abinit_input, 'r')
        lines = f.readlines()
        f.close()
        
        data = remove_comments(lines, "#")
        data = remove_comments(lines, "!")
        data = " ".join(data).split()
        
        # Alrighty. Get the number and type of atoms, then generate our species
        # list.
        n = abinit_int(data, "natom", 1)
        ntypat = abinit_int(data, "ntypat", 1)
        znucl = abinit_int(data, "znucl", ntypat)
        # znucl has to be a list
        if type(znucl) == type(1):
            znucl = [znucl]
        typat = abinit_int(data, "typat", n)
        
        # Double brackets because species, positions, etc are all timestep-friendly.
        self.species = [[znucl[x-1] for x in typat]]
        
        acell = abinit_value(data, "acell", 3)
        if acell is None:
            acell = array([1.0, 1.0, 1.0])
        else:
            acell_unit = abinit_unit(data, "acell")
            print acell_unit
            if acell_unit == "ang":
                acell = ang2bohr(acell)
        scalecart = abinit_value(data, "scalecart", 3)
        if scalecart is None:
            scalecart = array([1.0, 1.0, 1.0])        
        rprim = abinit_array(data, "rprim", 9)
        if rprim is None:
            rprim = array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape((3,3))
        
        lat = zeros((3,3))
        for i in [0, 1, 2]:
            lat[i,:] = rprim[i,:] * acell[i]
            print rprim[i,:], acell[i], lat[i,:]
        for i in [0, 1, 2]:
            lat[:,i] = lat[:,i] * scalecart[i]
               
        self.lattice = [lat]
        
        # Try to get positions in the order: xred, xcart, xangst
        pos = abinit_array(data, "xred", n * 3)
        if pos is None:
            pos = abinit_array(data, "xcart", n * 3)
            if pos is None:
                pos = abinit_array(data, "xangst", n * 3)
                if pos is None:
                    raise ESCError("loadFromAbinit", "ERROR: Must have at least one of xred, xcart or xangst specified in the abinit input file.")
                else:
                    pos = ang2bohr(pos)
        else:
            # Use variable lat to convert xred to actual positions.
            pos = array(mat(pos) * mat(lat))
            
        self.positions = [[array(x) for x in pos.tolist()]]
        
    def loadFromXSF(self, xsf_file):
        """ atoms = Atoms.loadFromXSF(xsf_file)
        
        Internal, inits an Atoms object from an xsf file. Note that we can deal with
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
        
        data = remove_comments(lines, "#")
        
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
            if kw == "ANIMSTEPS":
                self.nsteps = int(data[s].split()[1])
            if kw == "CRYSTAL":
                self.is_crystal = True
            if kw == "PRIMVEC":
                a = ang2bohr(array([float(x) for x in data[s+1].split()[0:3]]))
                b = ang2bohr(array([float(x) for x in data[s+2].split()[0:3]])) 
                c = ang2bohr(array([float(x) for x in data[s+3].split()[0:3]]))
                self.lattice.append([a, b, c])
            if kw == "PRIMCOORD":
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
            if kw == "ATOMS":
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
        
            
