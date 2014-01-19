################################################################################
#
# esc_lib.py
#
# Library of electronic-structure related code
#
################################################################################
#
# Copyright 2012 Kane O'Donnell
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
# 1. Internally we *always* use atomic units (bohr, hartree, etc) and provide 
# converter routines to deal with other common units such as ang and eV.
#
# 2. In this version, a lot of the redundant stuff that can be done well in other
# public codes has been removed, e.g. NEXAFS stuff. If you really want it, go
# to esc_lib.old.py.
#
################################################################################


from __future__ import division
from numpy import array, zeros, sqrt, reshape, mat, pi, matrix, cos, sin, exp, arange, arccos, arctan2, complex, polyfit, poly1d, loadtxt, amin, amax, argmin, argmax, dot
from random import random as rand
from random import choice
from numpy.linalg import norm, inv
from numpy.fft import fftn, ifftn
#from scipy.interpolate import interp1d
#from scipy.integrate import quad
#from scipy.special import sph_harm
#from scipy.optimize import leastsq, curve_fit, fmin_slsqp
#from scipy.io import netcdf
import os

# Debugging flag - set to 1 to see debug messages.
DEBUG=1

# Element dictionaries

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
            print "(libesc.getElementZ) Warning: Element %s is not in the elements dictionary. Returning 0 for element Z." % elstr
            return 0
        else:
            for key, value in elements.items():
                if elstr.title() == value:
                    return key
                    
def vectorAxisComparison(a, b):
  """ result = vectorAxisComparison(a, b)
  
  For sorting lists of vectors. Two vectors are compared first by their z component, 
  then their y component, then their x component. If a < b in this scheme, returns
  a negative number. If a = b, returns 0. If a > b, returns positive number.
  
  """
  
  if a[2] - b[2]:
    comp = a[2] - b[2]
  elif a[1] - b[1]:
    comp = a[1] - b[1]
  elif a[0] - b[0]:
    comp = a[0] - b[0]
  else:
    return 0
  
  if comp < 0:
    return -1
  else:
    return 1
    
def sortListPair(a, b):
  """ a_sorted, b_sorted = sortListPair(a, b)
  
  Sorts two lists of vectors a and b so that the mapping of a[i] to b[i] remains
  the same even if the i changes.
  
  """
  
  if len(a) <= 1:
    return a, b
  else:
    less_a = []
    more_a = []
    less_b = []
    more_b = []
    same_a = []
    same_b = []
    ipiv = choice(range(len(a)))
    pivot = a[ipiv] 
    for va, vb in zip(a,b):
      cmp = vectorAxisComparison(va,pivot)
      if cmp < 0:
        less_a.append(va)
        less_b.append(vb)
      if cmp > 0:
        more_a.append(va)
        more_b.append(vb)
      if cmp == 0:
        same_a.append(va)
        same_b.append(vb)
    less_a, less_b = sortListPair(less_a, less_b)
    more_a, more_b = sortListPair(more_a, more_b)
    return more_a + same_a + less_a, more_b + same_b + less_b
    

def substringInList(substring, list):
  """ line = substringInList(substring, list)
  
  If substring is a substring of any of the lines in the passed list, return the line
  index. If not, return False. First instance only!
  
  """
  
  for i, line in enumerate(list):
    if substring in line:
      return i
  
  return False
  
  
def substringPositionsInList(substring, list):
  """ line_indices = substringPositionsInList(substring, list)
  
  Same as substringInList but returns ALL line indices. Slower because
  it always goes through the whole list.
  
  """
  
  line_indices = []
  for i, line in enumerate(list):
    if substring in line:
      line_indices.append(i)
  
  if len(line_indices) > 0:
    return line_indices
  else:
    return False

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
  
def integral_smoothing(data, period, pbc=False):
  """ smooth_data = integral_smoothing(data, period)
  
  Assumes data is two column with x values data[:,0] and y values data[:,1].
  Then the new y value at x is an integral of the old data within a single
  period centred on x. Period is in the same units as the data x-axis itself.
  
  Set pbs to True (default False) to integrate past the edge and wrap around 
  (ie periodic boundary conditions).
  
  The only condition imposed is that the period must be less than the span
  of x - I haven't written "multi-wrap-around" code here so it can't handle
  multiply-periodic spans.
  
  """
  
  # Check the period is less than the span of x.
  if period > abs(amax(data[:,0]) - amin(data[:,0])):
    print "ERROR (integral_smoothing): Period greater than span of x. This is not allowed."
    return None
  
  # Make an interpolater and a new data array. We also need to know the 
  # bounds of x to use in the case of periodic boundary conditions.
  smooth = zeros(data.shape)
  func = interp1d(data[:,0], data[:,1])
  xmin = amin(data[:,0])
  xmax = amax(data[:,0])
  
  lx = xmin
  rx = xmax
  
  halfp = float(period) / 2.0
  
  # For each point compute the integral.
  for i, x in enumerate(data[:,0]):
    smooth[i,0] = x
    if x - halfp < xmin:
      lx = xmin
    else:
      lx = x - halfp
    if x + halfp > xmax:
      rx = xmax
    else:
      rx = x + halfp
    if DEBUG:
      print "At point %f, integrating between %f and %f." % (x, lx, rx)
    # Integrate the bit between lx and rx
    this_int = quad(func, lx, rx)[0]
    if DEBUG:
      print "Result is %f." % this_int
    # If PBCs are used, need to add any bits that pass the boundaries.
    if pbc is True:
      if x - halfp < xmin:
        lx = xmax - (halfp + xmin - x)
        rx = xmax
        this_int += quad(func, lx, rx)[0]
        if DEBUG:
          print "Added PBC extra bit %f." % (this_int)
      if x + halfp > xmax:
        lx = xmin
        rx = xmin + halfp - xmax + x
        this_int += quad(func, lx, rx)[0]
        if DEBUG:
          print "Added PBC extra bit %f." % (this_int)
    smooth[i,1] = this_int / period
    
  return smooth
  
def read_seq_xy(filename, comment_delim="#", option=None):
  """ data = read_seq_xy(filename, comment_delim="#", option=None)
  
  Reads xy files where the individual data sets are listed in two columns
  sequentially. We assume (but don't check) that all the sets are the same
  size. The returned data has the form:
  
  data[i,j,k]
  
  where i runs over the sets, j runs over the points and k runs over the 
  columns of data, so 0 is the x axis, 1 is the y1 axis etc.
  
  Note that we can accommodate blank lines separating the sets.
  
  If option="castep_elnes" we correct for the silly output format and "k" is 4.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  header_lines = []
  cur_block = []
  data_blocks = []
  
  for i in range(len(lines)):
    try:
      x = float(lines[i].split()[0])
      cur_block.append([float(y) for y in lines[i].split()])
    except ValueError:
      if cur_block != []:
        data_blocks.append(cur_block)
        cur_block = []
      header_lines.append(lines[i])
      
  if cur_block != []:
    data_blocks.append(cur_block)
  
  data_blocks = array(data_blocks)
  
  if option=="castep_elnes":
    npoints = data_blocks.shape[1]
    data = zeros((data_blocks.shape[0], npoints/2, 4))
    for i in range(data_blocks.shape[0]):
      data[i,:,0] = data_blocks[i,0:npoints:2,0]
      data[i,:,1] = data_blocks[i,0:npoints:2,1]
      data[i,:,2] = data_blocks[i,1:npoints:2,0]
      data[i,:,3] = data_blocks[i,1:npoints:2,1]
    
    return data
  else:
    return array(data_blocks)

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
    
def abinit_read_bands(filename):
  """ bands, properties = abinit_read_bands(filename)
  
  Reads in a _EIG output from ABINIT and returns as an array. Also reads in
  the weightings (wtk) and the reduced coords of the kpoints.
  
  """
  
  f = open(filename, 'r')
  
  lines = f.readlines()
  f.close()
  nkpts = int(lines[0].split()[4])
  bands = []
  kpts = []
  wtks = []
  
  curbands = []
  for line in lines[1:]:
    if line.startswith(" kpt#"):
      if curbands != []:
        bands.append(curbands)
      curbands = []
      wtks.append(float(line.strip().split()[5].strip(',')))
      kpts.append([float(x) for x in line.strip().split()[7:10]])
    else:
      curbands = curbands + [float(x) for x in line.split()]
      
  bands.append(curbands)
   
  props = {}
  props["kpts"] = kpts
  props["wtks"] = wtks
   
  return bands, props

def abinit_read_gw_bands(filename):
  """ bands, properties = abinit_read_gw_bands(filename)
  
  Same as abinit_read_bands except reads the eigenvalues from an abinit
  _GW file. Also returns the kpoints and the eigenvalue corrections. Note that
  the _GW file contains the levels in eV. For consistency in accordance with our
  internal code policy we convert this to hartree.
  
  """
  
  f = open(filename, 'r')
  
  lines = f.readlines()
  f.close()
  nkpts = int(lines[0].split()[0])
  kpts = []
  bandindices = []
  bands = []
  corrs = []
  lines = lines[1:]
  
  while len(lines) > 0:
    curline = lines[0]
    kpts.append([float(x) for x in curline.split()])
    nvals = int(lines[1].split()[0])
    lines = lines[2:]
    curband = []
    curbandindices = []
    curcorrs = []
    for i in range(nvals):
      bits = lines[i].split()
      curbandindices += [int(bits[0])]
      curband += [float(bits[1])]
      curcorrs += [float(bits[2])]
    bands.append(curband)
    corrs.append(curcorrs)
    bandindices.append(curbandindices)
    lines = lines[i+1:]
    
  props = {}
  props["indices"] = bandindices
  props["corrs"] = corrs
  props["kpts"] = kpts
  
  bands = eV2hartree(bands)
  return bands, props
   
   
def castep_read_bands(filename):
  """ bands, properties = castep_read_bands(filename)
  
  Reads in the SEED.bands from CASTEP and returns as an array. Also reads
  in some of the key properties, such as the number of bands, number of electrons,
  fermi level.
  
  """
  
  f = open(filename, 'r')
  
  lines = f.readlines()
  nkpts = int(lines[0].split()[3])
  nspins = int(lines[1].split()[4])
  if nspins == 1:
    nelectrons = [float(lines[2].split()[3])]
    nbands = [int(lines[3].split()[3])]
    efermi = [float(lines[4].split()[5])]
  else:
    nelectrons = [float(x) for x in lines[2].split()[3:]]
    nbands = [int(x) for x in lines[3].split()[3:]]
    efermi = [float(x) for x in lines[4].split()[5:]]
  lattice = []
  lattice.append([float(x) for x in lines[6].split()])
  lattice.append([float(x) for x in lines[7].split()])
  lattice.append([float(x) for x in lines[8].split()])
  lattice = array(lattice)
  data = lines[9:]
  
  bandsdict = {}
  bands = []
  kptdict = {}

  for k in range(nkpts):
    # Read the kpoint line because they aren't necessarily in order.
    kptnum = int(data[0].split()[1])
    kpt = array([float(x) for x in data[0].split()[2:5]])
    kptdict[kptnum] = kpt
    data = data[1:]
    for s in range(nspins):
      data = data[1:]
      tmp = data[0:nbands[s]]
      bandsdict[kptnum] = array([float(x.strip()) for x in tmp])
      data = data[nbands[s]:]

  for k in range(nkpts):
    bands.append(bandsdict[k+1])
  
  bands = array(bands)

  props = {}  
  props["nbands"] = nbands
  props["nkpts"] = nkpts
  props["nspins"] = nspins
  props["efermi"] = efermi
  props["nelectrons"] = nelectrons
  props["kpts"] = kptdict
  props["lattice"] = lattice
  
  return bands, props

def castep_generate_band_path(kpoints, keypoints,recip_lattice):
  """ path, keyvalues = generate_band_path(kpoints, keypoints, recip_lattice)
  
  Generates a path in inverse angstrom units that can be used
  as an x-axis for a bandstructure. kpoints can be a dictionary (keyed by an integer
  1-based index, like the props['kpts'] property from castep_read_bands) or a list/array.
  
  keypoints is a list/array of kpoints that act as turning points. The output keyvalues
  indicates the value along the path where the keypoint is or would be located.
  
  You need the reciprocal space lattice to make these distances meaningful - they aren't
  even *relatively* correct if the primitive BZ is not a cube.
  """
  
  # Tolerance for how close two vectors can be to be considered the same.
  tol = 1.0e-6
  
  # If kpoints is a dictionary, make it a list.
  if type(kpoints) == type({}):
    kpts = []
    for k in kpoints.keys():
      kpts.append(kpoints[k])
  elif type(kpoints) == type([]) or type(kpoints) == type(array([])):
    kpts = kpoints
  else:
    print "(castep_generate_band_path) ERROR: Invalid input kpoints - must be array, list or dictionary."
    return None
    
  # Expand all kpts into cartesian coordinates using the recip_lattice.
  for i in range(len(kpts)):
    kpts[i] = kpts[i][0] * recip_lattice[0] + kpts[i][1] * recip_lattice[1] + kpts[i][2] * recip_lattice[2]
  
  # The same for the keypoints.
  for i in range(len(keypoints)):
    keypoints[i] = keypoints[i][0] * recip_lattice[0] + keypoints[i][1] * recip_lattice[1] + keypoints[i][2] * recip_lattice[2]  
  
  # Algorithm is trickier than it appears: idea is to generate linear functions
  # from the keypoints and continually ask if a given kpoint is on the current linear
  # function. If it is, add it's distance along the linear function to the path, otherwise
  # generate a new segment linear function.
  
  # First task is to make segments from the keypoints.
  segments = []
  segment_lengths = []
  cumulative_segment_lengths = []
  seg_starts = []
  seg_ends = []
  k_low = keypoints[0]
  for k in keypoints[1:]:
    segments.append(k - k_low)
    segment_lengths.append(norm(k - k_low))
    cumulative_segment_lengths.append(sum(segment_lengths[:-1]))
    seg_ends.append(k)
    seg_starts.append(k_low)
    k_low = k
  
  if DEBUG:
    print "Segments are:"
    for s in segments:
      print s
    print "Segment lengths are:"
    for s in segment_lengths:
      print s
    print "Cumulative segment lengths are:"
    for s in cumulative_segment_lengths:
      print s
    print "Segment starts are:"
    for s in seg_starts:
      print s
    print "Segment ends are:"
    for s in seg_ends:
      print s

  # Now cycle through the kpoints. Three cases: either it is one of the keypoints,
  # between the two (on the segment) or off the segment. 
  path = []
  keyvalues = []
  current_segment = 0
  current_kpoint = 0
  finished = False
  k = kpts[0]
  while not finished:
    #print "Current k is: ", current_kpoint, k, " Current segment:", current_segment
    if norm(k - seg_starts[current_segment]) < tol:
      # k is the starting point of this segment. Add to path, increment kpt.
      #print "kpoint is the starting point of a segment!"
      path.append(cumulative_segment_lengths[current_segment])
      keyvalues.append(path[-1])
      current_kpoint += 1
    elif norm(k - seg_ends[current_segment]) < tol:
      # k is the ending point of this segment. Add to path, increment kpoint AND segment.
      #print "kpoint is the ending point of a segment!"
      path.append(cumulative_segment_lengths[current_segment] + norm(k - seg_starts[current_segment]))
      keyvalues.append(path[-1])
      current_kpoint += 1
      current_segment += 1
    else:
      # Check if it's on the current segment.
      v = segments[current_segment]/norm(segments[current_segment])
      kklow = (k - seg_starts[current_segment]) / norm(k - seg_starts[current_segment])
      #print "Vector is neither a beginning or endpoint."
      #print "v is: ", v, "kklow is ", kklow, "Norm is: ", norm(v - kklow)
      if norm(v - kklow) < tol:
        # k is on the segment. Add to path, increment kpoint.
        #print "kpoint is on a segment!"
        path.append(cumulative_segment_lengths[current_segment] + norm(k - seg_starts[current_segment]))
        current_kpoint += 1
      else:
        # k is not on the segment. Can't add to path, increment segment.
        #print "kpoint is NOT on this segment."
        # The keyvalue, even though it isn't part of the path, must be the cumulative
        # sum of the path lengths so far, as we calculated earlier.
        current_segment += 1
        keyvalues.append(cumulative_segment_path_lengths[current_segment])
    
    # Check to see if we're finished.
    if current_kpoint == len(kpts):
      finished = True
    else:
      k = kpts[current_kpoint]
  
  return path, keyvalues
        
def elk_read_bands(filename):
  """ path, bands = elk_read_bands(filename="BAND.OUT")
  
  Reads in the BAND.OUT file from Elk (can optionally pass
  a different filename) and returns as an array. The first
  array column is the path parameter and the subsequent columns
  are the eigenvalues.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  sets = []
  curset = []
  curpath = []
  have_path = False
  for line in lines:
    bits = line.split()
    if len(bits) != 2:
      if have_path:
        sets.append(curset)
      else:
        sets.append(curpath)
        sets.append(curset)
        have_path=True
      curset = []
      curpath = []
    else:
      curpath.append(float(bits[0]))
      curset.append(float(bits[1]))
  
  return array(sets).T
  
def abinit_read_dos(filename, properties=True):
  """ dos[, properties] = abinit_read_dos(filename, properties=True)
  
  Reads the _DOS files generated by Abinit. Properties is optional (default
  is True) and if True returns a dictionary of properties read from the _DOS
  file.
  
  The dos return is an array of all the columns in the data section of the file.
  What the columns mean depends on the type of DOS file - the value of
  header["columns"] gives the title row for the file.
  
  Current header contents: nsppol, nkpt, nband, fermi, iat, ratsph, columns.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  header = {}
  datastart = 0
  
  for i,l in enumerate(lines):
    bits = l.split()
    if l.strip().startswith("#"):
      # Line is part of header, search for properties.
      if "nsppol" in l:
        header["nsppol"] = int(bits[3].rstrip(","))
        header["nkpt"] = int(bits[6].rstrip(","))
        if header["nsppol"] == 1:
          header["nband"] = int(bits[8])
        else:
          header["nband"] = (int(bits[8].rstrip(",")), int(bits[10]))
      
      if "Fermi energy" in l:
        header["fermi"] = float(bits[4])
      if "iat=" in l:
        header["iat"] = int(bits[5])
      if "ratsph" in l:
        header["ratsph"] = float(bits[6])
      if "energy(Ha)" in l:
        header["columns"] = l.strip()
    else:
      datastart = i
      break
  
  cols = len(lines[datastart].split())
  rows = len(lines[datastart:])
    
  dos = zeros((rows,cols))
  
  for i,l in enumerate(lines[datastart:]):
    for j,bit in enumerate(l.split()):
      dos[i,j] = float(bit)
      
  if properties:
    return dos, header
  else:
    return dos  

def abinit_collate_dos(files):
  """ dos_sum, column_title, fermi = abinit_collate_dos(files)
  
  Sums together dos components from separate _DOS_ATXXXX files from abinit.
  We assume column 0 is the energy column (not summed) and the remainder are
  all added together in the return array. The output column_title is just
  one of the file header["columns"] values (from the last file, in fact) and
  is used in abinit_write_collated_dos. Same for fermi.
  
  """
  
  multi_dos = []
  column_title = ""
  fermi = 0.0
  
  for filename in files:
     dos, header = abinit_read_dos(filename)
     column_title = header["columns"] # Just take the last one
     fermi = header["fermi"] # and here
     multi_dos.append(dos)
     
  dos_sum = zeros(multi_dos[0].shape)
  
  for dos in multi_dos:
    # If the shapes mismatch, the following is guaranteed to fail.
    dos_sum[:,1:] = dos_sum[:,1:] + dos[:,1:]
  
  # Just grab the energy scale from the first data set.
  dos_sum[:,0] = multi_dos[0][:,0]
  
  return dos_sum, column_title, fermi
  
def abinit_write_collated_dos(files, outfile="collated_dos.xy"):
  """ result = abinit_write_collated_dos(files, outfile="collated_dos.xy")
  
  Writes collated dos to the specified file. Format is just multi-column xy
  with the fermi and column_title outputs used as header lines.
  
  """
  
  f = open(outfile, 'w')
  
  dos, title, fermi = abinit_collate_dos(files)
  
  f.write("# Fermi energy: " + str(fermi) + " Ha\n")
  f.write(title + "\n")
  
  for i in range(dos.shape[0]):
    line = "\t".join([str(x) for x in dos[i,:]]) + "\n"
    f.write(line)
    
  f.close()
  
  return True
      
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
    
    if vals is None:
      return None
    elif len(vals) == 1:
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

def elk_value(data, keyword, nvalues):
    """ values = elk_values(data, keyword, nvalues)
    
    Returns nvalues values from an elk input datastream corresponding to
    the given keyword. Items are returned as strings - use wrapper functions to
    get floats and ints.
    
    """
    
    if keyword in data:
      start = data.index(keyword)+1
      try:
        return [float(x) for x in data[start:start+nvalues]]
      except ValueError:
        # value is not convertible to a float: must be a string instead.
        return data[start:start+nvalues]
    else:
      print "elk_value WARNING: keyword %s does not exist in this file." % keyword
      return None
      
def elk_array(data, keyword, nvalues, newshape=None):
    """ a = elk_array(data, keyword, nvalues, newshape=None)
    
    Convenience function wrapping elk_value but returning floats. Returns a 
    single value if nvalues is 1.
    
    """
    
    vals = elk_value(data, keyword, nvalues)
    
    if vals is None:
      return None
    elif len(vals) == 1:
      return float(vals[0])
    else:
      vals = array(vals)
      
    if newshape is None:
      return vals
    else:
      return vals.reshape(newshape)
      
def elk_int(data, keyword, nvalues):
    """ i = elk_int(data, keyword, nvalues)
    
    Convenience function returning integers. Will return a single int if nvalues
    is one, otherwise returns a *list* of ints (not array type).
    
    """
    
    vals = elk_value(data, keyword, nvalues)
    
    if vals is None:
      return None
    elif len(vals) == 1:
      return int(vals[0])
    else:
      return [int(x) for x in vals]
    
            
def chop128(in_string):
    """ out_string = chop128(in_string)
    
    The abinit input file format requires that each individual
    line be less than 132 characters long otherwise it ignores
    the rest of the input.
    
    """
    tail =  in_string.split('\n')[-1]
    head = in_string.split('\n')[:-1]
    
    if len(tail) <= 128:
      return in_string
    else:
      return "\n".join(head + [tail[0:128]] + [chop128(tail[128:])])
      
    
def write_cube(filename, positions, species, lattice, datagrid, timestep=0):
  """ succeeded = write_cube(filename, positions, species, lattice, datagrid, timestep=0)
  
  Writes a CUBE file containing the passed information. The datagrid must be a 3D array.
  
  We can't animate, so we must specify the timestep. By default this is the first timestep.
  
  """
  pos = positions[timestep]
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
  nx = datagrid.shape[0]
  ny = datagrid.shape[1]
  nz = datagrid.shape[2]
  f.write("%d %g %g %g\n" % (nx, lat[0][0] / nx, lat[0][1] / nx, lat[0][2] / nx))
  f.write("%d %g %g %g\n" % (ny, lat[1][0] / ny, lat[1][1] / ny, lat[1][2] / ny))
  f.write("%d %g %g %g\n" % (nz, lat[2][0] / nz, lat[2][1] / nz, lat[2][2] / nz))
  # Now list atomic number, charge, position (absolute) for each atom. We don't
  # actually deal with charge here, so set to 0.
  for p, s in zip(pos, spec):
    f.write("%d 0.000000 %g %g %g\n" % (s, p[0], p[1], p[2]))
    
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        f.write("%g " % datagrid[i,j,k])
        # Throw in a newline every now and then to keep the output readable.
        if k % 6 == 5:
          f.write("\n")
      f.write("\n")
  
  f.close()
  return True  
  
    
def write_cube_density(filename, positions, species, lattice, densities, timestep=0):
    """ succeeded = write_cube_density(filename, positions, species, lattice, density, timestep=0)
    
    Writes a CUBE file containing the passed information. Note that we
    can't deal with non-crystals so the lattice variable must be passed.
    
    The CUBE file format requires 3D volumetric data, in this case the density.
    
    Since we can't animate, we have to specify a timestep to use. By default this is the
    first timestep.
    
    """
    
    return write_cube(filename, positions, species, lattice, densities[timestep], timestep)
    
def write_xyz(filename, positions, species):
  """ succeeded == write_xyz(filename, positions, species)
  
  Writes an XYZ file, dimensions in angstroms and using the chemical symbol (not the Z
  number) as the identifier on each line. All positions are written so animation is supported.
  
  If len(species) = len(positions), the species array at each timestep is used, otherwise
  species[0] is used for all and there better not be changes in the number of atoms!
  
  """
  
  f = open(filename, 'w')
  
  apos = bohr2ang(positions)
  
  if len(apos) == len(species):
    for step, pos, spec in enumerate(zip(apos, species)):
      f.write(str(len(pos))+"\n")
      f.write("Created by libesc.py - step %d.\n" % step)
      for z, p in zip(spec, pos):
        if z in elements.keys():
          s = elements[z]
        else:
          s = "UKN"
        f.write("%s\t\t%10.5g\t%10.5g\t%10.5g\n" % (s, p[0], p[1], p[2]))
  else:
    spec = species[0]
    for step, pos in enumerate(apos):
      f.write(str(len(pos))+"\n")
      f.write("Created by libesc.py - Animation step %d.\n" % step)
      for z, p in zip(spec, pos):
        if z in elements.keys():
          s = elements[z]
        else:
          s = "UKN"
        f.write("%s\t\t%10.5g\t%10.5g\t%10.5g\n" % (s, p[0], p[1], p[2])) 
  
  f.close()
  return True
  
def write_molden(filename, positions, species, opt={}):
  """ succeeded = write_molden(filename, positions, species, opt={})
  
  Writes a .molden file with any animation sequence put in the [GEOMETRIES] section.
  Eventually there will be an option to output normal coordinates for vibrations.
  
  Options:
    "converged main geometry" : Boolean - if True (default), the FINAL set of positions
    is used for the main Atoms section on the basis that it is the converged geometry.
    If False, the FIRST set of positions is used. 
  
  """
      
  f = open(filename, 'w')
  apos = bohr2ang(positions)
  
  f.write("[Molden Format]\n")
  f.write("[Atoms] Angs\n")
  
  if "converged main geometry" in opt.keys():
    if not opt["converged main geometry"]:
      maingeom = apos[0]
      mainspec = species[0]
    else:
      maingeom = apos[-1]
      mainspec = species[-1]
  else:
    maingeom = apos[-1]
    mainspec = species[-1]    
  
  for j, (z, p) in enumerate(zip(mainspec, maingeom)):
    if z in elements.keys():
      s = elements[z]
    else:
      s = "UKN"
    f.write("%s\t%d\t%d\t%10.5g\t%10.5g\t%10.5g\n" % (s, j+1, z, p[0], p[1], p[2]))
  
  # Geometries
  f.write("[GEOMETRIES] XYZ\n")
  
  if len(apos) == len(species):
    for step, (pos, spec) in enumerate(zip(apos, species)):
      f.write(str(len(pos))+"\n")
      f.write("Created by libesc.py - step %d.\n" % step)
      for z, p in zip(spec, pos):
        if z in elements.keys():
          s = elements[z]
        else:
          s = "UKN"
        f.write("%s\t\t%10.5g\t%10.5g\t%10.5g\n" % (s, p[0], p[1], p[2]))
  else:
    spec = species[0]
    for step, pos in enumerate(apos):
      f.write(str(len(pos))+"\n")
      f.write("Created by libesc.py - Animation step %d.\n" % step)
      for z, p in zip(spec, pos):
        if z in elements.keys():
          s = elements[z]
        else:
          s = "UKN"
        f.write("%s\t\t%10.5g\t%10.5g\t%10.5g\n" % (s, p[0], p[1], p[2]))
  
  f.close()
  
  return True
  
def write_xsf(filename, positions, species, lattice=None, letter_spec=True):
    """ succeeded = write_xsf(filename, positions, species, lattice=None, letter_spec=True)
    
    Writes a XSF file containing the passed information. Can be animated,
    in crystal or molecular format, fixed cell or variable cell.
    
    If letter_spec is True (default), writes the species identifier as a letter
    rather than a nuclear Z value.
    
    NOTE: NEED TO ADD OPTIONAL FORCES!
    
    """
    
    if DEBUG: print len(positions)
    # Convert everything back to angstroms for XSF
    apos = bohr2ang(positions)
    if lattice is not None and lattice != []:
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
        elif alat is not None and len(alat) > 1:
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
          if letter_spec:
            f.write("%s    %g    %g    %g\n" % (elements[species[i][j]], apos[i][j][0], apos[i][j][1], apos[i][j][2]))
          else:
            f.write("%d    %g    %g    %g\n" % (species[i][j], apos[i][j][0], apos[i][j][1], apos[i][j][2]))
        
    f.close()
    return True
    
def write_pdb(filename, positions, species, lattice=None, opt=None, timestep=0):
  """ Write a PDB file.
  
  Some notes about compliance:
  
  1. Uses "ATOM" keyword for all atoms (never HETATM).
  2. If a lattice is passed, this is written in ABC format to CRYST1 but the space
  group is not calculated and is therefore void (by default P 1).
  3. Does not include any bonding information.
  
  """
  
  pos = positions[timestep]
  spec = species[timestep]
  
  pos = bohr2ang(pos)
  print lattice
  if lattice is not None and lattice != []:
    avec = bohr2ang(lattice[timestep])
    is_crystal = True
  else:
    is_crystal = False
    
  f = open(filename, 'w')
  f.write("REMARK written by esc_lib.py (Kane O'Donnell, kaneod on GitHub)\n")
  
  if is_crystal:
    # Get the lengths of the vectors and the angles between them.
    a = norm(avec[0])
    b = norm(avec[1])
    c = norm(avec[2])
    alpha = arccos(abs(dot(avec[0], avec[1]))/(a * b)) * 180 / pi
    beta = arccos(abs(dot(avec[1], avec[2]))/(b * c)) * 180 / pi
    gamma = arccos(abs(dot(avec[2], avec[0]))/(c * a)) * 180 / pi
    
    print "Found lengths and angles:"
    print "%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f" % (a, b, c, alpha, beta, gamma)
    
    f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2fP 1           1\n" % (a, b, c, alpha, beta, gamma))
    
  for i, (p, s) in enumerate(zip(pos,spec)):
    specname = elements[getElementZ(s)].center(4)
    f.write("ATOM  %5d %s   X     1    %8.3f%8.3f%8.3f                      %s  \n" % \
                  (i+1, specname, p[0],p[1],p[2],specname))
  
  f.write("END\n")
  f.close()
  return True

def write_aims(filename, positions, species, lattice, xtype="ang", opt=None, timestep=0):
  """ succeeded = write_aims(filename, positions, species, lattice, xtype="ang", opt=None, timestep=0)

    Writes a FHI-aims geometry.in file using the given positions, species and lattice. 

    FHI-aims allows either periodic or non-periodic boundary conditions. For periodic,
    specify xtype="frac" and provide lattice vectors. For non-periodic, specify xtype=
    "ang". No lattice vectors will be written in that case.
    
    opt is a dictionary where the key is the name of the option and the value is some
    additional data. Options implemented so far are:
    
    "constrain atoms": value is a list of atomic indices to constrain.
    
    "constrain species": value is a list of species to constrain.

  """

  pos = positions[timestep]
  spec = species[timestep]
  
  if "constrain atoms" in opt.keys():
    index_constraints = True
  else:
    index_constraints = False
   
  if "constrain species" in opt.keys():
    species_constraints = True
    constraint_list_species = [getElementZ(x) for x in opt["constrain species"]]
  else:
    species_constraints = False 

  if xtype == "frac":
    avec = lattice[timestep]
    pos = cart2reduced(pos, avec)
    avec = bohr2ang(avec)
  elif xtype == "ang":
    pos = bohr2ang(pos)
  else:
    print "write_aims ERROR: Must specify xtype=ang or frac."
    return False

  f = open(filename, 'w')
  f.write("# geometry.in written by esc_lib.py\n\n")
  if xtype == "ang":  
    for i in range(len(pos)):
      p = pos[i]
      s = spec[i]
      f.write("atom  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], elements[s]))
      if index_constraints and i in opt["constrain atoms"]:
        f.write("  constrain_relaxation .true.\n")
      elif species_constraints and s in constraint_list_species:
        f.write("  constrain_relaxation .true.\n")
  elif xtype == "frac":
    for l in avec:
      f.write("lattice_vector %4.8g %4.8g %4.8g\n" % (l[0], l[1], l[2]))
    
    f.write("\n")
    for s, p in zip(spec, pos):
      f.write("atom_frac  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], elements[s]))
      if index_constraints and i in opt["constrain atoms"]:
        f.write("  constrain_relaxation .true.\n")
      elif species_constraints and s in constraint_list_species:
        f.write("  constrain_relaxation .true.\n")
  f.close()

  return True
    
def write_castep(filename, positions, species, lattice, xtype="ang", opt=None, timestep=0):
  """ succeeded = write_castep(filename, positions, species=None, xtype="bohr", opt=None, timestep=0)
  
  Writes a CASTEP .cell file using the passed positions. Unlike the abinit case,
  species must NOT be None here. Options for xtype are "ang", "bohr" or "frac".
  We assume, as always, that EVERYTHING internal is in atomic units, so a 
  conversion is performed if "frac" is specified using the passed lattice
  vectors. Also, if fractional output is specified, the lattice vectors
  are output in angstroms, the CASTEP default length unit.
  
  opt is a dictionary that gives output options. Options are:
  
  'special atom' : n - index of atom that is special. Will have ":special" appended
                      to the Z number in the positions_abs/frac block so that
                      you can specify a separate species_pot pseudopotential
                      generation string or PSP file.
  
  """
  
  pos = positions[timestep]
  avec = lattice[timestep]
  spec = species[timestep]
  
  # Do conversions if necessary.
  if xtype == "ang":
    pos = bohr2ang(pos)
    avec = bohr2ang(avec)
  elif xtype == "frac":
    pos = cart2reduced(pos,avec)
    avec = bohr2ang(avec)
  
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
  for i, (s, p) in enumerate(zip(spec, pos)):
    if opt is not None and 'special atom' in opt:
      if opt["special atom"] == i+1:
        f.write("  %s %010e %010e %010e\n" % (elements[s]+":special", p[0], p[1], p[2]))
      else:
        f.write("  %s %010e %010e %010e\n" % (elements[s], p[0], p[1], p[2]))
    else:
      f.write("  %s %010e %010e %010e\n" % (elements[s], p[0], p[1], p[2]))
  if xtype == "frac":
    f.write("%endblock positions_frac\n")
  else:
    f.write("%endblock positions_abs\n")
  
  f.close()
  return True
  

def write_abinit(filename, positions, species=None, xtype="bohr", opt=None, timestep=0):
    """ succeeded = write_abinit(filename, positions, species=None, xtype="ang", timestep=0)
    
    Writes the passed positions in a format suitable to be copy and pasted
    into an abinit input file. If species are passed, the natom, ntypat,
    typat and znucl parts are also output. Options for xtype are "ang", for
    xangst output (Default) or "bohr" for xcart output.
    
    One needs to specify a timestep since the positions variable can be animated
    as can the species variable - the default is zero.
    
    opt is a dictionary of special options. Can have:
    
    'special atom' : n (index of a special atom, given a unique label and 
                    moved to the first position in the atoms list. For use
                    with conducti/core hole style calculations where the 
                    first atom is special according to abinit.
    
    Note: NEED TO ADD xred OUTPUT!
    
    """
    
    pos = positions[timestep]
    
    # If we have a special atom, have to swap it to position 0.
    
    if opt is not None and 'special atom' in opt:
      n = opt['special atom']
      pos0 = pos[0]
      posn = pos[n-1]
      pos[0] = posn
      pos[n-1] = pos0
    
    f = open(filename, 'w')
    
    if species is not None:
        spec = species[timestep]
        # If we have special atom, have to swap it's species to the first pos
        # and give it a special element type.
        if opt is not None and 'special atom' in opt:
          n = opt['special atom']
          spec0 = spec[0]
          specn = spec[n-1]
          spec[0] = specn
          spec[n-1] = spec0
          f.write("natom       %d\n" % len(spec))
          f.write("ntypat      %d\n" % (len(uniqify(spec[1:])) + 1))
          f.write("znucl       %s\n" % (str(spec[0])+" "+" ".join([str(x) for x in uniqify(spec[1:])])))
          spec_dict = {}
          typat = []
          for i,s in enumerate(uniqify(spec[1:])):
              spec_dict[s] = i+2
          for s in spec[1:]:
              typat.append(str(spec_dict[s]))
          typatstr = "1 " + " ".join(typat)
        else:         
          f.write("natom       %d\n" % len(spec))
          f.write("ntypat      %d\n" % len(uniqify(spec)))
          f.write("znucl       %s\n" % " ".join([str(x) for x in uniqify(spec)]))
          # Generate typat string
          spec_dict = {}
          typat = []
          for i,s in enumerate(uniqify(spec)):
              spec_dict[s] = i+1
          for s in spec:
              typat.append(str(spec_dict[s]))
          typatstr = " ".join(typat)
        # Abinit complains if a line in an input file > 132 characters in length
        # so we break every 100 characters if necessary.
        if len(typatstr) >= 132:
          typatstr = chop128(typatstr)
        f.write("typat       %s\n" % typatstr)
    if xtype == "bohr":
        f.write("xcart\n")
        for p in pos:
            f.write("    %010e %010e %010e\n" % (p[0], p[1], p[2]))
    if xtype == "ang":
        f.write("xangst\n")
        for p in bohr2ang(pos):
            f.write("    %010e %010e %010e\n" % (p[0], p[1], p[2]))
    
    f.close()
    return True        
        
def write_elk(filename, positions, species, is_crystal=False, lattice=None, timestep=0):
    """ succeeded = write_elk(filename, positions, species, is_crystal=False, lattice= None, timestep=0)
    
    Writes the passed positions and species in a format suitable to be copy and
    pasted into an elk input file. Since elk requires atomic units, we do not
    have a unit conversion option here. If is_crystal is True, the coordinates
    are written as reduced with respect to the specified lattice.
    
    As usual we can accommodate animated data with a timestep parameter, the
    default is zero (first timestep).
    
    Note: MOLECULAR INPUT TO ELK IS NOT WORKING IN 1.4.5
    
    """
    
    pos = positions[timestep]
    spec = species[timestep]
    
    f = open(filename, 'w')
    
    if not is_crystal:
      f.write("molecule\n")
      f.write("  .true.\n\n")
    else:
      if lattice is None:
        raise ESCError("write_elk", "ERROR: If is_crystal is False, you must specify a lattice")
      else:
        lat = lattice[timestep]
        pos = cart2reduced(pos, lat)
        f.write("avec\n")
        f.write("  %g %g %g\n" % (lat[0][0], lat[0][1], lat[0][2]))
        f.write("  %g %g %g\n" % (lat[1][0], lat[1][1], lat[1][2]))
        f.write("  %g %g %g\n" % (lat[2][0], lat[2][1], lat[2][2]))
        f.write("\n")
    
    f.write("atoms\n")
    f.write("  %d\n" % len(uniqify(spec)))          # Total number of atoms
    for s in uniqify(spec):
      f.write("  'REPLACE.in'\n")             # Replace in the output file
                                              # with the species filename
      f.write("  %d\n" % spec.count(s))    # Number of atoms of species s 
      
      for i,p in enumerate(pos):
        if spec[i] == s:
          f.write("  %g %g %g 0.000000 0.000000 0.000000\n" % (p[0], p[1], p[2]))
      
    f.write("\n")    
    f.close()
    return True

def recip_lattice(avec):
  """ bvec = recip_lattice(avec)
  
  Computes the reciprocal lattice vectors for a given crystal lattice. We
  assume avec is in the form generated for the Atoms class - [a1, a2, a3]
  where ai is the ith lattice vector (an array). 
  
  bvec is returned in the same format.
  
  """
    
  A = mat(avec).T # columns of A are the crystal lattice vectors
  B = (2 * pi * inv(A)) # Rows of B are the recip lattice vectors
  return [array(v) for v in B.tolist()]
    
def g_vectors(bvec, ngrid):
  """ gv, [gcx,gcy,gcz], gb = g_vectors(bvec, ngrid)
  
  Returns a list of g-vectors (Nx3 array, fortran order).
  
  gcx/y/z is a list of three arrays containing the actual recip
  vector indices (ie, they go negative as appropriate).
  
  gb is the biggest G-vector norm.
  
  """
  
  gcx = zeros((ngrid[0]))
  gcy = zeros((ngrid[1]))
  gcz = zeros((ngrid[2]))
  
  for i in range(0,int(ngrid[0]/2)+1):
    gcx[i] = i
  for i in range(int(ngrid[0]/2) + 1, ngrid[0]):
    gcx[i] = i - ngrid[0]
    
  for i in range(0,int(ngrid[1]/2)+1):
    gcy[i] = i
  for i in range(int(ngrid[1]/2) + 1, ngrid[1]):
    gcy[i] = i - ngrid[1]
    
  for i in range(0,int(ngrid[2]/2)+1):
    gcz[i] = i
  for i in range(int(ngrid[2]/2) + 1, ngrid[2]):
    gcz[i] = i - ngrid[2]
  
  biggest = 0.0
  gvs = []
  for k in range(ngrid[2]):
    tmpz = gcz[k] * bvec[2]
    for j in range(ngrid[1]):
      tmpy = gcy[j] * bvec[1]
      for i in range(ngrid[0]):
        gv = tmpz + tmpy + gcx[i] * bvec[0]
        if norm(gv) > biggest:
          biggest = norm(gv)
        gvs.append(gv)
        
  
  return array(gvs), [gcx, gcy, gcz], biggest
  
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
    is_crystal = True          # True if we expect to have lattice vectors
    lattice = []                # List of lists of 3 vectors.
    positions = []              # List of lists of atomic positions.
    forces = []                 # Same for forces...
    species = []                # and species.
    densities = []              # List of 3D arrays
    wfks = []
    filehook = 0                # For *really big files* it's better to
                                # keep a filehook and just read as necessary.
                                # An example is any of the larger NetCDF output
                                # files from Abinit.
    ngrid = []                  # size of real-space grid.
    test = []                   # Test output
    
    def clean(self):
      """ Atoms.clean()
      
      Internal - wipes out data in the instance object, DON'T USE!
      
      """
      
      self.is_crystal = True
      self.lattice = []
      self.positions = []
      self.forces = []
      self.species = []
      self.densities = []
      self.wfks = []
      self.ngrid = []
      if self.filehook:
        self.filehook.close()
        self.filehook = 0
    
    def __init__(self, filename, filetype="XSF", options=None):
        """ atoms = Atoms(filename, filetype="XSF")
        
        Creates an Atoms object from some input file. You have to
        specify the type of file here. Can be:
        
        "XSF" : XCrysden Structure Format, can also be animated axsf.
        
        "abinit,input" : Abinit input file.
        
        "abinit,density" : Abinit _DEN file.
        
        "abinit,wfk" : _WFK file from Abinit.
        
        "elk,input" : Elk input file.
        
        "castep,cell" : CASTEP .cell file.
        
        "aims,geometry" : geometry.in file from FHI-aims
        
        "aims,output" : Output file from FHI-aims.
        
        "abinit,NetCDF" : Abinit NetCDF output. Note: this importer is not very clever
        and will put in default values for the species if they cannot be found -
        default is all atoms are carbon.
        
        "ETSF" : Read from ETSF-formatted NetCDF output.
        
        "pdb" : Protein Data Bank format.
        
        """
        
        if filetype == "XSF":
            self.loadFromXSF(filename)
        elif filetype == "abinit,input":
            self.loadFromAbinit(filename)
        elif filetype == "abinit,density":
            self.loadFromAbinitDensity(filename)
        elif filetype == "abinit,wfk":
            self.loadFromAbinitWFK(filename)
        elif filetype == "elk,input":
            self.loadFromElk(filename)
        elif filetype == "abinit,NetCDF":
            self.loadFromNetCDF(filename)
        elif filetype == "ETSF":
            self.loadFromETSF(filename)
        elif filetype == "castep,cell":
            self.loadFromCastep(filename)
        elif filetype == "VASP,poscar":
            self.loadFromVASP(filename, options)
        elif filetype == "aims,geometry":
          self.loadFromAimsGeometry(filename)
        elif filetype == "aims,output":
          self.loadFromAimsOutput(filename)
        elif filetype == "pdb":
          self.loadFromPDB(filename)
        else:
          print "(esc_lib.Atoms.__init__) ERROR: File type %s not handled at present." % filetype
          return None
  
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
            
      # If we have a crystal lattice specified, convert it to vectors. Note that since
      # in a PDB the crystal is specified in ABC format (three lengths, three angles)
      # the conversion is not unique and needs to be done according to a standard. Here,
      # the CASTEP standard is used: "a" goes along x, b is in the xy plane and the whole
      # thing must have a positive cell volume with a.(bxc) [this fixes c].
      if crys is not None:
        lvec = []
        lengths = array([float(x) for x in data[crys].split()[1:4]])
        angles = array([float(x) for x in data[crys].split()[4:7]])
        angles *= pi/180.0
        lvec.append(array([lengths[0], 0.0, 0.0])) # "a"
        lvec.append(array([lengths[1] * cos(angles[2]), lengths[1] * sin(angles[2]), 0.0]))
        ca = lengths[2] * cos(angles[1])
        cb = lengths[2] * (cos(angles[0]) - cos(angles[1]) * cos(angles[2])) / sin(angles[2])
        cc = sqrt(lengths[2] ** 2 - ca ** 2 - cb ** 2)
        lvec.append(array([ca, cb, cc]))
        self.is_crystal = True
        if num_models is not None:
          self.lattice.append(num_models * ang2bohr(lvec))
        else:
          self.lattice.append(ang2bohr(lvec))
      
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
        pos.append(ang2bohr(array([float(x) for x in at[30:54].split()])))
        spec.append(getElementZ(at[76:78].strip()))
      
      if num_models is not None:
        chunksize = int(float(len(pos)) / num_models)
        for j in range(num_models):
          self.positions.append(pos[j*chunksize:(j+1)*chunksize])
          self.species.append(spec[j*chunksize:(j+1)*chunksize])
      else:
        self.positions.append(pos)
        self.species.append(spec)
      
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
          self.positions.append(ang2bohr(pos))
          self.species.append(spec)
          self.is_crystal = False
        else:
          # Molecule is a crystal. Get the lattice vectors then the initial pos and spec.
          idx = substringInList("| Unit cell:", list)
          for line in data[(idx+1):(idx+4)]:
            lvec.append(array([float(x) for x in line.split()[1:4]]))
          for line in data[(idx+6):(idx+6+numat)]:
            spec.append(getElementZ(line.split()[3]))
            pos.append(array([float(x) for x in line.split()[4:7]]))
          self.positions.append(ang2bohr(pos))
          self.species.append(spec)          
          self.lattice.append(ang2bohr(lvec))
          self.is_crystal = True
      else:
        # This out file is broken: Didn't even get to initial calculation.
        print "(Atoms.loadFromAimsOutput) ERROR: File is not a complete FHI-aims output. Returning None."
        return None
      
      # Now get all the coordinate updates.
      idxs = substringPositionsInList("Updated atomic structure:", data)
      for i in idxs:
        pos = []
        for line in data[(i+2):(i+2+numat)]:
          pos.append(array([float(x) for x in line.split()[1:4]]))
        self.positions.append(ang2bohr(pos))
           
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
        self.positions.append(ang2bohr(pos))
        self.species.append(spec)
        self.is_crystal = False
      else:
        lvec = ang2bohr(lvec)
        self.lattice.append(lvec)
        if is_atoms:
          pos = ang2bohr(pos)
          self.positions.append(pos)
        else:
          self.positions.append(reduced2cart(pos, lvec))
        self.species.append(spec)
        self.is_crystal = True
      
    def loadFromCastep(self, filename):
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
              self.lattice.append(ang2bohr([vec1, vec2, vec3]))
              i = i + 4
            elif data[i+5].split()[0] == "%endblock":
              units = data[i+1].split()[0].lower()
              vec1 = array([float(x) for x in data[i+2].split()])
              vec2 = array([float(x) for x in data[i+3].split()])
              vec3 = array([float(x) for x in data[i+4].split()])
              if units == "ang" or units == "angstrom":
                self.lattice.append(ang2bohr([vec1, vec2, vec3]))
              elif units == "bohr":
                self.lattice.append([vec1, vec2, vec3])
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
              self.positions.append(ang2bohr(pos))
            elif unit == "bohr":
              self.positions.append(pos)
            self.species.append(specs)
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
            self.species.append(specs)
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
          
          # Since we don't actually need the values, we just print them.
          if DEBUG:
            print "Found option: ", option, " with value ", value
        i = i + 1
      
      if postype == "fractional":
        self.positions.append(reduced2cart(pos, self.lattice[0]))
      
    def loadFromETSF(self, filename):
      """ atoms = Atoms.loadFromETSF(filename)
      
      Internal, inits an Atoms object from an ETSF NetCDF output. Can get this
      from Abinit using pawprtwf, for example.
      
      """
      
      self.clean()
      
      f = netcdf.netcdf_file(filename, 'r')
      pos = [array(x) for x in f.variables['reduced_atom_positions'].data]
      avec = [array(x) for x in f.variables['primitive_vectors'].data]
      
      self.positions.append(reduced2cart(pos,avec))
      self.lattice.append(avec)
      self.recip_lattice = [recip_lattice(avec)]
      
      n1 = f.dimensions['number_of_grid_points_vector1']
      n2 = f.dimensions['number_of_grid_points_vector2']
      n3 = f.dimensions['number_of_grid_points_vector3']
      self.ngrid = [n1, n2, n3]
      
      # Generate G-vectors (we need them for WF calculations)
      #self.g_vectors, self.g_coords,self.max_g = g_vectors(self.recip_lattice[0],self.ngrid)     
      znucl = f.variables['atomic_numbers'].data
      self.species = [[int(znucl[x-1]) for x in f.variables['atom_species'].data]]
      self.filehook = f
            
    def loadFromNetCDF(self, filename):
      """ atoms = Atoms.loadFromNetCDF(filename)
      
      Internal, inits an Atoms object from an Abinit NetCDF output.
      
      Note that we are *expecting* something like a _HIST file with limited
      variables. At minimum, we need to have xcart, natom and rprimd in the
      output. The _OUT.nc files do not have xcart at present (Abinit 6.12.1).
      
      Also note that since the _HIST file doesn't contain typat/znucl variables,
      we just set everything to carbon by default.
      
      """
      
      self.clean()
      
      f = netcdf.netcdf_file(filename, 'r')

      self.nsteps = len(f.variables['mdtime'].data)
      if DEBUG: print self.nsteps, f.variables['mdtime'].data  
      xcart = f.variables['xcart'].data
      rprimd = f.variables['rprimd'].data
      fcart = f.variables['fcart'].data
      natom = f.dimensions['natom']
      
      for i in range(self.nsteps):
        self.lattice.append(rprimd[i])
        self.positions.append([xc for xc in xcart[i]])
        self.forces.append([fc for fc in fcart[i]])
        self.species.append(natom * [6])
         
    def loadFromElk(self, filename):
        """ atoms= Atoms.loadFromElk(filename)
        
        Internal, inits an Atoms object from an elk input file.
        
        """
        
        self.clean()
        
        f = open(filename, 'r')
        
        # Construct our datastream
        
        lines = f.readlines()
        f.close()
        
        data = remove_comments(lines, "#")
        data = remove_comments(lines, "!")
        data = remove_comments(lines, ":")
        data = " ".join(data).split()
        
        # Read scales and primitive vectors
        s = elk_array(data, "scale", 1)
        if s is None:
          s = 1.0
        s1 = elk_array(data, "scale1", 1)
        if s1 is None:
          s1 = 1.0
        s2 = elk_array(data, "scale2", 1)
        if s2 is None:
          s2 = 1.0
        s3 = elk_array(data, "scale3", 1)
        if s3 is None:
          s3 = 1.0
        avec = elk_array(data, "avec", 9, newshape=(3,3))
        if avec is None:
          avec = array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape((3,3))
        avec = s * avec
        avec[0,:] = s1 * avec[0,:]
        avec[1,:] = s2 * avec[1,:]
        avec[2,:] = s3 * avec[2,:]
        self.lattice = [[array(x) for x in avec.tolist()]]
        
        # Parse the atoms block
        start = data.index("atoms")
        nspecies = int(data[start+1])
        curpos = start+2
        spec = []
        pos = []
        for i in range(nspecies):
          spfname = data[curpos].strip("'").strip('"').split(".")[0]
          ncurspec = int(data[curpos+1])
          spec = spec + ncurspec * [getElementZ(spfname)]
          for j in range(ncurspec):
            pstart = curpos+2 + j*6
            pos.append(array([float(x) for x in data[pstart:pstart+3]]))
          curpos = pstart + 6
        self.species.append(spec)
        
        # Need to check if the molecule flag is set: if not, convert
        # from reduced to cart coords.
        if "molecule" in data:
          if data[data.index("molecule") + 1].lower() == ".true.":
            self.is_crystal = False
          else:
            self.is_crystal = True
        else:
          self.is_crystal = True
        
        if self.is_crystal:
          pos = reduced2cart(pos, avec)
         
        self.positions.append(pos)
    
    def loadFromAbinitWFK(self, wfk_file):
      """ atoms = Atoms.loadFromAbinitWFK(dens_file)
     
      Internal, inits an Atoms object from a _WFK file written
      by abinit.
      
      """
     
      self.clean()
       
      # Libabitools does the work here (fortran library) because
      # we're potentially reading in a gigantic file.
       
      io.wavefunction(wfk_file)
      
      self.lattice = [[array(x) for x in io.rprimd.tolist()]]
      self.positions = [reduced2cart([array(x) for x in io.xred.T.tolist()], self.lattice[0])]
      self.species = [[int(io.znucltypat[x-1]) for x in io.typat]]
      
      
    def loadFromAbinitDensity(self, dens_file):
        """ atoms= Atoms.loadFromAbinitDensity(dens_file)
        
        Internal, inits an Atoms object from a _DEN file written
        by abinit.
        
        """
        
        self.clean()
        
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
        self.positions = [reduced2cart([array(x) for x in io.xred.T.tolist()], self.lattice[0])]
        self.species = [[int(io.znucltypat[x-1]) for x in io.typat]]
        self.densities = [dens2]
                                            
    def loadFromAbinit(self, abinit_input):
        """ atoms = Atoms.loadFromAbinit(abinit_input)
        
        Internal, inits an Atoms object from an abinit input file. We
        only allow default values for rprim - input files must contain:
        
        natom, ntypat, typat, znucl, acell, [xangst, xred, xcart]
        
        """
        
        # Destroy current data
        self.clean()
        
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
        if type(typat) == type(1):
          typat = [typat]
        
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
        
    def loadFromVASP(self, vasp_file, ions=None):
      """ atoms = Atoms.loadFromVASP(vasp_file, ions=[])
      
      Internal, inits an Atoms object from a VASP POSCAR or CONTCAR file. Very rudimentary
      for the present but good enough for converting from CONTCAR to other formats.
      
      Neither the POSCAR nor CONTCAR files actually specify what species are listed, 
      so you can give a list (strings like 'C' or numbers like 6) of the species to match
      those in the file, using the ions input.
      
      """
      
      # Destroy our current data
      self.clean()
      
      # We must be a crystal, and NO animations here.
      self.is_crystal = True
      self.nsteps = 1
      
      f = open(vasp_file)
      if not f:
        raise ESCError("File %s could not be opened - exiting.")
          
      lines = f.readlines()
      f.close()
      
      # The POSCAR/CONTCAR format is a line-position based format. First line is a 
      # comment, second line is the scale of the unit cell.
      
      scale = float(lines[1].split()[0])
      
      # Next three lines are the lattice vectors
      a = ang2bohr(array([float(x) for x in lines[2].split()[0:3]])) * scale
      b = ang2bohr(array([float(x) for x in lines[3].split()[0:3]])) * scale
      c = ang2bohr(array([float(x) for x in lines[4].split()[0:3]])) * scale
      
      # Now, depending on the version of VASP we start to get in a spot of bother,
      # so we need to look ahead. If there are two lines of numbers, the first line
      # gives the species, the second gives the species *counts*. If there is only
      # one line of numbers, it is just the species counts. Then the Selective Dynamics
      # line is OPTIONAL and then finally we must have either Direct or Cartesian. So
      # look ahead until we get either direct or cartesian (only the first letter matters).
      
      for i,l in enumerate(lines[5:]):
        if l[0] in ['d', 'D', 'c', 'C']:
          if DEBUG:
            print "Found Direct/Cart at line %d." % i
          # First mark whether it is Direct or Cartesian for later.
          if l[0] in ['d', 'D']:
            ptype = 'Direct'
          else:
            ptype = 'Cartesian'
          if i+5 == 6:
            # No Selective Dynamics, no species.
            # So, lines[5] is the set of ion counts.
            cline = 5
            iline = None
            pline = 7
          elif i+5 == 7:
            # Either Selective Dynamics or species are included (not both).
            pline = 8
            if lines[6][0] in ['S', 's']:
              # Selective Dynamics, so species NOT included.
              cline = 5
              iline = None
            else:
              # Species at line 5, counts at line 6
              iline = 5
              cline = 6
          elif i+5 == 8:
            pline = 9
            # Both selective Dynamics and species.
            iline = 5
            cline = 6
          break
      
      # The ions input overrides the iline if it is present. 
      ion_counts = [int(x) for x in lines[cline].split()]
      if ions is not None:
        species = []
        for i, ic in enumerate(ion_counts):
          species += ic * [getElementZ(ions[i])]
      elif iline is not None:
        ion_species = [int(x) for x in lines[iline].split()]
        for i, ic in enumerate(ion_counts):
          species += ic * [getElementZ(ion_species[i])]      
      else:
        species = []
        # If there are no ions specified, just start at H and work up the integers.
        for i, ic in enumerate(ion_counts):
          species += ic * [i+1]
      
      positions = []
      for i in range(sum(ion_counts)):
        pos = array([float(x) for x in lines[pline+i].split()[0:3]])
        if ptype == 'Direct':
          pos = pos[0] * a + pos[1] * b + pos[2] * c
        else:
          pos = ang2bohr(pos)
        positions.append(pos)
      
      self.species.append(species)
      self.positions.append(array(positions))
      self.lattice.append([a, b, c])      
          
                   
    def loadFromXSF(self, xsf_file):
        """ atoms = Atoms.loadFromXSF(xsf_file)
        
        Internal, inits an Atoms object from an xsf file. Note that we can deal with
        animated (AXSF) files just fine.
        
        """
        
        # Destroy our current data
        self.clean()
        
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
                
    def writeWFCube(self, filename, spin,kpt,band,spinor, option="density"):
      """ succeeded = Atoms.writeWaveFunctionCube(filename, spin, kpt, band, spinor, option="density")
      
      If this Atoms instance has a filehook (must point to an open and valid ETSF-formatted
      NetCDF file object), read a wavefunction and write to a CUBE file. 
      
      Note we can't really deal with complex numbers in the CUBE format so have to
      choose an option - density, mod, real or imag.
      
      """
      
      wf = self.getRealSpaceWF(spin,kpt,band,spinor)
      
      if option == "real":
        d = real(wf)
      elif option == "imag":
        d = imag(wf)
      elif option == "density":
        d = real(wf) ** 2 + imag(wf) ** 2
      elif option == "mod":
        d = sqrt(real(wf) ** 2 + imag(wf) ** 2)
        
      return write_cube(filename, self.positions, self.species, self.lattice, d)
                        
    def getRealSpaceWF(self, spin, kpt, band, spinor):
      """ wf = Atoms.getRealSpaceWaveFunction(spin, kpt, band, spinor)
      
      If this Atoms instance has a filehook (pointing to an open ETSF-formatted NetCDF file
      object), use it to read the realspace wavefunction for a single spin, kpt etc.
      
      
      """
      
      wf_real = self.filehook.variables['real_space_wavefunctions'][spin,kpt,band,spinor,:,:,:,0]
      wf_imag = self.filehook.variables['real_space_wavefunctions'][spin,kpt,band,spinor,:,:,:,1]
      
      return wf_real + 1j*wf_imag     
          
        
    def orderAtoms(self, order=None):
      """ new_order = Atoms.orderAtoms(order=None)
      
      Takes the positions and species and reorders them (for example, if the
      species alternate C, H, C, C, N, H, O, etc, can reorder to C C H H N O) 
      according to the passed list order. If order is None, the ordering is done
      by calling uniqify(self.species[0]), the order of the elements in the
      first timestep.
      
      """

      if order is None:
        my_order = uniqify(self.species[0])
      else:
        my_order = uniqify(order)
              
      for t in range(len(self.species)):
        newpos = []
        newspec = []
        pos = self.positions[t]
        spec = self.species[t]

        for o in my_order:
          for i in range(len(spec)):
            if getElementZ(spec[i]) == getElementZ(o):
              newspec.append(spec[i])
              newpos.append(pos[i])
        
        self.positions[t] = newpos
        self.species[t] = newspec
        
      return my_order
      
    def listSpecies(self,spec):
      """ spec_idx = Atoms.printSpecies(spec)
      
      Returns a list of indices of all atoms of type spec. Spec can be given
      as a Z number or an abbreviation (H, Mg, Ca, etc). Note these are zero-
      based indices that can be used within the Atoms object itself, and not
      the 1-based atom indices used by Abinit, for example.
      
      """
      
      spec_idx = []
      
      for i in range(len(self.species[0])):
        if self.species[0][i] == getElementZ(spec):
          spec_idx.append(i)
          
      return spec_idx
        
    
    def generateCASTEPConstraints(self, constrained_atoms):
      """ text_constraints = Atoms.generateCASTEPConstraints(constrained_atoms)
      
      Generates the %BLOCK IONIC_CONSTRAINTS text to fix the listed atoms and 
      returns the text. That's all. CASTEP-style constraints. Note that the
      timestep is irrelevant here so we always take from timestep 0.
      
      The list of constrained atoms is 1-based, not 0-based, to match with
      the style used in CASTEP.
      
      """
      
      text_constraints = ["%block ionic_constraints"]
      count = 0
      
      for i in constrained_atoms:
        spec = elements[self.species[0][i-1]] # i-1, not i, because 1-based.
        # To find the ion number, we count from the first occurrence of the same
        # species.
        ion = i - self.species[0].index(self.species[0][i-1])
        text_constraints.append("  %d  %s  %d 1.0 0.0 0.0" % (count+1, spec, ion))
        text_constraints.append("  %d  %s  %d 0.0 1.0 0.0" % (count+2, spec, ion))
        text_constraints.append("  %d  %s  %d 0.0 0.0 1.0" % (count+3, spec, ion))
        count += 3
      
      text_constraints.append("%endblock ionic_constraints")
      
      return "\n".join(text_constraints)
      
    def writeXSF(self, filename):
      """ success = Atoms.writeXSF(filename)
      
      Member function wrapper for esc_lib.write_xsf.
      
      """
      
      return write_xsf(filename, self.positions, self.species, self.lattice)
      
    def writeAbinit(self, filename, xtype="ang", opt=None, timestep=0):
      """ success = Atoms.writeAbinit(filename, xtype="ang", opt=None, timestep=0)
      
      Member function wrapper for esc_lib.write_abinit.
      
      """
      
      return write_abinit(filename, self.positions, self.species, xtype, opt, timestep)

    def writeCastep(self, filename, xtype="ang", opt=None, timestep=0):
      """ success = Atoms.writeCastep(filename, xtype="ang", opt=None, timestep=0)
      
      Member function wrapper for esc_lib.write_castep.
      
      """
      
      return write_castep(filename, self.positions, self.species, self.lattice, xtype, opt, timestep)  

    def writeAims(self, filename, xtype="ang", opt=None, timestep=0):
      """ success = Atoms.writeAims(filename, xtype="ang", opt=None, timestep=0)
      
      Member function wrapper for esc_lib.write_aims

      """

      return write_aims(filename, self.positions, self.species, self.lattice, xtype, opt, timestep)
    
    def writePDB(self, filename, opt=None, timestep=0):
      """ Write this Atoms object to a PDB file.
      
      No options are implemented as yet.
      
      """
      
      if self.is_crystal:
        return write_pdb(filename, self.positions, self.species, self.lattice, opt, timestep)
      else:
        return write_pdb(filename, self.positions, self.species, opt=opt, timestep=timestep)
      
      
    def writeXYZ(self, filename):
      """ success = Atoms.writeXYZ(filename)
      
      Simple wrapper to write an XYZ file from the atomic positions.
      
      """
      
      return write_xyz(filename, self.positions, self.species)
      
    def writeMolden(self, filename, opt=None):
      """ success = Atoms.writeMolden(filename, opt=None)
      
      Simple wrapper to write a .molden file from the atomic positions. Options not 
      implemented.
      
      """
      
      return write_molden(filename)


            
