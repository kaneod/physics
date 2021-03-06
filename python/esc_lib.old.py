################################################################################
#
# esc_lib.py
#
# Library of electronic-structure related code.
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
################################################################################


from __future__ import division
from numpy import array, zeros, sqrt, reshape, mat, pi, matrix, cos, sin, exp, arange, arccos, arctan2, complex, polyfit, poly1d, loadtxt, amin, amax, argmin, argmax
from random import random as rand
from numpy.linalg import norm, inv
from numpy.fft import fftn, ifftn
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.special import sph_harm
from scipy.optimize import leastsq, curve_fit, fmin_slsqp
#from scipy.optimize import minimize
#from libabitools import io, wave, spectra
from scipy.io import netcdf
import os

# Debugging flag - set to 1 to see debug messages.
DEBUG=1

# Element dictionaries

elements = { 1 : "H", 2 : "He", 3 : "Li",  4 : "Be", 5 : "B", 6 : "C", 7 : "N" , 8 : "O", 9 : "F", 29 : "Cu", 14 : "Si" , 13 : "Al", 16: "S", 83 : "Bi"}
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

def gl_smear(x, y, xs, gw=None, lw=None, cutoff=10):
  """ xs, ys = gl_smear(x, y, xs, gw=None, lw=None,cutoff=10)
    
  Using a set of reference values given by x and y, create a spectrum ys
  over the range xs by applying a Gaussian/Lorentzian smearing algorithm:
  
  ys(xs) = sum(x) [y(x) * {gaussian(x-xs, gw) + lorentzian(x-xs, lw)}]
  
  The gw and lw parameters are the width characteristic for the gaussian
  and lorentzian respectively. If gw or lw is None, a dynamic smearing will
  be used:
  
  gw(xs[i]) = 0.25 * (xs[i+1] - xs[i-1])
  
  The cutoff parameter can be used to set a low-pass filter on the y data - if
  y > cutoff, it won't be contributed to the spectrum.
  
  NOTE: This gl_smear is incorrect - needs a weighting factor between the gaussian
  and lorentzian to maintain normalization (amongst other things). 
  
  """
  
  # Constants to speed up function evaluation
  g_pre = 1.0 / sqrt(2 * pi)
  l_pre = 1.0 / pi
  
  # Need some functions here
  def gaussian(p,w):
    return g_pre / w * exp(-0.5 * (p ** 2) /  (w ** 2))
    
  def lorentzian(p,w):
    return l_pre * w / (p ** 2 + w ** 2)
    
  ys = zeros((len(xs)))
  dg = []
  dl = []
  
  # Set smearing width
  if gw is None:
    for i in range(len(xs)):
      if i == 0:
        dg.append(abs(0.5 * (xs[1] - xs[0])))
      elif i == len(xs) - 1:
        dg.append(abs(0.5 * (xs[-1] - xs[-2])))
      else:
        dg.append(abs(0.25 * (xs[i+1] - xs[i-1])))
  else:
    dg = [gw] * len(xs)
    
  if lw is None:
    for i in range(len(xs)):
      if i == 0:
        dl.append(abs(0.5 * (xs[1] - xs[0])))
      elif i == len(xs) - 1:
        dl.append(abs(0.5 * (xs[-1] - xs[-2])))
      else:
        dl.append(abs(0.25 * (xs[i+1] - xs[i-1])))
  else:
    dl = [lw] * len(xs)
    
  # Calculate combined smeared function
  for i, xsi in enumerate(xs):
    for xi,yi in zip(x,y):
      if yi < cutoff:
        ys[i] = ys[i] + yi * (gaussian(xsi - xi, dg[i]) + lorentzian(xsi - xi, dl[i]))
  
  return ys
    
def rotate_positions(positions, filename, fileopt=0):
  """ new_positions = rotate_positions(positions, filename, fileopt=0)
    
    Rotates specified coordinates around a specified axis and returns the new positions.
    Can be used, for example, to rotate portions of a molecule. The rotations come
    from a file in the form:
    
    i1 i2 i3 theta
    ...
    ...
    
    (if fileopt = 0) where i1 is the atomic index to be rotated and the axis is 
    formed as the vector between atoms i2 and i3. Rotation is CLOCKWISE 
    viewed along the axis.
    
    If fileopt is 1, the format is i1, x1, x2, x3, o1,o2, o3, theta, where xi 
    give the actual axis and the oi give the origin of the vector xi for the 
    rotation. The units of the axis x don't matter, vector o must be given
    in atomic units (bohr).
    
    Comments marked with # or ! are ignored in the rotation file.
    
    NOTE: the indices i1, i2, etc are 1-based, to align with XCrysden, not
    0-based. We convert inside this routine.
    
  """
  
  new_positions = positions
  
  f = open(filename)
  lines = f.readlines()
  lines = remove_comments(lines, "#")
  lines = remove_comments(lines, "!")
  
  indices = []
  axes = []
  origins = []
  angles = []
  
  for line in lines:
    bits = line.split()
    indices.append(int(bits[0]) - 1) # remember to subtract 1 from the indices
                                     # since python is 0-based!
    if fileopt == 0:
      a = int(bits[1]) - 1
      b = int(bits[2]) - 1
      axes.append(positions[b] - positions[a])
      origins.append(positions[a])
      angles.append(float(bits[3]))
    elif fileopt == 1:
      axes.append(array([float(x) for x in bits[1:4]]))
      origins.append(array([float(x) for x in bits[4:7]]))
      angles.append(float(bits[7]))
    else:
      print "rotate_positions: ERROR - fileopt must be 0 or 1, not %d" % fileopt
  
  for i,a,o,t in zip(indices, axes, origins,angles):
    new_positions[i] = rotate(positions[i], a, o, t)
  
  return new_positions
  
def rotate(pos, axis, origin, angle):
  """ new_vec = rotate(pos, axis, origin, angle)
    
    Rotates a coordinate pos by the specified angle, CLOCKWISE looking along
    axis. Origin of the rotation axis vector is given by the origin parameter,
    the position coordinate pos is assumed to have an origin of (0,0,0) and is 
    converted to a vector with respect to origin before the calculation and 
    restored afterwards.
    
  """
  
  # Find the vector from origin to pos.
  v = pos - origin
    
  # Convert angle to radians
  r = float(angle) * pi / 180
    
  # Make sure axis is normalized
  ux, uy, uz = axis/norm(axis)
    
  # Matrices
  UU = matrix([[ux * ux, ux * uy, ux * uz],
               [ux * uy, uy * uy, uy * uz],
               [ux * uz, uz * uy, uz * uz]])
  
  UX = matrix([[0, -uz, uy],
               [uz, 0, -ux],
               [-uy, ux, 0]])
               
  I = matrix([[1.0, 0.0, 0.0],
              [0.0, 1.0, 0.0],
              [0.0, 0.0, 1.0]])
              
  R = cos(r) * I + sin(r) * UX + (1.0 - cos(r)) * UU
  
  # Apply matrix to column vector of v and restore to a coordinate
  
  return array(R * matrix(v).T).flatten() + origin

def gaussian_convolute(x, y, width):
  """ convolved = gaussian_convolute(x, y, width)
  
  Convolute the data with a gaussian of width w.
  
  """
  
  yc = zeros(y.shape)
  
  for j in range(len(yc)):
    cursum = 0.0
    for i in range(len(x)):
      cursum += y[i] * (1.0 / (width * sqrt(2 * pi))) * exp(-1.0 * (x[j] - x[i]) ** 2 / (2.0 * width ** 2))
    yc[j] = cursum
    
  return yc
  
def normalize_integral(x, y):
  """ yn = normalize_integral(x, y)
  
  Normalizes y to give a unit integral when integrated over x.
  
  """
  
  intsum = 0.0
  for i in range(len(x)-1):
    intsum += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
  
  return y / intsum
  
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
  
def read_xy(filename, comment_delim="#"):
  """ data = read_xy(filename, comment_delim="#")
  
  Reads a multi-column xy file, remove comments, returns as an array
  
  NOTE: This is redundant now that numpy has a loadtxt() method: should replace
  all instances with loadtxt.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  lines = remove_comments(lines, comment_delim)
  f.close()
  
  data = zeros((len(lines), len(lines[0].split())))
  
  for i, line in enumerate(lines):
    for j,bit in enumerate(line.split()):
      data[i,j] = float(bit)
      
  return data
  
def sum_atoms(data, start, finish, option="LeaveFirstColumn"):
  """ summed_data = sum_atoms(data, start, finish)
  
  Assuming an array of the format data[i,j,k], sum over the range in i
  and return an array in the other two indices. 
  
  Option "LeaveFirstColumn" prevents summing over the first column and simply
  copies the first set of values (for i = start). 
  
  """
  
  s = zeros((data.shape[1],data.shape[2]))
  
  if option == "LeaveFirstColumn":
    cs = 1
    s[:,0] = data[start,:,0]
  else:
    cs = 0
    
  for i in range(start,finish):
    s[:,cs:] = s[:,cs:] + data[i,:,cs:]
    
  return s
  
  
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
  
def read_ACF(filename):
  """ charges, total_charge = read_ACF(filename)
  
  Reads the ACF.dat bader output and returns the charge on each atom
  along with the total charge.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  total_charge = float(lines[-1].split()[3])
  
  charges = []
  for line in lines[2:-4]:
    charges.append(float(line.split()[4]))
    
  return array(charges), total_charge

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
  
  return bands, props

def elk_parse_bands(filename):
  """ path, bands = elk_parse_bands(filename="BAND.OUT")
  
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
  
def elk_write_bands(outfile="elk-bands.xy", infile="BAND.OUT"):
  """ result = elk_write_bands(outfile="elk-bands.xy", infile="BAND.OUT")
  
  Reads in the BAND.OUT file (another filename can optionally be specified
  in infile) and writes to a multi-column, tab-delimited xy file for plotting.
  
  """
  
  data = elk_parse_bands(infile)
  
  f = open(outfile, 'w')
  
  for i in range(data.shape[0]):
    f.write("\t".join([str(x) for x in data[i,:]]) + "\n")
    
  f.close()
  
  return True
  
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
          if letter_spec:
            f.write("%s    %g    %g    %g\n" % (elements[species[i][j]], apos[i][j][0], apos[i][j][1], apos[i][j][2]))
          else:
            f.write("%d    %g    %g    %g\n" % (species[i][j], apos[i][j][0], apos[i][j][1], apos[i][j][2]))
        
    f.close()
    return True

def write_aims(filename, positions, species, lattice, xtype="ang", opt=None, timestep=0):
  """ succeeded = write_aims(filename, positions, species, lattice, opt=None, timestep=0)

    Writes a FHI-aims geometry.in file using the given positions, species and lattice. 

    FHI-aims allows either periodic or non-periodic boundary conditions. For periodic,
    specify xtype="frac" and provide lattice vectors. For non-periodic, specify xtype=
    "ang". No lattice vectors will be written in that case.

  """

  pos = positions[timestep]
  spec = species[timestep]

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
    for s, p in zip(spec, pos):
      f.write("atom  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], elements[s]))
  elif xtype == "frac":
    for l in avec:
      f.write("lattice_vector %4.8g %4.8g %4.8g\n" % (l[0], l[1], l[2]))
    
    f.write("\n")
    for s, p in zip(spec, pos):
      f.write("atom_frac  %4.8g %4.8g %4.8g %s\n" % (p[0], p[1], p[2], elements[s]))
  
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

def integrate_grids(a, b):
  """ result =  integrate_grids(a, b)
  
  Performs a triple integration of a, b (complex) 
  over their own grids.
  
  Note: this is a LOT slower than the Fortran version in
  libdft.
  
  """
  
  # Check the shapes are the same
  if (a.shape != b.shape):
    print "(integrate_grids) ERROR: a and b don't have the same shape!"
    return None
  
  # Construct integrand a.conj() * b
  
  igrnd = a.conj() * b
  
  isum = 0.0+0j
  nx = a.shape[0]
  ny = a.shape[1]
  nz = a.shape[2]
  
  for k in range(nz):
    for m in range(ny):
      for i in range(nx):
        isum += igrnd[i,m,k]
  
  return isum/(nx*ny*nz)
  
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
        
class Spectra:
  """ Spectra class: designed to deal with common manipulations on
  spectra output by electronic structure codes like Abinit.
  
  Print Spectra.__init__.__doc__ to see the various ways you can
  construct a Spectra object.
  
  """
  
  def __init__(self, seed, atoms,source):
    """ spectra = Spectra(seed, atoms, source)
    
    Create a Spectra object. For conducti output, use source="conducti". For 
    nexspec output, use source="nexspec". Should also add conventional CASTEP
    ELNES task output here (will be source="castep"). 
    
    For nexspec output, set atoms to a number indicating the species you want to
    read in. For conducti output, atoms should be a range (indicating the 1-based atom
    indices from the original abinit input file.
    
    If the nexspec output is a series of separate runs with files of the form:
    
    SEED_X_S_1_1_1.nexafs
    
    where X is an atom identifier (usually means the location of the core hole)
    and S is the species (the 1s could foreseeably vary but typically won't), you
    can specify source="nexspec-corehole" to parse the atomic index correctly. 
    
    """
    
    if source == "conducti":
      opt = {'seed name' : seed, 'atoms' : atoms}
      self.loadConductiPawCore(opt)
    elif source == "nexspec":
      self.loadNexspec(seed, atoms)
    elif source == "nexspec-corehole":
      self.loadNexspec(seed, atoms, option="corehole")
      
  def loadNexspec(self, seed, species, option=""):
    """ Spectra.loadNexspec(seed, species, option="")
    
    Internal: inits a Spectra object from nexspec output based on the CASTEP ELNES
    task.
    
    NOTE: We can't yet deal with the situation where there are multiple core level
    spectra per atom (ie, 1s, 2s, 2p levels and so on). The current algorithm only
    stores the highest nlm-combination or whatever comes last as listed by listdir.
    
    The option input can presently be:
    
    option="" (default) - file name format is SEED_S_X_B_C.nexafs, read the atom
                          index from X and match the species to S.
    
    option="corehole" - file name format is SEED_X_S_A_B_C.nexafs, read the 
                        atom index from X and match the species to S.
    
    """
    
    self.data = {}
    self.atoms=[]
    
    # Cycle over the contents of the directory containing the nexspec .nexafs outputs
    # and pick out/load the files starting with seed_atoms
    head,tail = os.path.split(seed)
    if head is "":
      # We are in the current directory
      files = os.listdir(os.getcwd())
    else:
      files = os.listdir(head)
    
    for filename in files:
      if filename.startswith(seed) and filename.endswith(".nexafs") and int(filename.split("_")[-4]) == species:
        if option == "corehole":
          self.atoms.append(int(filename.split("_")[1]))
          self.data[int(filename.split("_")[1])] = loadtxt(filename)
        elif option == "":
          self.atoms.append(int(filename.split("_")[2]))
          self.data[int(filename.split("_")[2])] = loadtxt(filename)
    
    self.cmpts = zeros((len(self.atoms), self.data[self.atoms[0]].shape[0],7))
    
    for i,a in enumerate(self.atoms):
      # Now: we have to stick to our rule of always storing and dealing in Hartree
      # atomic units. The nexspec output is in eV, so convert.
      self.data[a][:,0] = eV2hartree(self.data[a][:,0])
      self.cmpts[i,:,:] = self.data[a][:,0:7]
    
  def loadConductiPawCore(self, opt):
    """ Spectra.loadConductiPawCore(opt)
    
    Internal: inits a Spectra object from Conducti PAW core-level outputs.
    
    Note: We can't deal with multiple core-levels in the same file (yet). These will
    be ignored in the read.
    
    """
    
    self.data = {}
    self.atoms = opt['atoms']
    self.spectra_type = "conducti_paw_core"
      
    for i in self.atoms:
      filename = opt['seed name']+str(i)
      self.data[i] = loadtxt(filename)
      
    # Construct a cmpts array (natoms,npoints,2) from data.
    
    self.cmpts = zeros((len(self.atoms), self.data[self.atoms[0]].shape[0],7))
    
    for i,a in enumerate(self.atoms):
      self.cmpts[i,:,:] = self.data[a][:,0:7]
      
  def spectrumAXYZ(self, atom, evec):
    """ spectrum = Spectra.spectrumAXYZ(atom, evec)
    
    Generates the spectrum for a given electric field direction
    and atom.
    
    The evec vector is normalized before use. The spectrum is returned as a
    two-column array with the second column the spectrum and the first column
    the independent variable.
    
    """
    
    # Construct spectrum from components
    cmpts = self.data[atom]
                    
    return spectra.spectrum_axyz(cmpts, evec)
    
  def spectrumXYZ(self, evec):
    """ spectrum = Spectra.spectrumXYZ(evec)
    
    Generates the spectrum across all atoms for a given e-field direction.
    
    See spectrumAXYZ for more details.
    
    """
    
    return spectra.spectrum_xyz(self.cmpts, evec)
    
  def spectrumTP(self, theta, phi):
    """ spectrum = Spectra.spectrumTP(theta,phi)
    
    Generates the spectrum across all atoms for a given e-field direction
    expressed as (theta, phi) in spherical coordinates with theta the polar
    angle with respect to the original z axis.
    
    """
    
    return spectra.spectrum_tp(self.cmpts, theta, phi)
    
  def spectrumComponents(self):
    """ components = Spectra.spectrumComponents()
    
    Returns an array of the same shape as any of the Spectra.data[i]
    but with the 1:6 components (inclusive) summed over all sets i.
    
    """
    
    components = zeros(self.data[self.atoms[0]].shape)
    components[:,0] = self.data[self.atoms[0]][:,0]
    
    for a in self.atoms:
      for i in [1,2,3,4,5,6]:
        components[:,i] = components[:,i] + self.data[a][:,i]
    
    # Note we ignore the 7th column from conducti, it's just an arbitary combination
    # of components set from conducti.
    
    components = components[:,0:7] # Cut off the last column
    return components
    
  def writeSpectrumComponents(self,filename):
    """ succeeded = Spectra.writeSpectrumComponents(filename)
    
    Writes out the output of Spectra.spectrumComponents() to disk
    with the given filename.
    
    """
    
    f = open(filename, 'w')
    cmpts = self.spectrumComponents()
    
    f.write("# Spectrum components from conducti and esc_lib.py\n")
    f.write("#\n")
    f.write("# E XX YY ZZ XY XZ YZ\n")
    
    for c in cmpts:
      f.write("\t".join([str(x) for x in c])+"\n")
      
    f.close()
    
  
  def spectrumATP(self, atom, theta, phi):
    """ spectrum = Spectra.spectrumATP(atom, theta, phi)
    
    Generates a spectrum for *EFIELD* angles theta and phi with respect
    to the coordinate system of the original cell used to compute the 
    spectrum. Note that these angles are NOT those used in experiments as the 
    incident light is at right angles to the e-field and the polarization
    remains a free angle.
    
    Here theta is the polar angle (with respect to z axis) and phi is the
    azimuthal angle in the xy plane.
    
    We assume theta and phi are in degrees and do the conversion.
    
    """
    
    cmpts = self.data[atom]
                    
    return spectra.spectrum_atp(cmpts, theta,phi)
    
  def optimizePeak(self,peak, maximize=True):
    """ fit_result = Spectra.optimizePeak(peak)
    
    Finds the theta,phi that maximizes the given peak. Alternatively if
    maximize=False, we minimize instead. We return the output of the
    scipy.optimize minimize method.
    
    """
    
    i = spectra.closest_index_to_energy(self.cmpts[0,:,0],peak)
    if maximize:
      s = -1.0
    else:
      s = 1.0
      
    def optfunc(p):
      return s * spectra.spectrum_tp(self.cmpts,p[0],p[1])[i,1]
      
    return minimize(optfunc,[rand()*180,rand()*360])
  
  def spectrumRandom(self, samples=1000):
    """ spectrum = Spectra.spectrumRandom(samples=1000)
    
    Uses random numbers to generate angles for theta and phi, and sums
    over the spectra generated for each random angle. Optionally set
    samples to decide how many random angle sets are used.
    
    """
    
    spectrum = zeros((self.cmpts.shape[1], 2))
    
    for i in range(samples):
      spectrum[:,1] = spectrum[:,1] + self.spectrumTP(rand()*180, rand()*360)[:,1]
    
    spectrum[:,0] = self.cmpts[0,:,0]
    spectrum[:,1] = spectrum[:,1] / samples
    return spectrum
    
  def fitXYZ(self, exp_data, energy_range=None, energy_offset=0.0):
    """ a, b, c, I0, data = Spectra.fitXYZ(exp_data, energy_range=None, energy_offset=0.0)
    
    Takes an experimental spectrum (Nx2-shaped array, first column is energy,
    second column is intensity) and finds the best linear combination of 
    orthogonal spectra I0 * (aX+bY+cZ) to fit. Note that this is physically
    incorrect but seems to be what people do in the papers...
    
    See the docstring for Spectra.bestFitToExperiment to find out more about
    the parameters. The only difference is that here we do internal unit
    conversions, assuming the exp_data is in eV.
    
    """
    
    if energy_range is None:
      estart = 0
      eend = self.cmpts.shape[1]
    else:
      estart = spectra.closest_index_to_energy(self.cmpts[0,:,0], energy_range[0])
      eend = spectra.closest_index_to_energy(self.cmpts[0,:,0], energy_range[1])      

    edat = exp_data.copy()   # So we don't mess with the input variable
    edat[:,0] = eV2hartree(edat[:,0] - energy_offset)
    print "Exp_data bounds: ", amin(edat[:,0]), amax(edat[:,0])
    
    y = interp1d(edat[:,0], edat[:,1], bounds_error=False, fill_value=0.0)
    
    # Get our three orthogonal spectra
    X = spectra.spectrum_xyz(self.cmpts, [1,0,0])
    Y = spectra.spectrum_xyz(self.cmpts, [0,1,0])
    Z = spectra.spectrum_xyz(self.cmpts, [0,0,1])
    
    yexp = array([y(x) for x in self.cmpts[0,estart:eend+1,0]])
    def fitfunc(p):
      # Function is p[0] * (p[1] * X + p[2] * Y + p[3] * Z)
      return norm(p[0] * (p[1] * X[estart:eend+1,1] + p[2] * Y[estart:eend+1,1] + p[3] * Z[estart:eend+1,1])-yexp)
      
    r = fmin_slsqp(fitfunc, [1000.0, 0.3, 0.3, 0.3], bounds=[[0.1,1e8],[0.0,1.0],[0.0,1.0],[0.0,1.0]], full_output=True)
    cdat = spectra.spectrum_xyz(self.cmpts, [r[0][1],r[0][2],r[0][3]])
    data = zeros((cdat.shape[0], 7))
    data[:,0] = cdat[:,0]
    data[:,1] = array([y(x) for x in self.cmpts[0,:,0]])
    data[:,2] = r[0][0] * (r[0][1] * X[:,1] + r[0][2] * Y[:,1] + r[0][3] * Z[:,1])
    data[:,3] = r[0][0] * r[0][1] * X[:,1]
    data[:,4] = r[0][0] * r[0][2] * Y[:,1]
    data[:,5] = r[0][0] * r[0][3] * Z[:,1]
    data[:,6] = r[0][0] * cdat[:,1]
    return r[0][1], r[0][2], r[0][3], r[0][0], data
    
    
    
  def bestFitToExperiment(self, exp_data, energy_range=None, energy_offset=0.0):
    """ theta, phi, I0, data = Spectra.bestFitToExperiment(exp_data, energy_range=None, energy_offset=0.0)
    
    Takes an experimental spectrum (Nx2-shaped array, first column is energy,
    second column is intensity) and perform a least squares fit over 
    the specified energy range. If the energy range is omitted (or None), 
    the whole range of the computed data is used. Note that the experimental
    data is interpolated and set to zero outside the experimental energy range.
    
    The energy_offset parameter can be used to impose a rigid x-axis shift of
    the computed values before fitting. If None, the optimal x-axis shift is
    computed as a fitting parameter. By default the shift is 0.0, not None.
    
    Note that we do *not* do unit conversions here, these must be accounted
    for before entry into the routine. The data output is a Nx3 array, first 
    column is the energy scale, second is the experimental and third is the 
    fitted computed data.
    
    """
    
    if energy_range is None:
      estart = 0
      eend = self.cmpts.shape[1]
    else:
      estart = spectra.closest_index_to_energy(self.cmpts[0,:,0], energy_range[0])
      eend = spectra.closest_index_to_energy(self.cmpts[0,:,0], energy_range[1])
    
    if DEBUG:
      print "Fitting computational data between indices %d and %d." % (estart, eend)
      print "Exp_data range is", amin(exp_data[:,0]), amax(exp_data[:,0])
    y = interp1d(exp_data[:,0], exp_data[:,1], bounds_error=False, fill_value=0.0)
    ran = self.spectrumRandom()
    
    yexp = array([y(x-energy_offset) for x in self.cmpts[0,estart:eend+1,0]])
    
    def fitfunc(p):
      # p[0] = theta, p[1] = phi, p[2] = I0, p[3] = mixing fraction
      spec = spectra.spectrum_tp(self.cmpts, p[0],p[1])
      return norm(p[2] * (p[3] * spec[estart:eend+1,1] + (1-p[3]) * ran[estart:eend+1,1]) - yexp)

    r = fmin_slsqp(fitfunc, [45.0, 45.0, 100, 0.5], bounds=[[0.0,180.0],[0.0,360.0], [0.1, 1e5], [0.0, 1.0]], full_output=True, iter=1000, iprint=2, acc=1.0e-9, epsilon=0.001)
    cdat = spectra.spectrum_tp(self.cmpts, r[0][0], r[0][1])
    cdat[:,0] = cdat[:,0] - energy_offset
    cdat[:,1] = r[0][2] * cdat[:,1]
    rdat = r[0][2] * self.spectrumRandom()[:,1]
    # Construct an output array with columns energy, experimental intensity,
    # calculated best fit intensity, random component, oriented component.
    data = zeros((cdat.shape[0], 5))
    data[:,0] = cdat[:,0]
    data[:,1] = array([y(x-energy_offset) for x in self.cmpts[0,:,0]])
    data[:,2] = r[0][3] * cdat[:,1] + (1.0 - r[0][3]) * rdat 
    data[:,3] = (1.0 - r[0][3]) * rdat
    data[:,4] = r[0][3] * cdat[:,1]
    return r[0][0], r[0][1],r[0][2],r[0][3], data
       
  def bestFitToFile(self, filename, energy_range=None,exp_offset=0.0, comp_offset=0.0):
    """ theta, phi, I0, data = Spectra.bestFitToFile(filename, energy_range=None, energy_offset=0.0)
    
    Same as bestFitToExperiment but loads automatically from a file. The
    pretty safe assumption here is that the experimental energy axis is in eV
    whereas the stored spectra here are all in Hartree, so the spectra are
    returned in Hartree. There are two offsets allowed, with exp_offset applied
    to the X-axis of the experimental data and assumed to be in eV, and
    comp_offset passed directly to bestFitToExperiment.
    
    """
    
    # Load experimental data, offset and convert to Ha.
    expdat = loadtxt(filename)
    expdat[:,0] = eV2hartree(expdat[:,0] - exp_offset)
    
    return self.bestFitToExperiment(expdat, energy_range, comp_offset)
    
class Atom:
  """ Atom class: a single atom, as represented by, for example, a 
  PAW data set.
  
  Create using Atom("my_atom.paw"), where the input file is the
  abinit-formatted pseudopotential file. 
  
  """
  
  def __init__(self, pawfile):
    """ atom = Atom(pawfile)
    
    Creates an atom object by reading a abinit PAW dataset.
    
    """
    
    # In the interests of future-proofing, offload to an internal.
    self.loadFromAbinit(pawfile)
  
  def loadFromAbinit(self, pawfile):
    """ Atom.loadFromAbinit(pawfile)
    
    Internal: initializes atom data from abinit PAW file.
    
    """
    
    f = open(pawfile)
    lines = f.readlines()
    f.close()
    
    # Grab the important bits from the header
    
    if "All-electron core wavefunctions" in lines[0]:
    
      self.element = lines[0].split()[9]
      self.generation = lines[0].split("-")[3].strip()
      
      self.method = int(lines[1].split()[0])
      self.nspinor = int(lines[1].split()[1])
      self.nsppol = int(lines[1].split()[2])
      
      self.zatom = float(lines[2].split()[0])
      self.zion = float(lines[2].split()[1])
      self.pspdat = int(lines[2].split()[2])
      
      self.pspcod = int(lines[3].split()[0])
      self.pspxc = int(lines[3].split()[1])
      self.lmax = int(lines[3].split()[2])
      
      self.pspfmt = lines[4].split()[0]
      self.creator = int(lines[4].split()[1])
      
      self.basis_size = int(lines[5].split()[0])
      self.lmn_size = int(lines[5].split()[1])

      self.orbitals = [int(x) for x in lines[6].split()[0:self.basis_size]]
      self.number_of_meshes = int(lines[7].split()[0])
      
      self.mesh_info = []
      self.mesh = []
      for i in range(self.number_of_meshes):
        meshbits = lines[8+i].split()
        tmpdict = {}
        tmpdict["type"] = int(meshbits[1])
        tmpdict["size"] = int(meshbits[2])
        tmpdict["rad_step"] = float(meshbits[3])
        tmpdict["log_step"] = float(meshbits[4])
        self.mesh_info.append(tmpdict)
        j = arange(1,tmpdict["size"]+1)
        if tmpdict["type"] == 1:
          # Linear grid?
          self.mesh.append(tmpdict["rad_step"] * j)
        elif tmpdict["type"] == 2:
          self.mesh.append(tmpdict["rad_step"] * (exp(tmpdict["log_step"] * (j - 1)) - 1)) 
        
        self.r_max = float(lines[8+self.number_of_meshes].split()[0])
        data = lines[9+self.number_of_meshes:]
                
    else:
      self.element = lines[0].split()[5]
      self.generation = lines[0].split("-")[1].strip()
    
      self.zatom = float(lines[1].split()[0])
      self.zion = float(lines[1].split()[1])
      self.pspdat = int(lines[1].split()[2])
    
      self.pspcod = int(lines[2].split()[0])
      self.pspxc = int(lines[2].split()[1])
      self.lmax = int(lines[2].split()[2])
      self.lloc = int(lines[2].split()[3])
      self.mmax = int(lines[2].split()[4])
      self.r2well = float(lines[2].split()[5])
    
      self.pspfmt = lines[3].split()[0]
      self.creator = int(lines[3].split()[1])
    
      self.basis_size = int(lines[4].split()[0])
      self.lmn_size = int(lines[4].split()[1])
    
      self.orbitals = [int(x) for x in lines[5].split()[0:self.basis_size]]
      self.number_of_meshes = int(lines[6].split()[0])
    
      self.mesh_info = []
      self.mesh = []
      for i in range(self.number_of_meshes):
        meshbits = lines[7+i].split()
        tmpdict = {}
        tmpdict["type"] = int(meshbits[1])
        tmpdict["size"] = int(meshbits[2])
        tmpdict["rad_step"] = float(meshbits[3])
        tmpdict["log_step"] = float(meshbits[4])
        self.mesh_info.append(tmpdict)
        j = arange(1,tmpdict["size"]+1)
        if tmpdict["type"] == 1:
          # Linear grid?
          self.mesh.append(tmpdict["rad_step"] * j)
        elif tmpdict["type"] == 2:
          self.mesh.append(tmpdict["rad_step"] * (exp(tmpdict["log_step"] * (j - 1)) - 1))
      
      self.r_cut = float(lines[7+self.number_of_meshes].split()[0])
      self.shape_type = int(lines[8+self.number_of_meshes].split()[0])
      self.rshape = float(lines[8+self.number_of_meshes].split()[1])
      data = lines[9+self.number_of_meshes:]
    
    self.data = []
    # Now go through and grab all the bits!
    is_data = True
    while (is_data):
      # Header of the chunk specifies the length.
      chunk = {}
      chunk["title"] = data[0].split("=====")[1].strip()
      chunk["comment"] = data[0].split("=====")[2].strip()
      print chunk["title"]
      # There are a few special chunks that have slightly different formatting
      if "Dij0" not in chunk["title"] and "Rhoij0" not in chunk["title"]:
        chunk["mesh index"] = int(data[1].split()[0])
        if "VHntZC" in chunk["title"]:
          try:
            self.vloc_format = int(data[1].split()[1])
          except ValueError:
            # Some PAW files don't specify the vloc_format
            self.vloc_format = 1 
        elif "Core wave functions" in chunk["title"]:
          chunk["core n"] = int(data[2].split()[0])
          chunk["core l"] = int(data[2].split()[1])
          chunk["core s"] = int(data[2].split()[2])
          chunk["core E"] = float(data[3].split()[0]) / 2 # Convert immediately to Ha from Ry.
          chunk["core occ"] = float(data[3].split()[1])
          data = data[2:]
        data = data[2:]
        size = self.mesh_info[chunk["mesh index"] - 1]["size"]
      else:
        size = self.lmn_size * (self.lmn_size + 1) * 0.5
        chunk["mesh index"] = -1
        data = data[1:]
      chunk["data"] = []
      got_data = False
      while not got_data:
        for val in [float(x) for x in data[0].split()]:
          chunk["data"].append(val)
        if len(chunk["data"]) == size:
          data = data[1:]
          got_data = True
        else:
          data = data[1:]
      chunk["data"] = array(chunk["data"])
      self.data.append(chunk)
      if len(data) == 0:
        is_data = False
    
    if "paw" in self.pspfmt:
      # Expand rhoij and dij matrices into a triangular matrix.
      for i, d in enumerate(self.data):
        if "Rhoij0" in d['title']:
          self.iRhoij0 = i
        elif "Dij0" in d['title']:
          self.iDij0 = i
    
      self.Dij0 = zeros((self.lmn_size, self.lmn_size))
      self.Rhoij0 = zeros((self.lmn_size, self.lmn_size))
    
      for row in range(1, self.lmn_size+1):
        for column in range(1, row+1):
          self.Dij0[row - 1, column - 1] = self.data[self.iDij0]['data'][(row  - 1 + (row -1 )** 2) / 2 + column - 1]
          self.Rhoij0[row - 1, column - 1] = self.data[self.iRhoij0]['data'][(row - 1 + (row - 1) ** 2) / 2 + column - 1]
          
    # Make interpolators for the radial data so we don't need to construct
    # them every time we need them.
    self.interpolators = {}
    for i, d in enumerate(self.data):
      if d['mesh index'] != -1:
        self.interpolators[i] = interp1d(self.mesh[d['mesh index'] - 1], d['data'])
    
  def radial2Cart(self, data_index, r, r0):
    """ psi = Atom.radial2Cart(data_index, r, r0)
    
    Given an atomic centre r0, calculate the actual value of the requested function
    in realspace at position r.
    
    r, r0 are expected to be 3-element numpy arrays.
    
    """
    
    # Check we're within the PAW cutoff radius.
    if "paw" in self.pspfmt:
      rcut = self.r_cut
    elif "core" in self.pspfmt:
      rcut = self.r_max
    if norm(r - r0) > rcut:
      return 0
    if self.data[data_index]['mesh index'] == -1:
      print "Requested data does not represent a radial function."
      return None
    else:
      #ur = self.data[data_index]['data']
      #rr = self.mesh[self.data[data_index]['mesh index'] - 1]
      u = self.interpolators[data_index]
      
      # Use the title of the data to figure out which orbital
      # this is (hence the angular momentum)
      
    l = self.orbitals[int(self.data[data_index]['title'].split()[-1])-1]
    
    # Get the spherical coordinates of our displacement vector  
    s = r - r0
    rho = norm(s)
    # Need a different treatment if rho == 0.
    if rho > 0:
      theta = arccos(s[2]/rho)
      phi = arctan2(s[1], s[0])
    
      # We aren't doing this in the presence of a magnetic field
      # so we can just use m=0 for our spherical harmonic.
      return u(rho) / rho * sph_harm(0,l,theta,phi)
      
    else:
      theta = 0.0
      phi = arctan2(s[1],s[0]) # Doesn't really matter what this is
      
      # Use a 3rd order polynomial fit to points close to the origin to extrapolate.
      y = [u(x)/x for x in [0.001, 0.002, 0.003]]
      z = polyfit([0.001, 0.002, 0.003], y, 3)
      return poly1d(z)(0.0) * sph_harm(0,l,theta,phi)
    
  def expandOntoGrid(self, data_index, ngrid, avec, r0):
    """ grid = Atom.expandOntoGrid(, data_index, ngrid, avec, r0)
    
    Expand the specified dataset onto a grid with grid
    dimensions ngrid and unit vectors avec[i]/ngrid[i]. 
    
    Place the atomic sphere at position r0.
    
    """
    
    grid = zeros(ngrid, dtype="complex")
    
    for i in range(ngrid[0]):
      for j in range(ngrid[1]):
        for k in range(ngrid[2]):
          r = (1.0 * i) /ngrid[0] * avec[0] + (1.0 * j)/ngrid[1] * avec[1] + (1.0 * k)/ngrid[2] * avec[2]
          grid[i,j,k] = self.radial2Cart(data_index, r, r0)
    
    # Renormalize.
    #inorm = integrate_grids(grid,grid)
    #return grid/sqrt(inorm)
    return grid         
        
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
        
        "abinit" : Abinit input file.
        
        "abi_density" : Abinit _DEN file.
        
        "elk" : Elk input file.
        
        "castep" : CASTEP .cell file.
        
        "NetCDF" : Abinit NetCDF output. Note: this importer is not very clever
        and will put in default values for the species if they cannot be found -
        default is all atoms are carbon.
        
        "ETSF" : Read from ETSF-formatted NetCDF output.
        
        """
        
        if filetype == "XSF":
            self.loadFromXSF(filename)
        elif filetype == "abinit":
            self.loadFromAbinit(filename)
        elif filetype == "abi_density":
            self.loadFromAbinitDensity(filename)
        elif filetype == "abi_wfk":
            self.loadFromAbinitWFK(filename)
        elif filetype == "elk":
            self.loadFromElk(filename)
        elif filetype == "NetCDF":
            self.loadFromNetCDF(filename)
        elif filetype == "ETSF":
            self.loadFromETSF(filename)
        elif filetype == "castep":
            self.loadFromCastep(filename)
        elif filetype == "VASP":
            self.loadFromVASP(filename, options)
        else:
          print "(esc_lib.Atoms.__init__) ERROR: File type %s not handled at present." % filetype
          return None
    
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
    
#     def gradWF(self, spin, kpt, band, spinor):
#       """ gx, gy, gz = Atoms.gradWF(spin,kpt,band,spinor)
#       
#       If this Atoms instance has a filehook to a ETSF WF file, return
#       the gradient of a specific wavefunction using a FFT method.
#       
#       The returned grids are the components of the complex gradient
#       with respect to the cartesian axes.
#       
#       """
#       
#       psi = self.getRealSpaceWF(spin,kpt,band,spinor)
#       reduced_k = self.filehook.variables['reduced_coordinates_of_kpoints'][kpt]
#       K = reduced2cart(reduced_k, self.recip_lattice[0])
#       
#       print "Got psi and K"
#       psiG = fftn(psi)
#       print "Got past fftn!"
#       
#       speed.grad_wf(K, self.g_vectors, psiG)
#       grad = speed.grd[:]
#       
#       print "Got past grad!"
#       gx = grad[:,:,:,0]
#       gy = grad[:,:,:,1]
#       gz = grad[:,:,:,2]
#       
#       return ifftn(gx), ifftn(gy), ifftn(gz)
                
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
                
    def getBondLengths(self, cutoff=3.0, animstep=0, give_species=False):
        """ bonds = Atoms.getBondLengths(cutoff=3.0, animstep=1, give_species=False)
        
        Returns a list of bond specs [i, j, length] for all pairwise
        distances less than the cutoff distance (default: 3.0 Bohr) for
        the specified animation step (default: first step, ie 0).
        
        If give_species=True, the bond spec includes the species abbreviation:
        
        [i, Zi, j, Zj, length]
        
        """
        
        return getBondLengths(self.positions[animstep], self.species[animstep], cutoff, give_species)
        
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
        
    def autoUnLars(self, animstep=0, convtol=0.01, maxsteps=50, cutoff=2.2):
        """ new_pos, pos_hist = Atoms.autoUnLars(animstep=0, convtol=0.01, maxsteps=50, cutoff=2.2)
        
        Attempts to fix positions generated by Lars automatically. This
        might not work very well. Converges positions to convtol.
        
        Why? Because Dr Lars Thomsen likes the gunslinging approach to chemical
        construction. He laughs in the face of the 2nd Law of Thermodynamics.
        
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
    
    def generateConstraints(self, constrained_atoms):
      """ text_constraints = Atoms.generateConstraints(constrained_atoms)
      
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
      
    def rotateAtoms(rotation_file, fileopt=0, timestep=0):
      """ success = Atoms.rotateAtoms(rotation_file, fileopt=0)
      
      Apply the rotations specified in rotation_file to an Atoms object positions.
      
      We can deal with animation steps here (Default is 0). See the rotate_positions
      documentation for details on fileopt.
      
      """
      
      pos = self.positions[timestep]
      self.positions[timestep] = rotate_positions(pos, rotation_file, fileopt)
      
      return True
      
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



            
