################################################################################
#
# specs_utils.py
#
# Read multi datasets exported as xy files from SPECSLab and generate
# interpolaters.
#
################################################################################
#
# Copyright 2012 Kane O'Donnell
#
#    This library is free software: you can redistribute it and/or modify
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
#
#
################################################################################

from numpy import array, zeros
from scipy.interpolate import interp1d
from esc_lib import remove_comments

def read_xy(filename):
  """ data = specs_utils.read_xy(filename)
  
  At present, just read and return the xy data from a SPECS xy output. We can
  make a few educated guesses about the data based on the header but we don't
  return parsed header values just yet. 
  
  The returned list, data, is of the form:
  
  data = [set1, set2, set3, ...]
  
  where each set is a multi-dim array where the rows (first index) are the
  values and the columns denote different outputs. Typically first column is
  energy, second is counts, third is scaling.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  data = []
  
  # Locate the start of each region using the "Region:" tag.
  regions = []
  for i,line in enumerate(lines):
    if "Region:" in line.split():
      regions.append(i)
  
  # Locate the start of the data in each region.
  dtstarts = []
  for r in regions:
    i = 1
    located = False
    while not located:
      if lines[r+i].strip().startswith("#"):
        i = i+1
      else:
        dtstarts.append(r+i)
        located = True
  
  # For each region, read the number of values and infer the number of columns.
  # Note: the columns read is a dirty hack...
  dtnvals = []
  dtncols = []
  for r,d in zip(regions, dtstarts):
    for line in lines[r:d]:
      if "Values/Curve" in line:
        dtnvals.append(int(line.split()[2]))
      elif "ColumnLabels" in line:
        dtncols.append(len(line.split()) - 2)
   
  raw = remove_comments(lines, "#")
  raw = array([float(x) for x in " ".join(raw).split()])
  
  for v, c in zip(dtnvals, dtncols):
    curset = raw[0:v*c].reshape((v,c))
    raw = raw[v*c:]
    data.append(curset)
  
  return data
    
def make_interps(filename, columnx=0, columny=1, kind='linear'):
  """ interps = specs_utils.make_interps(filename, columnx=0, columny=1, kind='linear')
  
  Generates an interpolator for columny as a function of columnx for each of the
  datasets in the given file. Returns a list of scipy interpolators.
  
  Can specify the kind (linear, nearest, zero, slinear, quadratic, cubic) if
  desired.
  
  """
  
  data = read_xy(filename)
  interps = []
  
  for dset in data:
    interps.append(interp1d(dset[:,columnx], dset[:,columny], kind))
   
  return interps
  

 
