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

DEBUG = True

def read_xy(filename):
  """ data = specs_utils.read_xy(filename)
  
  At present, just read and return the xy data from a SPECS xy output. We can
  make a few educated guesses about the data based on the header but we don't
  return parsed header values just yet. 
  
  The returned list, data, is of the form:
  
  data = [set1, set2, set3, ...]
  
  where each set is a 10-element list [main, chan1, chan2, ..., chan9] and the
  external channel bits are empty arrays if external channel data is not
  present. 
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  data = []
  
  # Check whether we have extended channels.
  if lines[8].split()[4] == "yes":
    have_channels=True
    channels = []
  
  # Locate the start of each region using the "Region:" tag.
  # We also look for channels here - there is no hassle separating them
  # because we know there are always only 9 sets if they are present.
  regions = []
  for i,line in enumerate(lines):
    if "Region:" in line.split():
      regions.append(i)
    elif "External Channel Data Cycle:" in line:
      channels.append(i)
  
  # Be tricky with the channels.
  channels = array(channels).reshape((-1,9)).tolist()
  if DEBUG: print channels
  
  if DEBUG: print "Found %d regions." % len(regions)
  
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
        if DEBUG: print "%s values in region starting at line %d." % (line.split()[2], r)
      elif "ColumnLabels" in line:
        dtncols.append(len(line.split()) - 2)
        if DEBUG: print "%d columns for this region." % (len(line.split())-2)  
  
  raw = remove_comments(lines, "#")
  raw = array([float(x) for x in " ".join(raw).split()])
  
  # If we have channels, need to figure out how many columns in each,
  # then interleave with the dtnvals and dtncols lists.
  if have_channels:
    cdtnvals = []
    cdtncols = []
    for v,d,c in zip(dtnvals, dtncols, channels):
      cdtnvals = cdtnvals + 10 * [v]
      cdtncols.append(d)
      for i in [0,1,2,3,4,5,6,7,8]:
        cdtncols.append(len(lines[c[i]+1].split()) - 2)    
    dtnvals = cdtnvals
    dtncols = cdtncols
  
  # Read and store data.
  for v, c in zip(dtnvals, dtncols):
    curset = raw[0:v*c].reshape((v,c))
    raw = raw[v*c:]
    data.append(curset)
  
  if have_channels:
    step = 10  
    sets = [data[i:i+step] for i in range(0, len(data), step)]
  else:
    sets = []
    for d in data:
      sets.append([d] + 9 * [array([])])
    
  return sets
    
def make_interps(filename, columnx=0, columny=1, kind='linear'):
  """ interps, data = specs_utils.make_interps(filename, columnx=0, columny=1, kind='linear')
  
  Generates an interpolator for columny as a function of columnx for each of the
  datasets in the given file. Returns a list of scipy interpolators and the
  original data.
  
  Can specify the kind (linear, nearest, zero, slinear, quadratic, cubic) if
  desired.
  
  """
  
  data = read_xy(filename)
  interps = []
  
  for dset in data:
    curint = []
    for d in dset:
      curint.append(interp1d(d[:,columnx], d[:,columny], kind, bounds_error=False, fill_value=0.0))
    interps.append(curint)
   
  return interps, data
  

 
