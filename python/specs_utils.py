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

def line_index(lines, searchstr):
  """ line_idx = line_index(lines, searchstr):
  
  Returns the line index of the first occurence of searchstr (or None).
  
  """
  
  for i, l in enumerate(lines):
    if searchstr in l:
      return i
  
  return None

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
  
  # Blank lines are completely useless in this file format: remove them now.
  lines = remove_comments(lines, "@")
  
  # Check if we have extended channels.
  if lines[line_index(lines, "External Channel Data:")].split()[4] == "yes":
    have_channels = True
    channels = []
  
  # Locate commentblocks and other important identifiers
  comms_begin = []
  comms_end = []
  regions = []
  clabels = []
  col_count = []
  in_comment = False
  
  for i, line in enumerate(lines):
    if line.startswith("#"):
      if not in_comment:
        comms_begin.append(i)
        in_comment = True
      if "Region:" in line:
        regions.append(i)
      if "ColumnLabels:" in line:
        clabels.append(i)
        col_count.append(len(line.split()) - 2)
    elif len(line.strip()) == 0:
      if DEBUG: print "Blank at %d" % i
    else:
      if in_comment:
        comms_end.append(i - 1)
        in_comment = False
        
  if DEBUG: print comms_begin
  if DEBUG: print comms_end
  if DEBUG: print regions
  if DEBUG: print clabels
  if DEBUG: print col_count
  
  # Decide which sets go in which region.
  set_regions = len(clabels) * [0]
  for i, c in enumerate(clabels):
    for j, r in enumerate(regions):
      if c > r:
        set_regions[i] = j
  
  if DEBUG: print set_regions 

  # Figure out how long each data set is using the comment block sizes.
  set_lengths = len(clabels) * [0]
  set_lengths[-1] = len(lines) - comms_end[-1] - 1
  for i in range(len(set_lengths) - 1):
    set_lengths[i] = comms_begin[i+1] - comms_end[i] - 1
    
  if DEBUG: print set_lengths
  
  
def old_read_xy(filename):
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
  
  if DEBUG: print regions
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
      if DEBUG: print lines[r+i].strip()
      if lines[r+i].strip().startswith("#") or (len(lines[r+i].strip()) == 0):
        if DEBUG: print "Line is a comment."
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
  

 
