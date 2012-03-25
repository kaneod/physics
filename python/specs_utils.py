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

DEBUG = False

def line_index(lines, searchstr):
  """ line_idx = line_index(lines, searchstr):
  
  Returns the line index of the first occurence of searchstr (or None).
  
  """
  
  for i, l in enumerate(lines):
    if searchstr in l:
      return i
  
  return None

def read_xy(filename):
  """ group = specs_utils.read_xy(filename)
  
  At present, just read and return the xy data from a SPECS xy output. We can
  make a few educated guesses about the data based on the header but we don't
  return parsed header values just yet. 
  
  The returned group is of the form:
  
  group = [region1, region2, region3, ...]
  
  where each region is a list [set1, set2, set3]. There may possibly be only
  one set in each region if there are no extended channels.
  
  """
  
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  
  # Blank lines are completely useless in this file format: remove them now.
  lines = remove_comments(lines, just_blanks=True)
  
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
    else:
      if in_comment:
        comms_end.append(i - 1)
        in_comment = False
        
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
  
  
  # Generate a raw datastream and organize into regions/sets to return.
  raw = remove_comments(lines, "#")
  raw = array([float(x) for x in " ".join(raw).split()])
  
  group = []
  
  for v,c,r in zip(set_lengths, col_count, set_regions):
    curset = raw[0:v*c].reshape((v,c))
    raw = raw[v*c:]
    if len(group) < r+1:
      group.append([curset])
    else:
      group[r].append(curset)
    if DEBUG: print "Adding data of shape %d x %d to group %d" % (v,c,r)
    
  return group
  
def check_xy(filename):
  """ is_xy = specs_utils.check_xy(filename)
  
  Checks that a file is actually a specs-written xy file.
  
  """
  
  f = open(filename, 'r')
  
  if f:
    for line in f:
      # This could theoretically run over the whole file but in practice
      # it should break and leave at the first line.
      if ("Created by:" in line) and ("SpecsLab2" in line):
        f.close()
        return True
  else:
    print "(specs_utils.check_xy) ERROR: File %s could not be opened for reading." % filename
    return False

def header_xy(filename, separate_regions=False):
  """ group_name, region_names, all_comments = specs_utils.header_xy(filename, separate_regions=False)
  
  Returns the group and region labels for a given XY file, plus all lines
  that start with comments for post-processing.
  
  If separate_regions is True, the regions list is separated into individual lists for each group.
  
  """
  
  if check_xy(filename):
    f = open(filename, 'r')
  
  group_names = []
  comments = []
  region_names = []
  cur_regions = []
  for line in f:
    if line.startswith("#"):
      comments.append(line)
    if "Group:" in line:
      group_names.append(line.partition("Group:")[2].strip())
      if separate_regions:
        region_names.append(cur_regions)
        cur_regions = []
    if "Region:" in line:
      cur_regions.append(line.partition("Region:")[2].strip())
  
  # Hack!
  if separate_regions:
    region_names.append(cur_regions)
    region_names = region_names[1:]
  else:
    region_names = cur_regions
  f.close()
  
  return group_names, region_names, comments
     
    
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
  

 
