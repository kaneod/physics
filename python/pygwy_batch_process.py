################################################################################
#
#
# Kane's Great Gwyddion Batch Processor!
#
#
################################################################################
#
# Edit the "Options" section below to get going.
#
# The gwyddion pygwy interface is incomplete, which means they've left out
# several functions from the underlying C libraries like gwymath and so on.
# The only one needed here is the function for median line correction, but if
# further features are required it's quite easy to duck into the Gwyddion source
# and figure out how they do something (go GPL!) so let me know if you need it.
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

import gwy
import gwyutils
from os.path import join

################################################################################
#
# Options
#
################################################################################

# For each entry in the file called file_list, we do the processing.
file_list = "/path/to/my/files.txt"

# Put the root path to the folder containing all the files here.
path_root = "/path/to/my"

# Which data set to use in each file.
data_set = 1 # 0 = TraceUp, 1 = RetraceUp, 2 = TraceDown, 3 = RetraceDown

# Colour ranges, min and max.
colour_min = 30.0
colour_max = 90.0

# Output file type. Set this to "png" for sanity, but technically all
# Gwyddion-supported file extensions are possible. What I'm saying here
# is that if you try anything else, you may regret it. On my work system 
# (Ubuntu) setting "tiff" causes Gwyddion to hang. :o)
file_type = "png"

# Colour palette. Default is Gwyddion.net to almost match MATRIX.
palette = "Gwyddion.net"

##############
#
# modules/process/linecorrect.c - line_correct_median
#
##############
def line_correct_median(data):
  """ corrected = line_correct_median(data)
  
  Does a median line correction on the passed gwy.DataField.
  
  """
  
  xres = data.get_xres()
  yres = data.get_yres()
  line = gwy.DataLine(xres,1.0,False)
  modi = gwy.DataLine(yres,1.0,False)
  
  for i in range(yres):
    data.get_row(line,i)
    median = line.get_median()
    modi.set_val(i,median)
  
  median = modi.get_median()
  
  for i in range(yres):
    data.area_add(0, i, xres, 1, median - modi.get_val(i))
    
  return data
  
################################################################################
#
# Script starts here
#
################################################################################

fl = open(file_list, 'r')
files = " ".join(fl.readlines()).split() # Ha, do that in Matlab!
fl.close()

saved_once = False

for f in files:

  mtrx_file = join(path_root, f)
  out_file = join(path_root, f) + "." + file_type
  
  container = gwy.gwy_app_file_load(mtrx_file)
  datafields = gwyutils.get_data_fields_dir(container)
  
  # Stupid bloody /0/data/nonsense
  key = "/"+str(data_set)+"/data"
  data = datafields[key]
  
  # Plane level
  a, bx, by = data.fit_plane()
  data.plane_level(a, bx, by)
  
  # Median line correction
  data = line_correct_median(data)
  
  # Fix zero
  data.add(-1.0 * data.get_min())
  
  # Set colour range
  base = "/"+str(data_set)+"/base"
  container.set_int32_by_name(base+"/range-type", 1)
  container.set_double_by_name(base+"/min", colour_min)
  container.set_double_by_name(base+"/max", colour_max)
  container.set_string_by_name(base+"/palette", palette)
  
  # Do silly things required to save
  gwy.gwy_app_data_browser_reset_visibility(container, gwy.VISIBILITY_RESET_SHOW_ALL)
  gwy.gwy_app_data_browser_select_data_field(container, data_set)
  if saved_once:
    gwy.gwy_file_save(container, out_file, gwy.RUN_NONINTERACTIVE)
  else:
    gwy.gwy_file_save(container, out_file, gwy.RUN_INTERACTIVE)
    saved_once = True
  gwy.gwy_app_data_browser_reset_visibility(container, gwy.VISIBILITY_RESET_HIDE_ALL)
