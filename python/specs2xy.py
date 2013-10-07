#!/usr/bin/env python
################################################################################
#
# specs2xy.py
#
# Converts a SPECSLab2 .xml file to .xy of energy and counts for each file (inside
# folders for each group).
#
# Usage: specs2xy.py FILE
#
################################################################################
#
# Copyright 2013 Kane O'Donnell
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
# 1. Uses specs.py.
#
# 2. Does not respect existing folders - will overwrite.
#
# 3. ONLY outputs the total counts, the extended channels and individual 
#    channeltron outputs are ignored.
#
################################################################################

import specs
import argparse
import os
import os.path
import subprocess
import sys
from numpy import *

parser = argparse.ArgumentParser(description="Convert SpecsLAB2 .xml to .xy")
parser.add_argument('file', help="Name of the .xml file to convert.")
args = parser.parse_args()

t = specs.SPECS(args.file)

print "specs2xy.py"
print ""
print "Written by Kane O'Donnell (kane.odonnell@gmail.com)"
print ""
print "Working on file: ", args.file

if t:
  # The file is valid, so strip the ".xml", make a directory. We overwrite if it exists.
  base_path = os.path.split(args.file)[0]
  cwd = os.getcwd()
  unpack_path = os.path.join(base_path, args.file.split(".xml")[0].strip())
  print "Creating unpack directory: ", unpack_path
  try:
    os.mkdir(unpack_path)
  except OSError:
    # Most probable cause is the directory already exists.
    # Right then. Delete that path, and try again.
    try:
      subprocess.check_output(["rm", "-r", unpack_path])
      subprocess.check_output(["mkdir", unpack_path])
    except CalledProcessError:
      # Ok still can't make the directory! I give up.
      print "ERROR: Could not create directory for unpacking. Exiting..."
      sys.exit(1)
      
  os.chdir(unpack_path)
  unpack_path = os.getcwd()
  
  # Now loop over groups. Need to make unique folders for each group.
  for g in t.groups:
    index = 0
    test_name = g.name.strip()
    while os.path.exists(test_name):
      index += 1
      test_name = g.name.strip()+"-"+str(index)
    os.mkdir(test_name)
    os.chdir(test_name)
    
    print "Working in group ", g.name.strip(), " - unpacking into folder ", test_name
    
    # Loop over regions within this group.
    for r in g.regions:
      # Again, make sure we have a unique name for the region to make a file.
      index = 0
      test_name = r.name.strip()+".xy"
      while os.path.exists(test_name):
        index += 1 
        test_name = r.name.strip()+"_"+str(index)+".xy"
      
      print "Region: ", r.name.strip(), " - writing as ", test_name
      
      # Make a header consistent with sinspect
      h = '#"'
      h += 'Analyzer mode:{}'.format(r.scan_mode)
      h += ', Dwell time:{}'.format(r.dwell_time)
      h += ', Pass energy:{}'.format(r.pass_energy)
      h += ', Lens mode:{}'.format(r.analyzer_lens)
      if r.scan_mode=='FixedAnalyzerTransmission':
        h += ', Excitation energy:{}'.format(r.excitation_energy)
      elif r.scan_mode=='ConstantFinalState':
        h += ', Kinetic energy:{}'.format(r.kinetic_energy)
      h += '"\n#"'
      
      x = []
      # Decide which x-axis we need.
      if r.scan_mode =="FixedAnalyzerTransmission":
        x.append(r.binding_axis)
        h += 'Binding Energy (eV)    Counts"'
      elif r.scan_mode == "ConstantFinalState":
        x.append(r.excitation_axis)
        h += 'Excitation Energy (eV)    Counts"'
      elif r.scan_mode == "FixedEnergies":
        x.append(r.time_axis)
        h += 'Time (s)    Counts"'
      else:
        x.append(r.kinetic_axis)
        h += 'Kinetic Energy (eV)    Counts"'
      
      x.append(r.counts)
      x = array(x).transpose()
      with open(test_name, 'w') as f:
        print >> f, h
        savetxt(f, x, fmt='%1.8g')
        
    
    # Make sure we return to the unpack directory before we process the next group.
    os.chdir(unpack_path)
else:
  print "ERROR: ElementTree cannot recognize this as an xml file - are you sure it's SPECSLab .xml output?"
  sys.exit(1)


      