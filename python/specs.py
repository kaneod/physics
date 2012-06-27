################################################################################
#
# specs.py
#
# Parse the XML output of SPECSLab v2 into a python object.
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
# 1. I hate XML!
#
# 2. This library is not complete - plan is to incrementally add more member
#    elements to the classes as necessary.
#
################################################################################

from __future__ import division
import xml.etree.ElementTree as ET
from numpy import array, linspace, zeros, ceil
from scipy.interpolate import UnivariateSpline

DEBUG = 1

class SPECS:
  """ Represent a SPECSLab .xml output as a python object. Construct with:
    
      specs_obj = specs.SPECS(my_xml_file)
      
  """
  
  def __init__(self, filename):
    """ Constructor, takes the xml file path. """
  
    tree = ET.ElementTree(file=filename)
    self.xmlroot = tree.getroot()
    
    # For convenience, store groups as a list and provide a member function
    # to access by name - same for regions.
    self.groups = []
    for group in list(self.xmlroot[0]):
      self.groups.append(SPECSGroup(group))
    

class SPECSGroup:
  """ Encapsulates a "RegionGroup" struct from the SPECS XML format. """
  
  def __init__(self, xmlgroup):
    
    self.name = xmlgroup[0].text
    
    if DEBUG:
      print "======================= ", self.name, " ========================"
    
    self.regions = []
    for region in list(xmlgroup[1]):
      self.regions.append(SPECSRegion(region))
    
class SPECSRegion:
  """ Encapsulates a "RegionData" struct from the SPECS XML format. """
  
  def __init__(self, xmlregion):
    
    self.name = xmlregion[0].text
    self.num_cycles = int(xmlregion[7].attrib['length'])
    
    self.raw_counts = []
    self.scaling_factors = []
    self.extended_channels = []
    for elem in xmlregion.iter('sequence'):
      if elem.attrib['type_name'] == "CountsSeq":
        tmp = array([int(x) for x in elem[0].text.split()])
        self.raw_counts.append(tmp)
      elif elem.attrib['name'] == "scaling_factors":
        tmp = array([float(x) for x in elem[0].text.split()])
        self.scaling_factors.append(tmp)
      elif elem.attrib['type_name'] == "YCurveSeq":
        # Here we need to iterate through the sub-element YCurveSeq
        # to grab the extended channels if they are present.
        for channel in elem.iter('sequence'):
          if channel.attrib['name'] == "data":
            tmp = array([float(x) for x in channel[0].text.split()])
            self.extended_channels.append(tmp)
    
    # Iterate over all the elements in the RegionDef struct.
    for elem in xmlregion[1].iter():
      if elem.attrib['name'] == "scan_mode":
        self.scan_mode = elem[0].text
      elif elem.attrib['name'] == "dwell_time":
        self.dwell_time = float(elem.text)
      elif elem.attrib['name'] == "scan_delta":
        self.scan_delta = float(elem.text)
      elif elem.attrib['name'] == "excitation_energy":
        self.excitation_energy = float(elem.text)
      elif elem.attrib['name'] == "pass_energy":
        self.pass_energy = float(elem.text)
      elif elem.attrib['name'] == "kinetic_energy":
        self.kinetic_energy = float(elem.text)
      elif elem.attrib['name'] == "values_per_curve":
        self.values_per_curve = int(elem.text)
      elif elem.attrib['name'] == "effective_workfunction":
        self.effective_workfunction = float(elem.text)
    
    # The kinetic energy and binding energy axes:
    ke_upper = self.kinetic_energy + (self.values_per_curve - 1) * self.scan_delta
    self.kinetic_axis = linspace(self.kinetic_energy, ke_upper, self.values_per_curve)
    self.binding_axis = self.excitation_energy - self.kinetic_axis
    
    # MCD head and tail are the extra elements added to the beginning and
    # end of the scan.
    self.mcd_head = int(xmlregion[2].text)
    self.mcd_tail = int(xmlregion[3].text)
    
    # Detector channel offsets are in the 5th subelement of the region.
    self.detector_channel_shifts = []
    self.detector_channel_positions = []
    self.detector_channel_gains = []
    for elem in xmlregion[4].iter('struct'):
      if elem.attrib['type_name'] == "Detector":
        for subelem in elem.iter():
          if "name" in subelem.attrib.keys():
            if subelem.attrib['name'] == "position":
              self.detector_channel_positions.append(float(subelem.text))
            if subelem.attrib['name'] == "shift":
              self.detector_channel_shifts.append(float(subelem.text))
            if subelem.attrib['name'] == 'gain':
              self.detector_channel_gains.append(float(subelem.text))
    
    self.detector_channel_shifts = array(self.detector_channel_shifts)
    self.detector_channel_positions = array(self.detector_channel_positions)
    self.detector_channel_gains = array(self.detector_channel_gains)
    
    # Use the pass energy to calculate detector calibration.
    self.detector_channel_offsets = self.pass_energy * self.detector_channel_shifts
    
    # Use interpolators to overlay the 9 detector channels and sum.
    self.counts = zeros((self.values_per_curve))
    self.channel_counts = zeros((self.values_per_curve, len(self.detector_channel_offsets)))
    
    num_detectors = len(self.detector_channel_offsets)
    for c in self.raw_counts:
      if DEBUG:
        print self.name, len(c)
      for i,s in enumerate(self.detector_channel_offsets):
        y = c[i:len(c):num_detectors]
        x = self.kinetic_axis + self.detector_channel_offsets[i]
        # Note use of mcd_head and mcd_tail here.
        if self.mcd_tail == 0:
          s = UnivariateSpline(x, y[self.mcd_head:len(y)],k=1)
        else:
          s = UnivariateSpline(x, y[self.mcd_head:-self.mcd_tail], k=1)
        self.counts += array([s(t) for t in x])
        self.channel_counts[:,i] = array([s(t) for t in x])
      
    # Extract the comment from the parameter list.
    for elem in xmlregion[9].iter("struct"):
      if elem[0].text == "Comment":
        self.comment = elem[1].text
    
    
    
    
        
