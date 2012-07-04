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
# 3. There are functions for things like Shirley background subtraction at the
#    end of the file after the classes. These are tuned to work with the outputs
#    of SPECS but will work more generally if necessary.
#
################################################################################

from __future__ import division
import xml.etree.ElementTree as ET
from numpy import array, linspace, zeros, ceil, amax, amin, argmax, argmin, abs
from numpy import polyfit, polyval, seterr
from numpy.linalg import norm
from scipy.interpolate import UnivariateSpline

DEBUG = 1

# We do not allow divide by zeros at all: raise an error if it happens.
seterr(divide='raise')

################################################################################
#
# CLASSES
#
################################################################################

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
        for ycurve in elem:
          if "Extended Channel" in ycurve[0].text:
            for channel in ycurve.iter('sequence'):
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
    
    # Excitation axis (for NEXAFS)
    exc_upper = self.excitation_energy + (self.values_per_curve - 1) * self.scan_delta
    self.excitation_axis = linspace(self.excitation_energy, exc_upper, self.values_per_curve)
    
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
    
    # Trim the extended channels if they are present.
    for i in range(len(self.extended_channels)):
      if self.mcd_tail == 0:
        c = self.extended_channels[i]
        self.extended_channels[i] = c[self.mcd_head:len(c)]
      else:
        self.extended_channels[i] = self.extended_channels[i][self.mcd_head:-self.mcd_tail]

    # If there are extended channels, reshape them into an array.
    if self.extended_channels:
      if DEBUG:
        print "Extended channels: ", len(self.extended_channels)
        print "Extended channel data length: ", len(self.extended_channels[0])
      tmparr = zeros((len(self.extended_channels[0]), len(self.extended_channels)))
      for i,tmpex in enumerate(self.extended_channels):
        tmparr[:,i] = tmpex
      self.extended_channels = tmparr
    else:
      self.extended_channels = None
          
    # Extract the comment from the parameter list.
    for elem in xmlregion[9].iter("struct"):
      if elem[0].text == "Comment":
        self.comment = elem[1].text
    
################################################################################
#
# FUNCTIONS
#
################################################################################

def preedge_calculate(x,y):
  """ P = specs.preedge_calculate(x,y)
  
  Calculates the best-fit linear pre-edge for a dataset (x,y). Finds the biggest peak,
  then finds the pre-edge region using a sequence of linear fits starting from the end
  point.
  
  """

  # First ensure the energy values are *decreasing* in the array,
  # if not, reverse them. 
  if x[0] < x[-1]:
    is_reversed = True
    x = x[::-1]
    y = y[::-1]
  else:
    is_reversed = False
    
  # Locate the biggest peak.
  maxidx = abs(y - amax(y)).argmin()
  
  # Find the gradient of every possible linear fit between the lowest binding energy
  # and the biggest peak.
  grads = []
  for i in range(2, len(x) - maxidx):    
    # Best linear fit to the last i values
    xs = x[-i:]
    ys = y[-i:]
    p = polyfit(xs,ys,1)
    grads.append(p[0])
    
  # Differentiate the gradient array.
  dgrads = []
  for i in range(len(grads)-1):
    dgrads.append(grads[i+1] - grads[i])
  dgrads = array(dgrads)
  
  # Find the minimum index of the absolute of the gradient of gradients.
  mingrad = abs(dgrads).argmin()
  
  # Make a best linear fit from this number of pre-edge points, generate linear
  # pre-edge.
  p = polyfit(x[-mingrad:], y[-mingrad:], 1)
  
  if is_reversed:
    return polyval(p,x)[::-1]
  else:
    return polyval(p,x)
  
def shirley_calculate(x,y, tol=1e-6, maxit=20):
  """ S = specs.shirley_calculate(x,y, tol=1e-6)
  
  Calculate the best auto-Shirley background S for a dataset (x,y). Finds the biggest peak 
  and then uses the minimum value either side of this peak as the terminal points of the
  Shirley background.
  
  The tolerance sets the convergence criterion, maxit sets the maximum number
  of iterations.
  
  """
  
  # First ensure the energy values are *decreasing* in the array,
  # if not, reverse them.  
  if x[0] < x[-1]:
    is_reversed = True
    x = x[::-1]
    y = y[::-1]
  else:
    is_reversed = False
    
  # Locate the biggest peak.
  maxidx = abs(y - amax(y)).argmin()
  
  # It's possible that maxidx will be 0 or -1. If that is the case,
  # we can't use this algorithm, we return a zero background.
  if maxidx == 0 or maxidx >= len(y)-1:
    print "specs.shirley_calculate: Boundaries too high for algorithm: returning a zero background."
    return zeros(x.shape)
  
  # Locate the minima either side of maxidx.
  lmidx = abs(y[0:maxidx] - amin(y[0:maxidx])).argmin()
  rmidx = abs(y[maxidx:] - amin(y[maxidx:])).argmin() + maxidx
  xl = x[lmidx]
  yl = y[lmidx]
  xr = x[rmidx]
  yr = y[rmidx]
  
  # Max integration index
  imax = rmidx - 1

  # Initial value of the background shape B. The total background S = yr + B,
  # and B is equal to (yl - yr) below lmidx and initially zero above.
  B = zeros(x.shape)
  B[:lmidx] = yl - yr
  Bnew = B.copy()
  
  it = 0
  while it < maxit:  
    # Calculate new k = (yl - yr) / (int_(xl)^(xr) J(x') - yr - B(x') dx')
    ksum = 0.0
    for i in range(lmidx,imax):
      ksum += (x[i] - x[i+1]) * 0.5 * (y[i] + y[i+1] - 2 * yr - B[i] - B[i+1])
    k = (yl - yr) / ksum
    # Calculate new B
    for i in range(lmidx,rmidx):
      ysum = 0.0
      for j in range(i,imax):
        ysum += (x[j] - x[j+1]) * 0.5 * (y[j] + y[j+1] - 2 * yr - B[j] - B[j+1])
      Bnew[i] = k * ysum
    # If Bnew is close to B, exit.
    if norm(Bnew - B) < tol:
      B = Bnew.copy()
      break
    else:
      B = Bnew.copy()
    it += 1
  
  if it >= maxit:
    print "specs.shirley_calculate: Max iterations exceeded before convergence." 
  if is_reversed:
    return (yr+B)[::-1]
  else:    
    return yr + B    
    
    
        
