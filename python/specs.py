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
# 1b. This is version 2 of the library - I have tried to be more version-
#    agnostic in the XML parsing here so we aren't as sensitive to file format.
#
# 2. This library is not complete - plan is to incrementally add more member
#    elements to the classes as necessary.
#
# 3. There are functions for things like Shirley background subtraction at the
#    end of the file after the classes. These are tuned to work with the outputs
#    of SPECS but will work more generally if necessary.
#
# 4. The fundamental issue to deal with here is that there are different
#    versions of the SPECS XML format (1.3 and 1.6 are the most prevalent as of
#    the 1st of October 2012) and they store the data in slightly different
#    ways. So we often don't pull things directly from the XML tree by position
#    but rather we search for them or iterate through a list and check props 
#    before we act on a given element.
#
################################################################################

from __future__ import division
import xml.etree.ElementTree as ET
from numpy import array, linspace, zeros, ceil, amax, amin, argmax, argmin, abs
from numpy import polyfit, polyval, seterr, trunc, mean
from numpy.linalg import norm
from scipy.interpolate import interp1d

DEBUG = 1
OPTION = 2

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
    
    # The version impacts on properties of the document so we need to read it
    # here.
    self.xmlversion = self.xmlroot.get('version')
    
    # For convenience, store groups as a list and provide a member function
    # to access by name - same for regions.
    self.groups = []
    for group in list(self.xmlroot[0]):
      # All the subelements will be individual groups (called a RegionGroup in
      # SPECS parlance) but we must check in case the file format changes.
      if group.get('type_name') == "RegionGroup":
        self.groups.append(SPECSGroup(group))
    

class SPECSGroup:
  """ Encapsulates a "RegionGroup" struct from the SPECS XML format. """
  
  def __init__(self, xmlgroup):
    
    self.name = xmlgroup[0].text
    
    if DEBUG:
      print "======================= ", self.name, " ========================"
    
    self.regions = []
    for region in list(xmlgroup[1]):
      if region.get('type_name') == "RegionData":
        self.regions.append(SPECSRegion(region))
    
class SPECSRegion:
  """ Encapsulates a "RegionData" struct from the SPECS XML format. """
  
  def __init__(self, xmlregion):
    
    self.name = xmlregion[0].text
    self.num_cycles = int(xmlregion[7].attrib['length'])
    
    self.raw_counts = []
    self.scaling_factors = []
    self.extended_channels = []
    
    # First grab the counts for this region: this is the most important part.
    # The counts can't be used directly as they incorporate all nine channels
    # in a single array and need to be chopped and aligned first.
    # Improvement from v1: we search directly for the named sequence rather 
    # than iterating generally.
    for elem in xmlregion.findall(".//sequence[@type_name='CountsSeq']"):
      self.raw_counts.append(array([int(x) for x in elem[0].text.split()]))
    
    # Scaling factors for the counts.
    for elem in xmlregion.findall(".//sequence[@name='scaling_factors']"):
      self.scaling_factors.append(array([float(x) for x in elem[0].text.split()]))
    
    # Look for Extended Channels in a YCurveSeq set.
    for ycs in xmlregion.findall(".//sequence[@type_name='YCurveSeq']"):
      for ycurve in ycs:
        if "Extended Channel" in ycurve[0].text:
          for channel in ycurve.iter('sequence'):
            if channel.attrib['name'] == "data":
              tmp = array([float(x) for x in channel[0].text.split()])
              self.extended_channels.append(tmp)
              
    # Grab the transmission function. This is *not* to be trusted but SPECS 
    # might implicitly use it for display within the SPECS program itself so
    # we need to read it.
    trans = xmlregion.find(".//sequence[@name='transmission']")
    self.transmission = array([float(x) for x in trans[0].text.split()])
    
    # Iterate over all the elements in the RegionDef struct.
    # Note: should ONLY BE ONE of these, so use find rather than findall.
    rdef = xmlregion.find(".//struct[@type_name='RegionDef']")
    for elem in rdef:
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
    self.mcd_head = int(xmlregion.find(".//*[@name='mcd_head']").text)
    self.mcd_tail = int(xmlregion.find(".//*[@name='mcd_tail']").text)
    
    # Get the detector information for the energy position of each channeltron.
    self.detector_channel_shifts = []
    self.detector_channel_positions = []
    self.detector_channel_gains = []
    
    detectors = xmlregion.find(".//sequence[@type_name='DetectorSeq']")
    for elem in detectors:
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
    
    num_detectors = len(self.detector_channel_offsets)
    
    # Now, we need to know the analyzer mode, because how we add the channeltron data
    # together depends on whether we are sweeping the kinetic energy in the analyzer or
    # not.
    scanmode = xmlregion.find(".//struct[@type_name='ScanMode']")
    self.scan_mode = scanmode[0].text
    
    # Calculate so and si (based on the SPECS document "Acquiring Data with
    # Multidetector systems"). Don't really need si or t.
    so = self.detector_channel_offsets[-1]
    #si = self.detector_channel_offsets[0]
    h = int(trunc(so / self.scan_delta + 0.5))
    #t = int(trunc(-si / self.scan_delta + 0.5))
    
    # Now use the h value to calculate the index offsets for each of the channels.
    # (This isn't used for ConstantFinalState)
    start_energies = []
    for i in range(num_detectors):
      start_energies.append(self.kinetic_energy - h * self.scan_delta + self.detector_channel_offsets[i])
    idxs = []
    for i in range(num_detectors):
      idxs.append(int(trunc((self.kinetic_energy - start_energies[i]) / self.scan_delta + 0.5)))    
    
    # We now need to separate the raw counts into channels and assign each counts value
    # to a nominal energy value again according to the SPECS document referenced above,
    # using the "Nearest-Neighbour" method.
    self.counts = zeros((self.values_per_curve))
    self.channel_counts = zeros((self.values_per_curve, len(self.detector_channel_offsets)))
  
    for c in self.raw_counts:
      tmp_channels = []   
      for i in range(num_detectors):
        tmp_channels.append(c[i::9])
      # IMPORTANT: If FixedAnalyzerTransmission or FixedRetardingRatio, we need to use
      # the nearest-neighbour method to align the channeltron energies. I have only
      # implemented the method for FixedAnalyzerTransission at the moment - the FRR 
      # implementation is different and rather more difficult and no one ever uses it.
      if self.scan_mode != "FixedAnalyzerTransmission":
        for i in range(num_detectors):
          self.counts += array(tmp_channels[i])
          self.channel_counts[:,i] += array(tmp_channels[i])
      else:
        for i in range(self.values_per_curve):
          for j in range(num_detectors):
            try:
              self.counts[i] += tmp_channels[j][i + idxs[j]]
              self.channel_counts[i, j] += tmp_channels[j][i + idxs[j]]
            except IndexError:
              print "SPECSRegion: Darn, an index error unpacking the channeltron data. This was not supposed to happen!"
                
    # Use the SPECS nearest-neighbour method to recombine the channeltron
    # spectra after extraction from the raw_counts array.
    #for c in self.raw_counts:
    #  tmp_channels = []
    #  tmp_indices = []
    #  tmp_E = self.kinetic_axis[0] # This is E1 from the SPECS document
    #  for i in range(num_detectors):
    #    # Pull out the individual channel spectrum
    #    tmp_channels.append(c[i::num_detectors])
    #    # Note we have to *reverse* the channel offsets here based on the sign
    #    # convention in the SPECS document.
    #    tmp_indices.append(trunc(h + self.detector_channel_offsets[::-1][i] / self.scan_delta + 0.5))
    #  if DEBUG:
    #    print "Starting indices are: ", tmp_indices
    #  # Now tmp_channels[j][k] is the kth value of the jth channeltron. 
    #  # We map these channeltrons together by using the index offsets in
    #  # the tmp_indices list which is incremented every j cycle.
    #  for i in range(self.values_per_curve):
    #    if DEBUG:
    #      print "---------------- Energy index %d %f ------------------" % (i, self.kinetic_axis[i])
    #    for j in range(num_detectors):
    #      if DEBUG:
    #        print "tmp index: ", tmp_indices[j], len(tmp_channels[j])
    #      self.counts[i] += tmp_channels[j][tmp_indices[j]]
    #      self.channel_counts[i,j] += tmp_channels[j][tmp_indices[j]]
    #      tmp_indices[j] += 1
       
    #for c in self.raw_counts:
    #  if DEBUG:
    #    print self.name, len(c)
    #  for i in range(num_detectors):
    #    # :: is not a mistake, shorthand for start:finish:step when we want to
    #    # finish at the end of the array.
    #    y = c[i::num_detectors]
    #    x = self.extended_axis + self.detector_channel_offsets[i]
    #    
    #    if DEBUG:
    #      print "Constructing interpolators: interp range is ", amin(x), amax(x)
    #      print "Kinetic energy axis range is ", amin(self.kinetic_axis), amax(self.kinetic_axis)
    #    # There are sometimes slight energy alignment problems with this axis
    #    # (the subsequent interpolation has strict boundary requirements for
    #    # numerical integrity reasons - we don't want to interpolate outside the
    #    # known range). We correct by forcing the highest energy in the axis to
    #    # be an exact multiple of the scan step size. This is a very small (3rd
    #    # or 4th decimal place) correction, well smaller than the resolution
    #    # of the instrument.roun
    #    shift = self.scan_delta * ceil(amax(x) / self.scan_delta) - amax(x)
    #    x += shift
    #    
    #    if DEBUG:
    #      print "Corrected interpolator range is ", amin(x), amax(x)  
    #    # We might get requested interpolators
    #    s = interp1d(x,y)      
    #    #s = interp1d(x, y, bounds_error=False, fill_value=0.0)
    #    # Now interpolate onto the kinetic_axis for the total counts (all
    #    # channels added together) and for the separate channels. The practical
    #    # effect is to energy-shift each channeltron and chop the energy
    #    # range so all channeltron spectra are over the same axis.
    #    self.counts += array([s(t) for t in self.kinetic_axis])
    #    self.channel_counts[:,i] += array([s(t) for t in self.kinetic_axis])
    
    # Trim the extended channels if they are present. There should not be any
    # calibration issue here - SPECS just treats the extended channels as if 
    # they are channeltrons and therefore gives them extra data points on either
    # side as indicated by MCD head and tail. One could leave the extra points
    # in without any hassle but you would then not match the actual excitation
    # or kinetic energy range specified by the end-user.
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

  # Make sure we've been passed arrays and not lists.
  x = array(x)
  y = array(y)
  
  # Sanity check: Do we actually have data to process here?
  if not (x.any() and y.any()):
    print "specs.preedge_calculate: One of the arrays x or y is empty. Returning zero background."
    return zeros(x.shape)

  # Next ensure the energy values are *decreasing* in the array,
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
    #p = polyfit(xs,ys,1)
    #grads.append(p[0])
    # Try a new algorithm that should be faster than polyfit
    xs = xs - mean(xs)
    ys = ys - mean(ys)
    grads.append((xs * ys).sum() / (xs * xs).sum())
    
  # Differentiate the gradient array.
  dgrads = []
  for i in range(len(grads)-1):
    dgrads.append(grads[i+1] - grads[i])
  dgrads = array(dgrads)
  
  # We might not have actually accumulated anything if the maximum is near the
  # edge (like in a survey scan - the SE background is very big). So, may have
  # to return a zero background.
  if not dgrads.any():
    print "specs.preedge_calculate: No pre-edge gradients. The spectrum must be very large at the low kinetic energy end. Returning zero background."
    return zeros(x.shape)
  
  # Find the minimum index of the absolute of the gradient of gradients.
  mingrad = abs(dgrads).argmin()
  
  # Make a best linear fit from this number of pre-edge points, generate linear
  # pre-edge.
  p = polyfit(x[-mingrad:], y[-mingrad:], 1)
  
  if is_reversed:
    return polyval(p,x)[::-1]
  else:
    return polyval(p,x)
  
def shirley_calculate(x,y, tol=1e-5, maxit=10):
  """ S = specs.shirley_calculate(x,y, tol=1e-5, maxit=10)
  
  Calculate the best auto-Shirley background S for a dataset (x,y). Finds the biggest peak 
  and then uses the minimum value either side of this peak as the terminal points of the
  Shirley background.
  
  The tolerance sets the convergence criterion, maxit sets the maximum number
  of iterations.
  
  """

  # Make sure we've been passed arrays and not lists.
  x = array(x)
  y = array(y)
  
  # Sanity check: Do we actually have data to process here?
  if not (x.any() and y.any()):
    print "specs.shirley_calculate: One of the arrays x or y is empty. Returning zero background."
    return zeros(x.shape)
  
  # Next ensure the energy values are *decreasing* in the array,
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
    if DEBUG:
      print "Shirley iteration: ", it 
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
    
    
        
