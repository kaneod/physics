
import xml.etree.ElementTree as ET
from numpy import *

class SPECS:

  def __init__(self, filename):
  
    # Idea here is to keep the xmlroot accessible but only
    # read from it as necessary.
    tree = ET.ElementTree(file=filename)
    self.xmlroot = tree.getroot()
    self.groups = []
    
    for group in list(self.xmlroot[0]):
      self.groups.append(SPECSGroup(group))
    

class SPECSGroup:
  """ Stores data about a group in a SPECS XML file. """
  
  def __init__(self, xmlgroup):
    
    self.name = xmlgroup[0].text
    self.regions = []
    
    for region in list(xmlgroup[1]):
      self.regions.append(SPECSRegion(region))
    
class SPECSRegion:
  """ Stores data about a region in a SPECS XML file. """
  
  def __init__(self, xmlregion):
    
    self.name = xmlregion[0].text
    self.num_cycles = int(xmlregion[7].attrib['length'])
    
    self.counts = []
    self.scaling_factors = []
    self.extended_channels = []
    for elem in xmlregion.iter('sequence'):
      if elem.attrib['type_name'] == "CountsSeq":
        tmp = array([int(x) for x in elem[0].text.split()])
        self.counts.append(tmp)
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
  

    
        