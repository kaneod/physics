#!/usr/bin/env python
################################################################################
#
# bamski
#
# Program to do basic nexafs data processing of .xy files that come from the
# specs_split unpacking of a SPECS .xml file.
#
################################################################################
#
# Copyright 2012 Kane O'Donnell
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# 
# NOTES
#
# 1. To normalize against a aey spectrum or against the photodiode, use
#    -n aey or -n photodiode respectively. TEY is the default norm channel.
#
################################################################################

from __future__ import division
from sys import exit
import tkFileDialog
import tkMessageBox
import os
import specs
import argparse

# Change this to specific later
from numpy import *
from scipy.interpolate import interp1d
from Tkinter import *
#from pylab import *

# Change to 0 to turn off debug statements.
DEBUG = 1

class Application(Frame):

  def __init__(self, master=None):
    
    parser = argparse.ArgumentParser(description="Auto-double normalize NEXAFS spectra.")
    parser.add_argument('-n', '--norm_channel', help="Choose TEY, AEY or photodiode channel for double normalization.", default="tey")
    parser.add_argument('-p', '--preedge', help="Choose either linear (default) for linear pre-edge subtraction, or none.", default="linear")
    parser.add_argument('-o', '--omit_aey_channels', nargs='+', type=int, default=[], help="Specify channels to OMIT from the aey total spectrum.")
    
    args = parser.parse_args()
    if args.norm_channel.lower() in ["aey", "photodiode", "tey"]:
      self.norm_channel = args.norm_channel.lower()
    else:
      print "ERROR: Option not recognized for norm channel. Must be tey, aey or photodiode. Exiting..."
      sys.exit(0)
    if args.preedge.lower() in ["none", "linear"]:
      self.preedge_treatment = args.preedge.lower()
    else:
      print "Option not recognized for preedge treatment. Must be linear or none."
    
    # Channel omission arguments: must be less than 9 of them, unique and they have to be
    # each integers >= 1 and <= 9.
    if list(unique(args.omit_aey_channels)) != args.omit_aey_channels:
      print "ERROR: Your channel omission list is not unique! Exiting..."
      sys.exit(0)
    if len(args.omit_aey_channels) > 0 and len(args.omit_aey_channels) < 9:
      # Check the channels are sensible.
      for c in args.omit_aey_channels:
        if c < 1 or c > 9:
          print "ERROR: Channel %d does not exist - must be in range [1,9]. Exiting..." % c
          sys.exit(0)
      self.omit_aey_channels = args.omit_aey_channels
    elif len(args.omit_aey_channels) >= 9:
      print "ERROR: Too many channels omitted from AEY - must be less than 9. Exiting..."
      sys.exit(0)
    else:
      self.omit_aey_channels = []
    Frame.__init__(self, master)
    self.pack()
    self.createWidgets()
    
  
  def createWidgets(self):
    """ Make the GUI and put it on the screen """
    
    self.btn_Root = Button(self)
    self.btn_Root["text"] = "Choose XML file..."
    self.btn_Root["command"] = self.set_root
    self.btn_Root.pack()
    
    self.do_norm = IntVar()
    self.chk_Norm = Checkbutton(self, text="Double-normalize", 
                                      variable=self.do_norm,
                                      command=self.hit_check)
    self.chk_Norm.select()
    
    self.chk_Norm.pack()
    
    self.btn_Norm = Button(self)
    self.btn_Norm["text"] = "Norm Ref file..."
    self.btn_Norm["command"] = self.set_norm
    self.btn_Norm.pack()
    
    self.btn_Go = Button(self)
    self.btn_Go["text"] = "GO!"
    self.btn_Go["command"] = self.go
    self.btn_Go.pack()
    
    self.btn_Quit = Button(self)
    self.btn_Quit["text"] = "Quit"
    self.btn_Quit["command"] = self.quit
    self.btn_Quit.pack()
    
    self.normfile = None
    self.rootfile = None
    self.defaultdir = None
    self.btn_Go.configure(state=DISABLED)
    
  def hit_check(self):
    """ User hit checkbutton: controls the NORMAL/DISABLED state of 
        some of the GUI elements. 
    """
    
    if self.do_norm.get():
      self.btn_Norm.configure(state=NORMAL)
    else:
      self.btn_Norm.configure(state=DISABLED)
      
    if self.rootfile and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
      self.btn_Go.configure(state=NORMAL)
    else:
      self.btn_Go.configure(state=DISABLED)
      
  def set_root(self):
    """ User clicked the Set Root button: ask user to select a  SPECS XML file. Some
        checking of appropriate NORMAL/DISABLED states is required.
    
    """
  
    if self.defaultdir == None:
      initdir="~/Desktop"
    else:
      initdir=self.defaultdir
    opts = {}
    opts["filetypes"] = [('All files', '.*'), ('XML files', '.xml')]
    
    self.rootfile = tkFileDialog.askopenfilename(title="Choose the SPECS XML file for processing...", initialdir=initdir)
    
    if self.rootfile and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
      self.defaultdir = os.path.split(self.rootfile)[0]
      self.btn_Go.configure(state=NORMAL)
    else:
      self.btn_Go.configure(state=DISABLED)
  
  def set_norm(self):
    """ User clicked the Set Norm button: ask the user to locate an .xy file.
        We don't check whether it is valid at this stage. Some checking of 
        appropriate NORMAL/DISABLED states is required for some of the buttons.
    
    """
    
    if self.defaultdir == None:
      initdir = "~/Desktop"
    else:
      initdir = self.defaultdir
      
    opts = {}
    opts["filetypes"] = [('All files', '.*'), ('.xy files', '.xy')]
    
    self.normfile = tkFileDialog.askopenfilename(title="Choose a .xy file to normalize against...",     
                                            initialdir=initdir, **opts)
    
    if self.rootfile and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
      self.defaultdir = os.path.split(self.normfile)[0]
      self.btn_Go.configure(state=NORMAL)
    else:
      self.btn_Go.configure(state=DISABLED)
  
  def go(self):
    """ User clicked the Go button: need to check whether the normalization
        file is ok or not, then process each subfolder of the root folder.
    
    """
    
    # First check the normfile, if it exists, is sensible (if we are 
    # double normalizing!)
    
    if self.do_norm.get():
      try:
        nfdata = loadtxt(self.normfile)
      except ValueError:
        tkMessageBox.showerror("Double-normalization file", "The file %s is not of the required .xy format." % self.normfile)
        return
    else:
      nfdata = None
    
    if DEBUG:
      print "Using root file: " + self.rootfile
      print "Current directory is: " + os.getcwd()
      
    
    self.specsobj = specs.SPECS(self.rootfile)
    
    ## Get the list of subfolders in the root folder.
    #self.olddir = os.getcwd()
    #os.chdir(self.rootfile)
    #allpaths = os.listdir(os.getcwd())                                     
    
    #for path in allpaths:
    #  if os.path.isdir(path):
    #    # Process this directory
    #    self.process_subfolder(os.path.join(self.rootfolder, path), nfdata)
    
    #os.chdir(self.olddir)
    
    if self.specsobj:
      self.process_file(nfdata)
      
    print "# ------------------------------------------------------------------"
    print "# "
    print "# Finished processing."
    print "# "
    print "# ------------------------------------------------------------------"
    
    
  def process_file(self, normfiledata):
    """ Process a given SPECS file. Iterate through the regions and if the scan 
        mode is ConstantFinalState we normalize it.
    
    """
    
    # Try to make a folder called rootfile-NEXAFS."
    dfname = self.rootfile + "-NEXAFS"
    if os.path.exists(dfname):
      print "Error: Unpack folder already exists. Delete it and then click Go again."
      return
    else:
      os.mkdir(dfname)
      olddir = os.getcwd()
      os.chdir(dfname)
    
    for g in self.specsobj.groups:
      
      if DEBUG:
        print "=====================", g.name, "====================="
      # Make a folder with the group name - need to make sure it's unique.
      have_name = False
      name_counter = 0
      while not have_name:
        if name_counter == 0:
          tmpname = g.name
        else:
          tmpname = g.name + "-" + str(name_counter)
        if not os.path.exists(tmpname):
          os.mkdir(tmpname)
          os.chdir(tmpname)
          have_name = True
          name_counter = 0
        else:
          # Darn, the directory prob already exists. 
          name_counter += 1
      
      for r in g.regions:
        
        newdata = None
        
        if r.scan_mode == "ConstantFinalState":
          newdata = self.process_region(r, normfiledata)
          print "Processing region: ", r.name
        else:
          print "Region %s is not NEXAFS - skipping." % (r.name)
        
        if newdata is not None:
          laey = None
          ltey = None
          lpey = None
          if self.preedge_treatment == "linear":
            laey = specs.preedge_calculate(newdata[:,0], newdata[:,9])
            ltey = specs.preedge_calculate(newdata[:,0], newdata[:,10])
            lpey = specs.preedge_calculate(newdata[:,0], newdata[:,11])
          
          if laey is None:
            laey = newdata[0,9] * ones(newdata[:,9].shape)
          if ltey is None:
            ltey = newdata[0,10] * ones(newdata[:,10].shape)
          if lpey is None:  
            lpey = newdata[0,11] * ones(newdata[:,11].shape)
          subdata = zeros((newdata.shape[0], newdata.shape[1] + 6))
          subdata[:,:12] = newdata
          subdata[:,12] = laey
          subdata[:,13] = ltey
          subdata[:,14] = lpey
          try:
            subdata[:,15] = newdata[:,9] - laey
          except TypeError:
            subdata[:,15] = newdata[:,9]
            print "Couldn't subtract a linear pre-edge from the aey spectrum for this region."
          try:
            subdata[:,16] = newdata[:,10] - ltey
          except TypeError:
            subdata[:,16] = newdata[:,10]
            print "Couldn't subtract a linear pre-edge from the tey spectrum for this region."
          try:
            subdata[:,17] = newdata[:,11] - lpey
          except TypeError:
            subdata[:,17] = newdata[:,11]
            print "Couldn't subtract a linear pre-edge from the pey spectrum for this region."
          
          # Individually normalize the components
          #try:
          #  subdata[:,15] = subdata[:,15] / subdata[-1,15]
          #except FloatingPointError:
          #  print "The AEY channel cannot be 1-normalized for this file."
          #try:
          #  subdata[:,16] = subdata[:,16] / subdata[-1,16]
          #except FloatingPointError:
          #  print "The TEY channel cannot be 1-normalized for this file."
          #try:
          #  subdata[:,17] = subdata[:,17] / subdata[-1,17]
          #except FloatingPointError:
          #  print "The PEY channel cannot be 1-normalized for this file."
          
          have_name = False
          name_counter = 0
          while not have_name:
            if name_counter == 0:
              r_write_name = r.name + ".xy"
            else:
              r_write_name = r.name + "-" + str(name_counter) + ".xy"
            if not os.path.exists(r_write_name):  
              have_name = True
              name_counter = 0
            else:
              # Darn, the directory prob already exists. 
              name_counter += 1
          
          savetxt(r_write_name, subdata)
      
      # If we didn't actually put any files in this folder, delete it.
      tmpdir = os.getcwd()
      os.chdir(dfname)
      if os.listdir(tmpdir) == []:
        os.rmdir(tmpdir)
    os.chdir(olddir)

      
  def process_region(self, region, normdata):
  
    if DEBUG:
      print "Inside process_region."
    
    # Pull out the channels we need. MUST have extended channels present for NEXAFS
    # work so fail if we don't find them. 
    
    if DEBUG:
      print "Region: ", region.name, region.extended_channels.shape
      
    if region.extended_channels is not None:
      energy = region.excitation_axis
      if len(self.omit_aey_channels) > 0:
        # Need to only select channels in the list.
        ic = list(set([1,2,3,4,5,6,7,8,9]) - set(self.omit_aey_channels))
        aey = zeros(region.counts.shape)
        for i in ic:
          aey += region.channel_counts[:,i-1]
      else:
        aey = region.counts
      tey = region.extended_channels[:,1]
      pey = region.extended_channels[:,5]
      i0 = region.extended_channels[:,2]
    else:
      if DEBUG:
        print "process_region: File not suitable for NEXAFS normalization: no I0 channel!"
      # Need to throw up a warning dialog here!
      return None
    
    # We have to start collecting outputs here. 
    newdata = [energy, aey, tey, pey, i0]
    
    # Normalize all the data streams by i0
    try:
      aey = aey / i0
    except ValueError:
      if DEBUG:
        print "process_region: Divide by zero in i0 normalization - aey is NOT normalized in this data set."
      # Warning or write to log?
    
    try:
      tey = tey / i0
    except ValueError:
      if DEBUG:
        print "process_region: Divide by zero in i0 normalization - tey is NOT normalized in this data set."
      # Warning or write to log?
      
    try:
      pey = pey / i0
    except ValueError:
      if DEBUG:
        print "process_region: Divide by zero in i0 normalization - pey is NOT normalized in this data set."
      # Warning or write to log?    
     
    newdata.append(aey)
    newdata.append(tey)
    newdata.append(pey)
    
    # Now deal with the normalization file if necessary. Same deal.
    if normdata is not None:
      if normdata.shape[1] == 26:
        nenergy = normdata[:,0]
        ntey = normdata[:,12]
        if len(self.omit_aey_channels) > 0:
          # Need to only select channels in the list.
          ic = list(set([1,2,3,4,5,6,7,8,9]) - set(self.omit_aey_channels))
          naey = zeros(normdata[:,1].shape)
          for i in ic:
            naey += normdata[:,i+1]
        else:
          naey = normdata[:,1]
        npd = normdata[:,15]
        ni0 = normdata[:,13]
      elif normdata.shape[1] == 16:
        nenergy = normdata[:,0]
        ntey = normdata[:,3]
        if len(self.omit_aey_channels) > 0:
          # We're fooked because the normalization file doesn't have the channeltrons.
          print "ERROR: Channels were requested to be omitted from the AEY normalization but the norm file doesn't have each separate channeltron. Exiting.."
          sys.exit(0)
        else:
          naey = normdata[:,1]
        npd = normdata[:,6]
        ni0 = normdata[:,4]
      else:
        if DEBUG:
          print "process_region: Normalization file not suitable for NEXAFS normalization: no I0 channel!"
        # Need to throw up a warning dialog here!
        return None      
      
      # Now decide which norm we use.
      if self.norm_channel == "tey":
        if DEBUG:
          print "Using TEY channel for normalization."
        ns = ntey
      elif self.norm_channel == "aey":
        if DEBUG:
          print "Using AEY channel for normalization."
        ns = naey
      elif self.norm_channel == "photodiode":
        if DEBUG:
          print "Using Photodiode channel for normalization."
        ns = npd
      else:
        if DEBUG:
          print "process_region: Unknown norm channel! Must be tey, aey or photodiode."
        return None
      
      try:
        ns = ns / ni0
      except ValueError:
        if DEBUG:
          print "process_region: Divide by zero in i0 normalization - cannot double normalize using this data."
        # Need a warning dialog...
      
      # It's critical that the energy scales are aligned for this kind of normalization
      # otherwise we get serious problems: assume the actual data scale is correct
      # and shift the norm data.
      nmin = argmin(ni0)
      dmin = argmin(i0)
      ecorrect = energy[dmin] - nenergy[nmin]
      nenergy += ecorrect
      # Now normalize by interpolating the norm data. We use a fill value of 
      # positive infinity to ensure that data that is out-of-range of the 
      # normalization energy spectrum is zeroed rather than falsely normalized.
      ni = interp1d(nenergy, ns, bounds_error=False, fill_value=inf)
      nie = array([ni(x) for x in energy])
      newdata.append(nie)

      if 0.0 in nie:
        if DEBUG:
          print "process_region: nie contains a zero - cannot double normalize with this data set."
      else:
        aey = aey / nie
        tey = tey / nie
        pey = pey / nie
        newdata.append(aey)
        newdata.append(tey)
        newdata.append(pey)
    else:
      newdata.append(zeros(aey.shape))
      newdata.append(zeros(aey.shape))
      newdata.append(zeros(aey.shape))
      newdata.append(zeros(aey.shape))
         
    return array(newdata).T
      
        
  
  # Build our little window with the options.
  
root = Tk()
app = Application(master = root)
app.mainloop()
root.destroy()
