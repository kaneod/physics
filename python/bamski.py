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
#
#
################################################################################

from __future__ import division
from sys import exit
import tkFileDialog
import tkMessageBox
import os

# Change this to specific later
from numpy import *
from scipy.interpolate import interp1d
from Tkinter import *

# Change to 0 to turn off debug statements.
DEBUG = 1

class Application(Frame):

  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.pack()
    self.createWidgets()
  
  def createWidgets(self):
    """ Make the GUI and put it on the screen """
    
    self.btn_Root = Button(self)
    self.btn_Root["text"] = "Choose root..."
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
    self.rootfolder = None
    self.btn_Go.configure(state=DISABLED)
    
  def hit_check(self):
    """ User hit checkbutton: controls the NORMAL/DISABLED state of 
        some of the GUI elements. 
    """
    
    if self.do_norm.get():
      self.btn_Norm.configure(state=NORMAL)
    else:
      self.btn_Norm.configure(state=DISABLED)
      
    if self.rootfolder and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
      self.btn_Go.configure(state=NORMAL)
    else:
      self.btn_Go.configure(state=DISABLED)
      
  def set_root(self):
    """ User clicked the Set Root button: ask user to select a folder. Some
        checking of appropriate NORMAL/DISABLED states is required.
    
    """
  
    self.rootfolder = tkFileDialog.askdirectory(title="Choose the root SPECS unpacked folder...", initialdir="~/Desktop")
    
    if self.rootfolder and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
      self.btn_Go.configure(state=NORMAL)
    else:
      self.btn_Go.configure(state=DISABLED)
  
  def set_norm(self):
    """ User clicked the Set Norm button: ask the user to locate an .xy file.
        We don't check whether it is valid at this stage. Some checking of 
        appropriate NORMAL/DISABLED states is required for some of the buttons.
    
    """
    
    if self.rootfolder is not None:
      initdir = self.rootfolder
    else:
      initdir = "~/Desktop"
      
    opts = {}
    opts["filetypes"] = [('All files', '.*'), ('.xy files', '.xy')]
    
    self.normfile = tkFileDialog.askopenfilename(title="Choose a .xy file to normalize against...",     
                                            initialdir=initdir, **opts)
    
    if self.rootfolder and ((self.do_norm.get() and self.normfile) or not self.do_norm.get()):
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
      print "Using rootfolder: " + self.rootfolder
      print "Current directory is: " + os.getcwd()
    # Get the list of subfolders in the root folder.
    self.olddir = os.getcwd()
    os.chdir(self.rootfolder)
    allpaths = os.listdir(os.getcwd())                                     
    
    for path in allpaths:
      if os.path.isdir(path):
        # Process this directory
        self.process_subfolder(os.path.join(self.rootfolder, path), nfdata)
    
    os.chdir(self.olddir)
       
    print "# ------------------------------------------------------------------"
    print "# "
    print "# Finished processing."
    print "# "
    print "# ------------------------------------------------------------------"
    
    
  def process_subfolder(self, subfolder, normfiledata):
    """ Process a given subfolder: for each file, if it ends in .xy try
        to load it. If it loads, send it off for processing. If not, 
        throw a warning dialog and carry on with the next file.
    """
  
    olddir = os.getcwd()
    os.chdir(subfolder)
    allpaths = os.listdir(subfolder)
    
    if DEBUG:
      print "Processing subfolder: " + subfolder
      
    for path in allpaths:
      if os.path.isfile(path) and path.endswith(".xy"):
        try:
          data = loadtxt(path)
        except ValueError:
          tkMessageBox.showwarning("Invalid Data File", "The file %s is not of the required .xy format. Skipping..." % os.path.join(os.getcwd(), path))
          break
        
        if DEBUG:
          print "Processing file: " + path
        newdata = self.process_xy(data, normfiledata)
        
        if newdata is not None:
          newfname = path.rstrip(".xy") + "-normed.xy"
          savetxt(newfname, newdata)
     
    os.chdir(olddir)
      
  def process_xy(self, data, normdata):
  
    if DEBUG:
      print "Inside process_xy."
    
    # What's in the file? There are only a few possibilities. If the data has
    # 25 columns, we have full output (channeltrons + extended). If the data
    # has 16 columns, we don't have the channeltrons (who cares), but if we
    # have 15 columns we can't normalize because we have the channeltrons but
    # not the extended channels. If we only have 6 columns, there is nothing
    # included except the raw data plus the analysis - not useful.
    
    if data.shape[1] == 25:
      energy = data[:,0]
      aey = data[:,1]
      tey = data[:,12]
      i0 = data[:,13]
      pey = data[:,16]
    elif data.shape[1] == 16:
      energy = data[:,0]
      aey = data[:,1]
      tey = data[:,3]
      i0 = data[:,4]
      pey = data[:,7]
    else:
      if DEBUG:
        print "process_xy: File not suitable for NEXAFS normalization: no I0 channel!"
      # Need to throw up a warning dialog here!
      return None
    
    # We have to start collecting outputs here. 
    newdata = [energy, aey, tey, pey, i0]
    
    # Normalize all the data streams by i0
    try:
      aey = aey / i0
    except ValueError:
      if DEBUG:
        print "process_xy: Divide by zero in i0 normalization - aey is NOT normalized in this data set."
      # Warning or write to log?
    
    try:
      tey = tey / i0
    except ValueError:
      if DEBUG:
        print "process_xy: Divide by zero in i0 normalization - tey is NOT normalized in this data set."
      # Warning or write to log?
      
    try:
      pey = pey / i0
    except ValueError:
      if DEBUG:
        print "process_xy: Divide by zero in i0 normalization - pey is NOT normalized in this data set."
      # Warning or write to log?    
     
    newdata.append(aey)
    newdata.append(tey)
    newdata.append(pey)
    
    # Now deal with the normalization file if necessary. Same deal.
    if normdata is not None:
      if normdata.shape[1] == 25:
        nenergy = normdata[:,0]
        ntey = normdata[:,12]
        ni0 = normdata[:,13]
      elif normdata.shape[1] == 16:
        nenergy = normdata[:,0]
        ntey = normdata[:,3]
        ni0 = normdata[:,4]
      else:
        if DEBUG:
          print "process_xy: Normalization file not suitable for NEXAFS normalization: no I0 channel!"
        # Need to throw up a warning dialog here!
        return None      
      
      try:
        ntey = ntey / ni0
      except ValueError:
        if DEBUG:
          print "process_xy: Divide by zero in i0 normalization - cannot double normalize using this data."
        # Need a warning dialog...
      
      # Now normalize by interpolating the norm data. We use a fill value of 
      # positive infinity to ensure that data that is out-of-range of the 
      # normalization energy spectrum is zeroed rather than falsely normalized.
      ni = interp1d(nenergy, ntey, bounds_error=False, fill_value=inf)
      nie = array([ni(x) for x in energy])
      newdata.append(nie)

      if 0.0 in nie:
        if DEBUG:
          print "process_xy: nie contains a zero - cannot double normalize with this data set."
      else:
        aey = aey / nie
        tey = tey / nie
        pey = pey / nie
        newdata.append(aey)
        newdata.append(tey)
        newdata.append(pey)
       
    return array(newdata).T
      
        
  
  # Build our little window with the options.
  
root = Tk()
app = Application(master = root)
app.mainloop()
root.destroy()
