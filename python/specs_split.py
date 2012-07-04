#!/usr/bin/env python

from Tkinter import *
import tkFileDialog
import csv
import specs
import os
import os.path
from numpy import zeros, savetxt
from numpy import __version__ as npyversion

DEBUG = 1

class Application(Frame):
  
  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.pack()
    self.createWidgets()
    
  def createWidgets(self):
    
    self.QUIT = Button(self)
    self.QUIT["text"] = "Quit"
    self.QUIT["command"] = self.quit
    
    self.QUIT.pack({"side" : "left"})
    
    self.choose_file = Button(self)
    self.choose_file["text"] = "Choose file..."
    self.choose_file["command"] = self.loadxml
    
    self.choose_file.pack({"side" : "right"})
  
  def loadxml(self):
    
    xmlfile = tkFileDialog.askopenfilename(title="Choose a SPECS .xml file.", initialdir="~/Desktop")
    
    if xmlfile != None:
      xmlcontents = specs.SPECS(xmlfile)
      
      # Make a folder: call it "xmlfile"-unpacked
      base_path = os.path.split(xmlfile)[0]
      cwd = os.getcwd()
      unpack_path = os.path.join(base_path, xmlfile+"-unpacked")
      os.mkdir(unpack_path)
      os.chdir(unpack_path)
      
      for g in xmlcontents.groups:
        if DEBUG:
          print "=================== OUTPUT: %s ====================" % (g.name)
        os.mkdir(g.name)
        os.chdir(g.name)
        for r in g.regions:
          if DEBUG:
            print "Constructing output for region %s..." % (r.name)
          # Decide how many outputs based on whether we have extended channels
          # and so on.
          have_nine_channels = False
          num_outputs = 2 # Bare output is energy and counts
          if r.channel_counts.shape[1] > 1:
            have_nine_channels = True
            num_outputs += r.channel_counts.shape[1]
            chan_out_idx = 2 + r.channel_counts.shape[1]
          if r.extended_channels is not None:
            # Output all extended channels plus counts / I0.
            ex_out_sidx = num_outputs
            num_outputs += r.extended_channels.shape[1] + 1
            ex_out_fidx = ex_out_sidx + 9
          # For now, automatically add the shirley and linear subtractions.
          # Eventually this will be an option.
          num_outputs += 3 # Linear preedge, shirley and counts - P - S.
          data = zeros((r.values_per_curve,num_outputs))
          hdr_string = ""
          if r.scan_mode == "FixedAnalyzerTransmission":
            data[:,0] = r.binding_axis
            x = r.binding_axis
            hdr_string += "Binding energy (eV)    "
          elif r.scan_mode == "ConstantFinalState":
            data[:,0] = r.excitation_axis
            x = r.excitation_axis
            hdr_string += "Photon energy (eV)    "
          else:
            data[:,0] = r.kinetic_axis
            x = r.kinetic_axis
            hdr_string += "Kinetic energy (eV)    "
          
          data[:,1] = r.counts
          y = r.counts
          hdr_string += "Counts    "
          if have_nine_channels:
            data[:,2:chan_out_idx] = r.channel_counts
            for i in range(2,chan_out_idx):
              hdr_string += "Channeltron #%d    " % (i-1)
          if r.extended_channels is not None:
            data[:,ex_out_sidx:ex_out_fidx] = r.extended_channels
            for i in range(1,ex_out_fidx-ex_out_sidx+1):
              hdr_string += "Extended Channel #%d    " % (i)
            # Need to check for divide-by-zero here
            try:
              data[:,ex_out_fidx] = r.counts / r.extended_channels[:,2]
              y = r.counts / r.extended_channels[:,2]
              hdr_string += "Counts/I0    "
            except FloatingPointError:
              print "Divide by zero in I0 channel: skipping I0 normalization."
              data[:,ex_out_fidx] = r.counts
              y = r.counts
              hdr_string += "Counts    "
          data[:,-3] = specs.preedge_calculate(x,y)
          data[:,-2] = specs.shirley_calculate(x,y - data[:,-3])
          data[:,-1] = y - data[:,-3] - data[:,-2]
          hdr_string += "Preedge    Shirley    Counts-Preedge-Shirley    "
          
          try:
            savetxt(r.name+".xy", data,header=hdr_string)
          except TypeError:
            # Numpy version not recent enough for the header option.
            print ("""Your numpy version (%s) is not recent enough to support the\nheader option in savetxt: upgrade to a version >= 1.7 to get\nthis functionality.""" % (npyversion))
            savetxt(r.name+".xy", data)
            
        os.chdir(unpack_path)
      
      os.chdir(cwd)
    
    print "===================================================================="
    print "=                                                                  ="
    print "=                    Finished processing                           ="
    print "=                                                                  ="
    print "===================================================================="
    
root = Tk()
app = Application(master = root)
app.mainloop()
root.destroy()
    
