#!/usr/bin/env python

from Tkinter import *
import tkFileDialog
import csv
import specs
import os
import os.path
from numpy import zeros, savetxt, amin
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
        # We want to make a folder matching the group name, but we might not
        # be able to because it might already exist (group names are not
        # required to be unique in SPECS). So we have to fiddle here a bit.
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
          num_outputs += 4 # Linear preedge, background
          # scaling factor, shirley and counts - L - S.
          data = zeros((r.values_per_curve,num_outputs))
          hdr_string = ""
          if r.scan_mode == "FixedAnalyzerTransmission":
            data[:,0] = r.binding_axis
            x = r.binding_axis
            hdr_string += "Binding_energy_(eV)    "
          elif r.scan_mode == "ConstantFinalState":
            data[:,0] = r.excitation_axis
            x = r.excitation_axis
            hdr_string += "Photon_energy_(eV)    "
          else:
            data[:,0] = r.kinetic_axis
            x = r.kinetic_axis
            hdr_string += "Kinetic_energy_(eV)    "
          
          data[:,1] = r.counts
          y = r.counts
          hdr_string += "Counts    "
          if have_nine_channels:
            data[:,2:chan_out_idx] = r.channel_counts
            for i in range(2,chan_out_idx):
              hdr_string += "Channeltron_#%d    " % (i-1)
          if r.extended_channels is not None:
            data[:,ex_out_sidx:ex_out_fidx] = r.extended_channels
            for i in range(1,ex_out_fidx-ex_out_sidx+1):
              hdr_string += "Extended_Channel_#%d    " % (i)
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
          # The order of operations here is:
          # 
          # 1. Calculate the linear pre-edge (gives L).
          # 2. Subtract the linear pre-edge (With next step gives y1).
          # 3. Divide by an intensity offset to normalize against changes
          #    in the overall intensity of the background, eg. from analyzer.
          #    We take this offset to be the value of the linear pre-edge at the
          #    lowest energy in the region. This could be more accurate but it 
          #    is pretty good. This gives y1.
          #    (We skip if the background is zero or contains zeros)
          # 4. Calculate Shirley background using this new spectrum (gives S).
          # 5. Subtract the Shirley background (gives y2).
          
          # Have to check for a zero-division here. If we encounter one,
          # don't normalize agaif,;lnst the pre-edge.
          try:
            L = specs.preedge_calculate(x,y)
            y1 = (y - L) / L[x.argmin()]
            S = specs.shirley_calculate(x,y1)
            y2 = y1 - S
            data[:,-4] = L
            data[:,-3] = L[x.argmin()]
            data[:,-2] = S
            data[:,-1] = y2
          except FloatingPointError:
            data[:,-4] = specs.preedge_calculate(x,y)
            data[:,-3] = 1.0
            data[:,-2] = specs.shirley_calculate(x,y - data[:,-4])
            data[:,-1] = y - data[:,-4] - data[:,-2]
          hdr_string += "Preedge    Scale_Division    Shirley    Counts-Preedge-Shirley    "
          
          # As before, there are possible filename conflicts because there is 
          # no requirement that SPECS region names be unique.
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
          
          try:
            savetxt(r_write_name, data,header=hdr_string)
          except TypeError:
            # Numpy version not recent enough for the header option.
            print ("""Your numpy version (%s) is not recent enough to support the\nheader option in savetxt: upgrade to a version >= 1.7 to get\nthis functionality.""" % (npyversion))
            savetxt(r_write_name, data)
            
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
    
