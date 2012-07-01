#!/usr/bin/env python

from Tkinter import *
import tkFileDialog
import csv
import specs
import os
import os.path
from numpy import *

class Application(Frame):
  
  def __init__(self, master=None):
    Frame.__init__(self, master)
    self.pack()
    self.createWidgets()
    
  def createWidgets(self):
    
    self.QUIT = Button(self)
    self.QUIT["text"] = "QUIT"
    self.QUIT["fg"] = "red"
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
        os.mkdir(g.name)
        os.chdir(g.name)
        for r in g.regions:
          data = zeros((r.values_per_curve,4))
          if r.scan_mode == "FixedAnalyzerTransmission":
            data[:,0] = r.binding_axis
            hdr_string = "Binding energy (eV)\tCounts\tI0\tCounts/I0"
          elif r.scan_mode == "ConstantFinalState":
            data[:,0] = r.excitation_axis
            hdr_string = "Photon energy (eV)\tCounts\tI0\tCounts/I0"
          else:
            data[:,0] = r.kinetic_axis
            hdr_string = "Kinetic energy (eV)\tCounts\tI0\tCounts/I0"
          data[:,1] = r.counts
          data[:,2] = r.extended_channels[2]
          data[:,3] = r.counts / r.extended_channels[2]
          
          savetxt(r.name+".xy", data,header=hdr_string)
        os.chdir(unpack_path)
      
      os.chdir(cwd)
    
root = Tk()
app = Application(master = root)
app.mainloop()
root.destroy()
    
