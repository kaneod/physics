#!/usr/bin/env python

from Tkinter import *
import tkFileDialog
import csv
import specs

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
    
    
root = Tk()
app = Application(master = root)
app.mainloop()
root.destroy()
    
