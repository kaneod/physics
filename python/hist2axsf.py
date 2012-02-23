# Converts an abinit HIST file to an animated XSF file for viewing.

################################################################################
#
# hist2axsf.py
#
# Converts an abinit HIST file to an animated XSF file for viewing.
#
################################################################################
#
# Created Jan2012 by Kane O'Donnell (Australian Synchrotron).
#
################################################################################
#
# Copyright 2012 Kane O'Donnell
#
#     This library is free software: you can redistribute it and/or modify
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
# 1. To get a HIST output file from abinit, you have to compile it with respect
# to netCDF and enable the trio output to include netcdf.
#
# 2. This really should be converted to a standalone script with CLI inputs
# rather than values you edit as below...
#
# 3. The functionality of hist2axsf should eventually also be absorbed into
# esc_convert rather than sitting here as it just duplicates some things already
# sorted out in esc_lib.
#
################################################################################

from numpy import *
from Scientific.IO import NetCDF

def bohr2ang(bohr):
    return bohr * 0.529177249
    
def ang2bohr(ang):
    return ang / 0.529177249

infilename = "/home/kane/Desktop/Scratch/poro_HIST"
outfilename = "/home/kane/Desktop/Scratch/out.axsf"

# Ideally want to read this automatically in the future, but for now, blah.
#ntypat = array([6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, \
#            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#ntypat = array(48*[6] + 4 * [7] + 38 * [1])
ntypat  = array([int(x) for x in "1 2 1 1 1 1 2 1 1 1 1 2 1 1 1 1 \
            2 1 1 1 1 2 1 1 1 1 2 1 1 1 1 2 \
            1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 \
            1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 \
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 \
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 \
            3 3".split()])
infile = NetCDF.NetCDFFile(infilename, 'r')
print "Opened NetCDF file %s for reading." % infilename
outfile = open(outfilename, 'w')
print "Opened AXSF file %s for writing." % outfilename
print ""

# Print some info
print "Variable names in this NetCDF are:"
print "\t".join(infile.variables.keys())
print ""
print "Number of MD steps in this file:"
steps = infile.variables['mdtime'].shape[0]
print steps
print ""

# Write AXSF header
outfile.write("ANIMSTEPS %d\n" % steps)
outfile.write("CRYSTAL\n")

# Get the variables we want to write in the animation
xcart = bohr2ang(infile.variables['xcart'].getValue())
rprimd = bohr2ang(infile.variables['rprimd'].getValue())
fcart = ang2bohr(infile.variables['fcart'].getValue())

# Loop over steps
for i in range(steps):
    
    # Write cell lattice at each step
    outfile.write("PRIMVEC %d\n" % (i+1))
    outfile.write("    %g    %g    %g\n" % (rprimd[i,0,0], rprimd[i,1,0], rprimd[i,2,0]))
    outfile.write("    %g    %g    %g\n" % (rprimd[i,0,1], rprimd[i,1,1], rprimd[i,2,1]))
    outfile.write("    %g    %g    %g\n" % (rprimd[i,0,2], rprimd[i,1,2], rprimd[i,2,2]))
    # Write atomic positions and forces at each step
    outfile.write("PRIMCOORD %d\n" % (i+1))
    outfile.write("%d    1\n" % xcart.shape[1])
    for j in range(xcart.shape[1]):
        outfile.write("%d    %g    %g    %g    %g    %g    %g\n" % (ntypat[j], \
                        xcart[i, j, 0], xcart[i, j, 1], xcart[i, j, 2], \
                        fcart[i, j, 0], fcart[i, j, 1], fcart[i, j, 2]))

outfile.write("\n")
outfile.close()
infile.close()
