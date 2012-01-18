# Converts an abinit HIST file to an animated XSF file for viewing.

from numpy import *
from Scientific.IO import NetCDF

def bohr2ang(bohr):
    return bohr * 0.529177249
    
def ang2bohr(ang):
    return ang / 0.529177249

infilename = "/home/kane/Desktop/Scratch/poro_HIST"
outfilename = "/home/kane/Desktop/Scratch/out.axsf"

# Ideally want to read this automatically in the future, but for now, blah.
ntypat = array([6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, \
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

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
