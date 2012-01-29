# Read in a raw position file, shift the positions by
# a fixed amount, write to a new raw positions file.

from numpy import array

# Input file name
infile = "/home/kane/Desktop/raw.dat"

# Output file name
outfile = "/home/kane/Desktop/shifted.dat"

# Shift vector
shift = array([4.0, 4.0, 2.0])

# Data array
data = []

# MAIN CODE
f = open(infile, 'r')
lines = f.readlines()
f.close()

for line in lines:
    data.append([float(x) for x in line.split()])
    
data = array(data) + shift
f = open(outfile, 'w')

for d in data:
    f.write("    %g    %g    %g\n" % (d[0], d[1], d[2]))

f.close() 

