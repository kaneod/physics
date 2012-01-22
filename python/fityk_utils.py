# fityk_utils.py
#
# Dirty hacks!

from numpy import array, zeros, max

def scale(infile, outfile=None, scale=None, shift=0.0, multi_offset=None):
    """ scaled = scale(infile, outfile=None, scale=None, shift=0.0)
    
    Opens infile, assumed to be a commentless multicolumn xy text file (ie, exactly what
    fityk exports). Scales and shifts according to the inputs, default is  normalize to 
    unit maximum. Default output name is the same as the infile with a .scaled suffix.
    
    If multi_offset is not None, will separately shift columns 1 and -1 by the value
    of multi_offset. This is so you can separate components of the peak fits from the
    data and actual fits.
    
    """
    
    f = open(infile, 'r')
    lines = f.readlines()
    f.close()
    
    # Make some safe assumptions about the contents - NO COMMENTS!
    rows = len(lines)
    columns = len(lines[0].split())
    
    data = zeros((rows, columns))
    
    for i, line in enumerate(lines):
        data[i,:] = array([float(x) for x in line.split()])
        
    # Use the *fit* column to find the max because the data might be a bit scattery.
    
    if scale is None:
        peak = max(data[:,2])
    else:
        peak = 1.0 / scale
    
    scaled = zeros((rows, columns))
    scaled[:,0] = data[:,0]
    
    for i in range(1, data.shape[1]):
        if multi_offset is None:
            scaled[:,i] = data[:,i] / peak + shift
        else:
            if i in [1, (data.shape[1]-1)]:
                scaled[:,i] = data[:,i] / peak + multi_offset
            else:
                scaled[:,i] = data[:,i] / peak + shift
            
    if outfile is None:
        f = open(infile+".scaled", 'w')
    else:
        f = open(outfile, 'w')
    
    for i in range(scaled.shape[0]):
        f.write(" ".join([str(x) for x in list(scaled[i,:])])+"\n")
    
    f.close()
    return scaled