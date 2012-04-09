
import esc_lib as el
import libdft as dft
from numpy import *
from pylab import *

reload(el) # Just in case...

print "Hello!"
print ""

core_wf_file = "C.GGA-PBE-corewf.abinit"
ae_wf_file = "o_AE_WFK-etsf.nc"

c = el.Atom(core_wf_file)
d = el.Atoms(ae_wf_file, "ETSF")

core_index = 0    # Atomic core level to use from the core_wf file.
site_index = 1    # Atomic site to place the core wavefunction on.
optical_axis = 1 

print "Using state %d from core wavefunction file: %s" % (core_index, core_wf_file)
print "Placing core state onto site %d from file: %s" % (site_index, ae_wf_file)
print "...",

# Expand the core wf onto our lattice. This is slow, so get it out of the way
# now.
c1s = c.expandOntoGrid(core_index, d.ngrid, d.lattice[0], d.positions[0][site_index])
print "done!"

print "Initializing Fortran optical module...",
avec = array(d.lattice[0])
bvec = array(el.recip_lattice(d.lattice[0]))
dft.optical.init(c1s, avec, bvec)
print "done!"

gw = 0.02 # Gaussian smearing, Ha
lw = 0.02 # Lorentzian smearing, Ha

def gaussian(x,w):
  return 1.0 / (w * sqrt(2*pi)) * exp(-x ** 2 / (2 * w ** 2))
  
def lorentzian(x,w):
  return (1.0/ pi) * w / (x ** 2 + w ** 2)

efs = d.filehook.variables['eigenvalues'][:] # spin, kpt, band
ei = c.data[core_index]['core E']
occf = d.filehook.variables['occupations'][:]
occi = c.data[core_index]['core occ']

E = []
I = []
Ig = []
for s in range(d.filehook.dimensions['number_of_spins']):
  for k in range(d.filehook.dimensions['number_of_kpoints']):
    for b in range(d.filehook.dimensions['max_number_of_states']):
      
      ef = efs[s,k,b]
      ff = occf[s,k,b]
      
      E.append(ef - ei)
      
      # Set wff to the current wavefunction.
      reduced_k = d.filehook.variables['reduced_coordinates_of_kpoints'][k]
      kpt = el.reduced2cart(reduced_k, d.recip_lattice[0])
      dft.optical.set_wff(d.getRealSpaceWF(s,k,b,0), kpt)
      
      # Calculate <core|grad|final>
      ome = dft.optical.optical_matrix_element(optical_axis)
      
      # Divide by energy difference squared
      # <i|del_x|f> / (ei - ef)
      o = ome / (ef - ei) ** 2
      
      # Multiply by the difference between the occupancy of the
      # final state and 2 (ie, 2 electrons) and divide out the
      # double factor.
      
      o = o * (2.0 - occf[s,k,b]) / 2.0
      
      Ig.append(ome)
      I.append(o)

ion()
#plot(E,I, 'rd')
#plot(E,Ig, 'bd')

# Create a dictionary of the E,I values
EIDict = {}
for (ee, ii) in zip(E,I):
  EIDict[ee] = ii

# Sort by keys
EE = []
II = []
for key in sorted(EIDict.iterkeys()):
  EE.append(key)
  II.append(EIDict[key])

#plot(EE,II, 'rd')
# Convolute with a gaussian
Ic = el.gaussian_convolute(array(EE), array(II), 0.01)
#plot(EE,Ic, 'b-')

# Construct a linear energy space, 1 hartree either end
energy = linspace(min(EE) - 1, max(EE) + 1, 1000)
intensity = zeros(energy.shape)

# For each energy, sum a gaussian/lorentzian from each eigenvalue.
for i, en in enumerate(energy):
  for eeig, ieig in zip(EE, II):
  
    # Artificial: ignore silly large values
    if ieig < 10:
      gauss = gaussian(en - eeig, gw)
      lorentz = lorentzian(en - eeig, lw)
      intensity[i] = intensity[i] + ieig * gauss * lorentz
      
plot(energy, intensity, 'r-')