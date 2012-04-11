
import esc_lib as el
import libdft as dft
from numpy import *
from pylab import *

reload(el) # Just in case...

print "Hello!"
print ""

core_wf_file = "/home/kane/Desktop/Calculations/atomic_data/carbon/C.GGA-PBE-corewf.abinit"
ae_wf_file = "/home/kane/Desktop/Calculations/porphyrins/pT-corehole/poro_AE_WFK-etsf.nc"

c = el.Atom(core_wf_file)
d = el.Atoms(ae_wf_file, "ETSF")

core_index = 0    # Atomic core level to use from the core_wf file.
site_index = 81    # Atomic site to place the core wavefunction on.

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
dft.optical.init(c1s, avec, bvec, d.positions[0][site_index], 1.1)
print "done!"

gw = 0.02 # Gaussian smearing, Ha
lw = 0.005 # Lorentzian smearing, Ha

efs = d.filehook.variables['eigenvalues'][:] # spin, kpt, band
ei = c.data[core_index]['core E']
occf = d.filehook.variables['occupations'][:]
occi = c.data[core_index]['core occ']

Eall = []
Iall = []
Igall = []
for oa in [1,2,3]: # The optical axes
  E = []
  I = []
  Ig = []
  for s in range(d.filehook.dimensions['number_of_spins']):
    for k in range(d.filehook.dimensions['number_of_kpoints']):
      for b in range(d.filehook.dimensions['max_number_of_states']):
        
        if abs(occf[s,k,b] - 2.0) < 1e-8:
          print "State s,k,b = %d %d %d is fully occupied, not calculating." % (s,k,b)
        else:
        
          print "Processing s,k,b: ", s, k, b
          ef = efs[s,k,b]
          ff = occf[s,k,b]
        
          E.append(ef - ei)
        
          # Set wff to the current wavefunction.
          reduced_k = d.filehook.variables['reduced_coordinates_of_kpoints'][k]
          kpt = el.reduced2cart(reduced_k, d.recip_lattice[0])
          dft.optical.set_wff(d.getRealSpaceWF(s,k,b,0), kpt)
          
          # Calculate <core|grad|final>
          ome = dft.optical.optical_matrix_element(oa)
         
          # Divide by energy difference squared
          # <i|del_x|f> / (ei - ef)
          o = ome / (ef - ei) ** 2
        
          # Multiply by the difference between the occupancy of the
          # final state and 2 (electrons) and divide out the
          # double factor.
        
          o = o * (2.0 - occf[s,k,b]) / 2.0
        
          Ig.append(ome)
          I.append(o)

  Eall.append(E)
  Iall.append(I)
  Igall.append(Ig)
  

#plot(E,I, 'rd')
#plot(E,Ig, 'bd')

print ""
print "Finished OME calculations."
print ""
print "Sorting data by energy...",
# Create a dictionary of the E,I values
all_data = []
for (E, I, Ig) in zip(Eall, Iall, Igall):
  EIDict = {}
  EIgDict = {}
  for (ee, ii, ig) in zip(E,I, Ig):
    EIDict[ee] = ii
    EIgDict[ee] = ig

  # Sort by keys
  EE = []
  II = []
  IG = []
  for key in sorted(EIDict.iterkeys()):
    EE.append(key)
    II.append(EIDict[key])
    IG.append(EIgDict[key])
    
  all_data.append([EE,II,IG])

print "done"
print "Writing to file nexafs_output_raw.dat...",
# Construct a linear energy space, 1 hartree either end
#energy = linspace(min(EE) - 1, max(EE) + 1, 4000)
#intensity = el.gl_smear(EE*27.211, II, energy, gw,lw)

f = open("nexafs_output_raw.dat", 'w')
f.write("# NEXAFS raw eigenvalue/OME data.\n")
f.write("# Columns:\n")
f.write("# energy ome_axis1 ome_axis2 ome_axis3 intens_axis1 intens_axis2 intens_axis3\n")
for i in range(len(all_data[0][0])):
  f.write("%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n" % (all_data[0][0][i], all_data[0][2][i], all_data[1][2][i], all_data[2][2][i], all_data[0][1][i], all_data[1][1][i], all_data[2][1][i]))

f.close()
print "done"
print ""
print "Finished calculation!"

