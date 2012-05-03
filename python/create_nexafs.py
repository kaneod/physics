
import esc_lib as el
import libdft as dft
from numpy import *
from pylab import *

reload(el) # Just in case...

gw = 0.02 # Gaussian smearing, Ha
lw = 0.005 # Lorentzian smearing, Ha

title="nc-pT"
core_wf_file = "/home/kane/Desktop/Calculations/atomic_data/carbon/C.GGA-PBE-corewf-neutral.abinit"
ae_wf_file = "/home/kane/Desktop/Calculations/porphyrins/nexafs/poro_AE_WFK-etsf.nc"

c = el.Atom(core_wf_file)
d = el.Atoms(ae_wf_file, "ETSF")

core_index = 0    # Atomic core level to use from the core_wf file.
# For nc-core
site_indices = [73, 78, 79, 75,66,68]
# SITES: [RCCC, RCCN, RCCH, LCCCO, LCCCI, LCHHH, LCCH, LCCCoppI] 
#site_indices = [54,55,42,58,60,64]     # for pT
## SITES: [RCCN, RCHC, RCCC, LCCC, LCHC, LCHHH]
  
print "Hello!"
print ""
print "Using state %d from core wavefunction file: %s" % (core_index, core_wf_file)

for site_index in site_indices:
  print "Placing core state onto atomic site %d from file: %s..." % (site_index, ae_wf_file),
  # Note: we always publically count atomic sites from 1. Need to adjust for python arrays.
  
  # Expand the core wf onto our lattice. This is slow, so get it out of the way
  # now.
  c1s = c.expandOntoGrid(core_index, d.ngrid, d.lattice[0], d.positions[0][site_index-1])
  print "done!"

  print "Initializing Fortran optical module...",
  avec = array(d.lattice[0])
  bvec = array(el.recip_lattice(d.lattice[0]))
  dft.optical.init(c1s, avec, bvec, d.positions[0][site_index-1], 1.3)
  print "done!"

  efs = d.filehook.variables['eigenvalues'][:] # spin, kpt, band
  ei = c.data[core_index]['core E']
  occf = d.filehook.variables['occupations'][:]
  occi = c.data[core_index]['core occ']

  Eall = []
  Iall = []
  Igall = []
  for oa in [1,2,3]: # The optical axes
    print "Calculating for optical axis %d." % (oa)
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
  print "Writing to file nexafs_output_raw_site%d.dat..." % (site_index),

  f = open(title+"_raw_site%d.dat" % (site_index), 'w')
  f.write("# NEXAFS raw eigenvalue/OME data for atomic site %d.\n" % (site_index))
  f.write("# Columns:\n")
  f.write("# energy ome intensity\n")
  for i,dt in enumerate(all_data):
    f.write("# Data set %d\n" % (i+1))
    for ie, ii, ig in zip(dt[0], dt[1], dt[2]):
      f.write("%.15g\t%.15g\t%.15g\n" % (ie, ig, ii))
    f.write("\n")

  f.close()
  
  print "done"
  print "Writing smoothed output to file nexafs_output_smooth_site%d.dat..." % (site_index),
  f = open(title+"_smooth_site%d.dat" % (site_index), 'w')
  f.write("# NEXAFS smoothed spectra for atomic site %d.\n" % (site_index))
  f.write("# Smoothing parameters - Gaussian width %f Ha - Lorentzian width %f Ha.\n" % (gw, lw))
  f.write("# Columns:\n")
  f.write("# energy (Ha) intensity (arb)\n")
  energy1 = linspace(min(all_data[0][0]) -1, max(all_data[0][0]) + 1, 4000)
  energy2 = linspace(min(all_data[1][0]) -1, max(all_data[1][0]) + 1, 4000)
  energy3 = linspace(min(all_data[2][0]) -1, max(all_data[2][0]) + 1, 4000)
  intensity1 = el.gl_smear(all_data[0][0], all_data[0][1], energy1, gw, lw)
  intensity2 = el.gl_smear(all_data[1][0], all_data[1][1], energy2, gw, lw)
  intensity3 = el.gl_smear(all_data[2][0], all_data[2][1], energy3, gw, lw)
  
  f.write("# Data set 1\n")
  for (ee, ii) in zip(energy1, intensity1):
    f.write("%.15g\t%.15g\n" % (ee, ii))
  f.write("\n")
  
  f.write("# Data set 2\n")
  for (ee, ii) in zip(energy2, intensity2):
    f.write("%.15g\t%.15g\n" % (ee, ii))
  f.write("\n")
  
  f.write("# Data set 3\n")
  for (ee, ii) in zip(energy3, intensity3):
    f.write("%.15g\t%.15g\n" % (ee, ii))
  f.write("\n")  
  
  f.close()
  
  print "done"
  print ""

print "Finished calculation!"

