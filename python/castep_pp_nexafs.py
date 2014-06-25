#!/usr/bin/env python

# This script does a bunch of routine post-processing tasks on CASTEP output files
# related to theoretical NEXAFS calculations. It assumes the following:
#
#   1. It's being run in a folder where a CASTEP calculation has just completed.
#   2. None of the resulting CASTEP outputs have been deleted.
#   3. The CASTEP task was ELNES so that .eels_mat files exist.
#   4. The executable "nexspec" is in your $PATH.
#
# The usage of the script is as follows:
# 
#   castep_pp_nexafs.py PROJECTOR MINBAND
#
# Here PROJECTOR is the (1-based) index of the core-level projector of interest (i.e. 
# the initial state of the NEXAFS transition to be considered). MINBAND is the minimum
# band index (1-based) to be included in the Fermi's Golden Rule sum. Both of these
# are passed onto nexspec.

# Tried to use modern interfaces here...
import glob
import subprocess
import argparse

# Change these settings as necessary
nexspec = "nexspec"

parser = argparse.ArgumentParser(description="Post-process CASTEP for NEXAFS outut.")
parser.add_argument('projector', help="1-based index of the core-level projector of interest (initial state).")
parser.add_argument('minband', help="1-based band index of the minimum band to include in the Fermi's Golden Rule sum.")
args = parser.parse_args()

# Remove the energies file (since we want to repeatedly append to it) if it exists.
subprocess.call(["rm", "energies.txt"])
# and again for the seeds file.
subprocess.call(["rm", "seeds.txt"])

# Get our seedlist
files = glob.glob("*.castep")

for f in files:
  
  # Get the seed and run nexspec
  seed = f.split(".castep")[0]
  print "Working with seed = %s" % seed
  subprocess.call("%s %s %s %s" % (nexspec, seed, args.projector, args.minband), shell=True)
  
  # Open the castep file associated with this seed and grab the species and relevant
  # energies
  c = open(seed+".castep", 'r')
  clines = c.readlines()
  c.close()
  ae_species = []
  ps_species = []
  elec_config = []
  ae_energies = []
  ps_energies = []
  total_energy = 0.0
  for i, line in enumerate(clines):
    if "Atomic calculation performed for" in line:
      # 5th word is the species with ':' appended, so remove last character.
      ae_species.append(line.split()[4][:-1])
      elec_config.append(" ".join(line.split()[5:]))
      # Two lines further on is the ae energy.
      ae_energies.append(float(clines[i+2].split()[9]))
    elif "Pseudo atomic calculation performed" in line:
      ps_species.append(line.split()[5])
      ps_energies.append(float(clines[i+2].split()[9]))
    elif "Final energy, E" in line:
      total_energy = float(line.split()[4])
  
  # Make lists unique by key
  ae_dict = {}
  config_dict = {}
  ps_dict = {}
  for (si, ei) in zip(ae_species, ae_energies):
    ae_dict[si] = ei
  for (si, ci) in zip(ae_species, elec_config):
    config_dict[si] = ci
  for (si, ei) in zip(ps_species, ps_energies):
    ps_dict[si] = ei
  
  # Now write/append the energies. Note that for the species, we overwrite assuming that
  # all are the same, for the energy we append instead. We also write the seed to a 
  # file so we know what order this stuff was processed in.
  ef = open("energies.txt", 'a')
  sf = open("species.txt", 'w')
  cf = open("configurations.txt", 'w')
  ff = open("seeds.txt", 'a')
  
  for k in ae_dict.keys():
    sf.write("%s %g %g\n" % (k, ae_dict[k], ps_dict[k]))
  
  for k in ae_dict.keys():
    cf.write("%s %s\n" %(k, config_dict[k]))
  
  ff.write("%s\n" % seed)
  ef.write("%.15g\n" % total_energy)
  
  ff.close()
  ef.close()
  sf.close()
  cf.close()