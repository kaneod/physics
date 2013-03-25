
import sys
from numpy import loadtxt, array, savetxt

total_energies = loadtxt(sys.argv[1])
neutral_energy = loadtxt(sys.argv[2])
atom_energies = loadtxt(sys.argv[3])

# Because total_energies can be a single line, have to explicitly check here.
if total_energies.shape is ():
  total_energies = array([total_energies])

Eae = atom_energies[0]
Eaestar = atom_energies[1]
Eps = atom_energies[2]
Epsstar = atom_energies[3]

Et = neutral_energy[0]

dEaea = Eaestar - Eae
dEpsa = Epsstar - Eps
dEca = dEaea - dEpsa

print "dEca = %f" % (dEca)

transitions = []

for Etstar in total_energies:
  dEv = Etstar - Et
  print "dEv = %f" % dEv
  transitions.append(dEv + dEca)

transitions=array(transitions)

savetxt("transition_energies.txt", transitions, fmt='%.6f')
