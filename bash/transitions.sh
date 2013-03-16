#!/bin/bash

# Extracts the total energies from each of the core hole calculations and puts them
# in a single file, does the same for the neutral calculation, then runs the python
# script that calculates the transition energies (written to transitions.txt). 
# Inputs are the CASTEP seed name and the name of a file containing the unique identifiers
# for each of the core hole calculations, so that the core hole seeds are constructed
# as ${seed}_${identifier} for each identifier. Normally I just give the identifiers numbers
# and the identifier file "atoms.list". All the ${seed}_${identifier}.castep files must
# be in the directory you run the script.

# In addition, must have a ${seed}_neutral.castep file available in the same directory.

# The python script requires that a file called atom_energies.txt be present that
# has four energies in this order: the AE energies of the neutral and core hole atom,
# and the PS energies of the neutral and core hole atom. This is the order they
# appear in the .castep header.

seed=$1
atlist=$2

if [ -e total_energies.txt ]
then
  rm total_energies.txt
fi

for i in `cat $atlist`
do
  cat ${seed}_${i}.castep | grep "Final energy," | awk '{print $5}' >> total_energies.txt
done

if [ -e ${seed}_neutral.castep ]
then
  if [ -e neutral.txt ]
  then
    rm neutral.txt
  fi

  cat ${seed}_neutral.castep | grep "Final energy," | awk '{print $5}' >> neutral.txt
else
  echo "Missing neutral castep file: See the transitions.sh file for details. Cannot continue."
  exit
fi

if [ -e atom_energies.txt ]
then
  python `dirname $0`/calculate_transitions.py total_energies.txt neutral.txt atom_energies.txt
else
  echo "atom_energies.txt does not exist - can't run the python script. See transitions.sh file for explanation."
  exit
fi
