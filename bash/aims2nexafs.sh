#!/bin/bash

# Usage: aims2nexafs.sh INPUT

# Need 3 pw.x inputs and an xspectra input.

# First is the ground state calculation
esc_switch.py -p "ground" -i aims -o qe $1 ${1%.*}.ground.qe

# Second is the FCH calculation
esc_switch.py -p "<<atom>>" -x 1.0 -i aims -o qe $1 ${1%.*}.fch.qe

#Third is the TP calculation
esc_switch.py -p "<<atom>>" -x 0.5 -i aims -o qe $1 ${1%.*}.tp.qe

# The XSpectra output is mostly the same for all calculations:
cat > ${1%.*}.xspectra.qe << EOF
&input_xspectra
    calculation='xanes_dipole',
    prefix='<<atom>>',
    outdir='./',
    xniter=2000,
    xcheck_conv=50,
    xepsilon(1)=<<ex>>,
    xepsilon(2)=<<ey>>,
    xepsilon(3)=<<ez>>,
    xcoordcrys=.false.,
    xiabs=1,
    x_save_file='xspectra.sav',
    xerror=0.02,
    xe0 = <<homo>>
/
&plot
		gamma_mode='variable',
    xnepoint=300,
  	gamma_energy(1) = 5.0,
  	gamma_energy(2) = 28.0,
  	gamma_value(1) = 0.25,
  	gamma_value(2) = 2.0,
    !xgamma=0.05,
    xemin=-10.0,
    xemax=35.0,
    terminator=.true.,
    cut_occ_states=.true.,
/
&pseudos
    filecore='<<corewf>>',
    r_paw(1)=3.2,
/
&cut_occ
  !cut_desmooth = 0.1
/
X X X 0 0 0
EOF


