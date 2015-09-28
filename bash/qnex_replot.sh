
# May need to customize this.
XS="mpirun -np 2 /opt/bin/xspectra.x -nk 2"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

echo "################# REPLOT XSPECTRA #####################"
echo "#"
echo "# Remember: atomic indices are the positional args."
echo "#"
echo "#######################################################"
echo ""

# If a replot file already exists, use it, otherwise generate one.
if [ -e replot.xspectra ]; then
	echo "Detected existing replot.xspectra - using."
else
	echo "No replot.xspectra detected - I will generate one for you."
	echo "Modify it as necessary and re-run this script."
	cat > replot.xspectra << EOF
&input_xspectra
  calculation='xanes_dipole',
  xonly_plot = .true.,
  x_save_file = '<<save>>',
  outdir='./',
  xniter=2000,
/
&plot
	gamma_mode='variable',
  xnepoint=300,
  gamma_energy(1) = 8.0,
  gamma_energy(2) = 28.0,
  gamma_value(1) = 0.3,
  gamma_value(2) = 2.0,
  !xgamma=0.05,
  xemin=-10.0,
  xemax=35.0,
  terminator=.true.,
  cut_occ_states=.true.,
/
&pseudos
  filecore='NotApplicable',
  r_paw(1)=3.2,
/
&cut_occ
!cut_desmooth = 0.1
/
2 2 1 0 0 0
EOF
	exit
fi

for i in $@; do
	echo "---------------- ATOM $i --------------------"
	echo "XSpectra calculations..."
	cp replot.xspectra base.in
	sed -i '' "s/<<atom>>/atom$i/" base.in
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/1/" tmp.in
	sed -i '' "s/<<ey>>/0/" tmp.in
	sed -i '' "s/<<ez>>/0/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.100.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.100.out
	mv xanes.dat xanes.atom${i}.100.dat	
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/1/" tmp.in
	sed -i '' "s/<<ey>>/1/" tmp.in
	sed -i '' "s/<<ez>>/0/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.110.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.110.out
	mv xanes.dat xanes.atom${i}.110.dat			
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/0/" tmp.in
	sed -i '' "s/<<ey>>/1/" tmp.in
	sed -i '' "s/<<ez>>/0/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.010.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.010.out
	mv xanes.dat xanes.atom${i}.010.dat		
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/0/" tmp.in
	sed -i '' "s/<<ey>>/1/" tmp.in
	sed -i '' "s/<<ez>>/1/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.011.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.011.out
	mv xanes.dat xanes.atom${i}.011.dat	
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/0/" tmp.in
	sed -i '' "s/<<ey>>/0/" tmp.in
	sed -i '' "s/<<ez>>/1/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.001.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.001.out
	mv xanes.dat xanes.atom${i}.001.dat
	cp base.in tmp.in
	sed -i '' "s/<<ex>>/1/" tmp.in
	sed -i '' "s/<<ey>>/0/" tmp.in
	sed -i '' "s/<<ez>>/1/" tmp.in
	sed -i '' "s/<<save>>/xspectra.atom${i}.101.sav/" tmp.in
	$XS -i tmp.in > xspectra.replot.atom${i}.101.out
	mv xanes.dat xanes.atom${i}.101.dat
done