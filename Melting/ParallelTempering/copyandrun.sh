#!/bin/bash

frompath='/home/kiran/Melting/BigCrystal/SingleRise/T12/21/ParallelTempering/'
topath='/home/kiran/Melting/BigCrystal/SingleRise/T12/21/ParallelTempering/'
file='run.sh'
file2='*.h'
file3='LD_woutput.cpp'
counter=1
Tcounter=12

#while [ $Tcounter -le 20 ]; do
#    counter=51
    while [ $counter -le 10 ]; do
	copytofolder=$topath$'/Run'$counter$'/'
	cp $frompath$file $copytofolder
	cd $copytofolder
	#tail -n 7020000 coordinates_2.xyz >> coordinates.xyz
	#tail -n 702 coordinates.xyz > last.xyz
	sbatch $file
	let counter=counter+1
    done
#    let Tcounter=Tcounter+2
#done
