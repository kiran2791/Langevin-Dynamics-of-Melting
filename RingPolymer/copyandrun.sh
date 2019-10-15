#!/bin/bash

frompath='/home/kiran/700_Noeq/'
topath='/home/kiran/700_Noeq/'
file='run.sh'
file2='*.h'
file3='LD_woutput.cpp'
counter=1
Hcounter=5

    copytofolder=$topath
    while [ $counter -le 45 ]; do
        cp $frompath$file $copytofolder$'Run'$counter$'/'
        cd $copytofolder$'Run'$counter$'/'
	#tail -n 702 coordinates.xyz > last.xyz
        sbatch $file
	#cp $frompath$file2 $copytofolder$'Run'$counter$'/'
    #cp $frompath$file3 $copytofolder$'Run'$counter$'/'
        let counter=counter+1
    done
