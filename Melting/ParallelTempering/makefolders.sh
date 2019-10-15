#!/bin/bash

path='/home/kiran/Melting/BigCrystal/SingleRise/T12/21/ParallelTempering/'
#hcounter=12
#while [ $hcounter -le 14 ];do
    runcounter=1
    while [ $runcounter -le 10 ]; do
	#mkdir $path$'T'$hcounter$'/'$'Run'$runcounter
        mkdir $path$'Run'$runcounter
	let runcounter=runcounter+1
    done
    #let hcounter=hcounter+2
#done
