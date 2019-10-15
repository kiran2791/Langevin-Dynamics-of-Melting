#!/bin/bash

frompath='/home/kiran/Melting/BigCrystal/SingleRise/T12/21/ParallelTempering/'
topath='/home/kiran/Melting/BigCrystal/SingleRise/T12/21/ParallelTempering/'
file='*.cpp'
file2='*.h'
file3='*.xyz'
#file4='*.dat'
counter=1
Hcounter=12

while [ $counter -le 10 ]; do
    #copytofolder=$topath$'T'$Hcounter$'/'
    #mkdir $copytofolder$'OrderParameter'
    #cd $copytofolder$'OrderParameter'
    #gfortran OPAveraging.f90
    #./a.out
    #counter=1
    #while [ $counter -le 30 ]; do
    cp $copytofolder$file $copytofolder$'Run'$counter$'/'
    cp $copytofolder$file2 $copytofolder$'Run'$counter$'/'
    cp $copytofolder$file3 $copytofolder$'Run'$counter$'/'
	#cp $copytofolder$'Run'$counter'/'$file $copytofolder$'OrderParameter/orderparameters_'$counter$'.dat'
    let counter=counter+1
    #done
    #let Hcounter=Hcounter+2
done
