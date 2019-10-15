#!/bin/bash

path='/home/kiran/700_Noeq/'
runcounter=16
while [ $runcounter -le 45 ]; do
    mkdir $path$'Run'$runcounter
    let runcounter=runcounter+1
done
