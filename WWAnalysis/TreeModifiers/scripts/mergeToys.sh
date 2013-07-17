#!/bin/bash

dirs='results/c1DNoMassErr results/c1DMassErr results/c2DMassErr'
channels='4mu 4e 2e2mu comb'

for dir in $dirs
  do
    cd $dir
    pwd=$PWD
    echo $pwd
    for channel in $channels 
      do
      echo "   merging toys for channel HZZ: " $channel
      hadd -k -f toys_$channel\_new\.root higgsCombineMASS_$channel.MultiDimFit.mH125.8*.root  
    done
    cd -
done

echo "done."
