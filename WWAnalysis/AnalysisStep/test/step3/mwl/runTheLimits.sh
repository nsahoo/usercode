#!/usr/bin/env bash

prefix=../inputCards/scenario4/hww-2.13fb
nParrallel=7
lumi=2.1
niceLumi=2.1fb
label=testLP # or maybe S1 for scenario 1?

mkdir -p out plots root
pids=
for mass in $masses; do 
    combine -M Asymptotic -v 2 -m $mass -n _${label}_comb $prefix.mH$mass.comb.txt 2> /dev/null | tee out/$mass.comb.out > /dev/null &
    pids="$! $pids"
    #uncomment me when you have the 2j cards
    #for j in _0j _1j _2j; do
    for j in 0j 1j 2j; do
        combine -M Asymptotic -v 2 -m $mass -n _${label}_comb_$j $prefix.mH$mass.comb_$j.txt 2> /dev/null | tee out/$mass.comb_$j.out > /dev/null &
        pids="$! $pids"
        while [ `ps | grep combine | wc -l` -ge $nParrallel ] ; do sleep 10; done
    done
done
wait $pids
mv *.root root/

rm limitSummary_${label}_comb.txt
for mass in $masses; do
    a=`grep "Limit:"                       out/$mass.comb.out | awk '{print $5}'`
    b=`grep "Median for expected limits:"  out/$mass.comb.out | awk '{print $5}'`
    c=`grep "Expected  2.5%:"              out/$mass.comb.out | awk '{print $5}'`
    d=`grep "Expected 16.0%:"              out/$mass.comb.out | awk '{print $5}'`
    e=`grep "Expected 84.0%:"              out/$mass.comb.out | awk '{print $5}'`
    f=`grep "Expected 97.5%:"              out/$mass.comb.out | awk '{print $5}'`
    echo $mass $a 99. $b $c $d $e $f >> limitSummary_${label}_comb.txt
done


for j in 0j 1j 2j; do
#for j in 0j 1j; do
    rm limitSummary_${label}_comb_$j.txt
    for mass in $masses; do
        a=`grep "Limit:"                       out/$mass.comb_$j.out | awk '{print $5}'`
        b=`grep "Median for expected limits:"  out/$mass.comb_$j.out | awk '{print $5}'`
        c=`grep "Expected  2.5%:"              out/$mass.comb_$j.out | awk '{print $5}'`
        d=`grep "Expected 16.0%:"              out/$mass.comb_$j.out | awk '{print $5}'`
        e=`grep "Expected 84.0%:"              out/$mass.comb_$j.out | awk '{print $5}'`
        f=`grep "Expected 97.5%:"              out/$mass.comb_$j.out | awk '{print $5}'`
        echo $mass $a 99. $b $c $d $e $f >> limitSummary_${label}_comb_$j.txt
    done
done;

                                                                                                          #don't peak!
root -b -l -q PlotLimit.C+"(\"limitSummary_${label}_comb.txt\",   \"limit_${niceLumi}_${label}_comb\",\"Hww, combined, L = $lumi fb^{-1}\",false,true)"
#for j in 0j 1j; do
for j in 0j 1j 2j; do
    root -b -l -q PlotLimit.C+"(\"limitSummary_${label}_comb_$j.txt\",\"limit_${niceLumi}_${label}_comb_$j\",\"Hww, 0-jet, L = $lumi fb^{-1}\",   false,true)"
done




