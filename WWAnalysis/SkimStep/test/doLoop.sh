#!/usr/bin/env bash

wwTag=`showtags | grep WWAnalysis | head -n 1 | awk '{print $1}'`

cat $1 | while read x; do

    [[ $x == ""   ]] && continue
    [[ "${x:0:1}" == "#" ]] && continue

    sampleNum=`echo $x | c 1`
    sampleName=`echo $x | c 2`
    dbsname=`echo $x | c 3`
    whatmode=`echo $x | c 4`
    tag=`echo $x | c 5`
    numEvents=`echo $x | c 6`
    evtsPer=20000
    fakes=None
    dbsnamenew=`echo $dbsname    | tr "/" "@"`
    echo dbsnamenew=$dbsnamenew
    dbsnamenew=`echo $dbsnamenew    | sed -e s/@/?$/g` 
    echo dbsnamenew=$dbsnamenew
    dbsnamenew=`echo $dbsnamenew | tr "?" "\\\"`
    echo dbsnamenew=$dbsnamenew
    dbsnamenew=`echo $dbsnamenew | tr "$" "/"` 
    echo dbsnamenew=$dbsnamenew

    echo sampleNum=$sampleNum
    echo sampleName=$sampleName
    echo dbsname=$dbsname
    echo whatmode=$whatmode
    echo tag=$tag
    echo numEvents=$numEvents
    echo evtsPer=$evtsPer
    echo fakes=$fakes
    echo dbsnamenew=$dbsnamenew
   


    mkdir -p $sampleNum.$sampleName

    cat latinosYieldSkim.py | \
        sed s/WHATMODE/$whatmode/g |  \
        sed s/RMMEFN/$sampleName.root/g | \
        sed s/RMMEFAKE/$fakes/g |  \
        sed s/RMMEGlobalTag/$tag/g > $sampleNum.$sampleName/$sampleName.py

    [ $whatmode = "DATA" ] && file=crabData.cfg || file=crabMC.cfg
    echo $sampleName
    cat $file | \
        sed s/RMMEDATASET/$dbsnamenew/g | \
        sed s/RMMENUM/$sampleNum/g | \
        sed s/RMMENOEVENTS/$numEvents/g | \
        sed s/RMMEEVTSPERJOB/$evtsPer/g | \
        sed s/RMMEPRELABEL/$wwTag/g | \
        sed s/RMMEFN/$sampleName/g > $sampleNum.$sampleName/crab.cfg

done
    
