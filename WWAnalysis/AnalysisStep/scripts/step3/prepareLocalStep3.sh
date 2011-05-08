#!/usr/bin/env bash

wwTag=`showtags | grep WWAnalysis | head -n 1 | awk '{print $1}'`
jsonFile="$CMSSW_BASE/src/WWAnalysis/Misc/Jsons/fullJson.json"
mass=${1:-160}
isLocal=${2:-False}
cat loopOverAll | while read x; do
    echo $x | grep "10[1-6][0-9][0-9][0-9]" > /dev/null 2>&1 
    if [[ $? -eq 0 ]] ; then
        echo $x | grep "10[1-6]$mass" > /dev/null 2>&1 || continue
    fi
    sampleNum=`echo $x | c 1`
    sampleName=`echo $x | c 2`
    numEvents=`echo $x | c 3`
    numJobs=`echo $x | c 4`
    isMC=`echo $x | c 5`
    tag=`echo $x | c 6`
    sampleType=`echo $x | c 7`
    dbsName=`echo $x | c 8`
    if [[ "$isLocal" == "True" ]]; then
	echo "quering dbs for sample " $sampleName
	inputFilesTmp=`dbsql2 "find file where dataset like $dbsName"|grep store | awk '{print "\""$1"\","}'`
	inputFiles=`echo $inputFilesTmp|tr -d " "`
    fi

    mkdir -p m$mass/$sampleNum.$sampleName

    cat createEventHistsFullFW_cfg.py | \
        sed s/RMMEFN/$sampleName.root/g |  \
        sed s/RMMEMC/$isMC/g |  \
        sed s/RMMENUM/$sampleNum/g | \
        sed s/RMMENAME/$sampleName/g | \
        sed s/RMMEMASS/$mass/g | \
        sed s/RMMESAMPLE/$sampleType/g | \
	sed s#RMMEINPUTFILES#$inputFiles#g | \
	sed s/RMMEISLOCAL/$isLocal/g | \
        sed s/RMMEGlobalTag/$tag/g > m$mass/$sampleNum.$sampleName/$sampleName.py

done
