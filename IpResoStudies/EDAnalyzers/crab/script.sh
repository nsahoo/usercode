#!/bin/bash
DIR="${CMSSW_BASE}/src/IpResoStudies/EDAnalyzers/"
#DIR=`$CMSSW_BASE/src/IpResoStudies/EDAnalyzers/`

echo "DIR: " $DIR

choice=$1



case $choice in
    0) echo "choice is 0. Setup MC crab jobs"
        ### for MC 
	attempt=4
	dataset=/MinBias/Spring10-START3X_V26A_356ReReco-v1/GEN-SIM-RECO 
	globalTag='START3X_V26A::All'
	taskName=RawResidualsMC-v$attempt
	outputFile="rawResidualsMC.root"
	;;

    1) echo "choice is 1"
	;;

    *) echo "default choice"
	;;
esac


cat $DIR/test/residuals.template.py |sed "s#SET_GLOBALTAG#$globalTag#" |sed "s#SET_OUTPUT#$outputFile#" > pset.$taskName.py
cat $DIR/crab/crab.template.cfg |sed "s#SET_DATASET#$dataset#" | sed "s#SET_TASK_NAME#$taskName#" > tmp.cfg
cat tmp.cfg| sed "s#SET_PSET#$PWD/pset.$taskName.py#" > crab.$taskName.cfg; rm tmp.cfg
crab -create -submit -cfg crab.$taskName.cfg