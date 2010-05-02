#!/bin/bash

### DOCUMENTATION ###
# task: 
#  --> resid
#  --> resol
#
# action:
#  --> submit
#  --> copy
######################

DIR="${CMSSW_BASE}/src/IpResoStudies/EDAnalyzers/"

#echo "DIR: " $DIR

task=$1
action=$2

case $task in
    resid)  ### preparing setup for MC RawResiduals tasks
	attempt=4
	dataset=/MinBias/Spring10-START3X_V26A_356ReReco-v1/GEN-SIM-RECO 
	globalTag='START3X_V26A::All'
	taskName=RawResidualsMC-v$attempt
	pset=residuals
	outputFile="rawResidualsMC.root"
	;;
    resol)  ### preparing setup for MC response and resolution tasks
	attempt=1
	dataset=/MinBias/Spring10-START3X_V26A_356ReReco-v1/GEN-SIM-RECO 
	globalTag='START3X_V26A::All'
	taskName=ResponseAndResolutionsMC-v$attempt
	pset=respAndResoOnMC
	outputFile="respAndResoMC.root"
	;;
    *) echo "ERROR: input type of task is not recognized. Choose between \"resid\" and \"resol\" "
	exit 1
	;;
esac


case $action in
    submit) echo "action is submit crab jobs"
	cat $DIR/test/residuals.template.py |sed "s#SET_GLOBALTAG#$globalTag#" | \
	    sed "s#SET_OUTPUT#$outputFile#" > pset.$taskName.py
	cat $DIR/crab/crab.template.cfg |sed "s#SET_DATASET#$dataset#" | sed "s#SET_TASK_NAME#$taskName#" > tmp.cfg
	cat tmp.cfg| sed "s#SET_PSET#$PWD/pset.$taskName.py#" > crab.$taskName.cfg; rm tmp.cfg
	crab -create -cfg crab.$taskName.cfg
	rm crab.$taskName.cfg; rm pset.$taskName.py;
	crab -submit -c $taskName
	;;
    

    copy) echo "action is copy: copying file from SE and merge them"
	hadoopPath="/hadoop/cms/store/user/mangano/IPStudies2010/mangano/MinBias/$taskName"
	if [ ! -d output/$taskName ] 
	then mkdir -p output/$taskName 
	fi
	scp   bmangano@uaf-6.t2.ucsd.edu:$hadoopPath/*/*.root output/$taskName
	cd output/$taskName
	listFiles=$(ls *.root)
	echo "list: " $listFiles
	hadd merged.root $listFiles >& /dev/null 
	rm $listFiles
	cd -
	;;


    *) echo "ERROR: input type of action is not recognized. Choose between \"submit\" and \"copy\" "
	exit 1
	;;
esac


