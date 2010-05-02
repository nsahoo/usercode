#!/bin/bash

### DOCUMENTATION ###
# task: 
#  --> resid
#  --> resol
#
# action:
#  --> submit
#  --> copy
#  --> macro
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
	#attempt=2
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
	cat $DIR/test/$pset.template.py |sed "s#SET_GLOBALTAG#$globalTag#" | \
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
	hadd merged.$outputFile $listFiles >& /dev/null 
	rm $listFiles
	cd -
	;;
    
    macro) echo "run macro"
	if [ ! -d outputMacro/$taskName ] 
	then mkdir -p outputMacro/$taskName 
	fi
	
	####### settings of the root macro project.C ########
        # prefix(1-2),category(1-3),type(1-4),projection(1-3)
	#####################################################

	input=$PWD/output/$taskName/merged.$outputFile
	echo "input: " $input
	cd outputMacro/$taskName
	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#200# | sed s#SET_PLOW#-2000# | sed s#SET_PHIGH#2000# | \
	    sed s#SETBINS#60#	> project.C
	echo "producing projections of rawD0 vs pt..."
	root.exe -b -q "project.C+(1,1,1,1)"

	echo "producing projections of rawDz vs pt..."
	root.exe -b -q "project.C+(1,1,2,1)"

	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#300# | sed s#SET_PLOW#-3000# | sed s#SET_PHIGH#3000# | \
	    sed s#SETBINS#50# | sed s#SETPTMIN#1.0# | sed s#SETPTMAX#1.4#   \
	    > project.C
	
	echo "producing projections of rawD0 vs eta..."
	root.exe -b -q "project.C+(1,1,1,2)"

	echo "producing projections of rawDz vs eta..."
	root.exe -b -q "project.C+(1,1,2,2)"
	
	cd -
	;;

    macro2) echo "run macro"
	if [ ! -d outputMacro/$taskName ] 
	then mkdir -p outputMacro/$taskName 
	fi
	
	####### settings of the root macro project.C ########
        # prefix(1-2),category(1-3),type(1-2),projection(1-3)
	#####################################################

	input=$PWD/output/$taskName/merged.$outputFile
	echo "input: " $input
	cd outputMacro/$taskName
	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#200# | sed s#SET_PLOW#-2000# | sed s#SET_PHIGH#2000# | \
	    sed s#SETBINS#60#	> project.C
	echo "producing projections of resoD0 vs pt..."
	root.exe -b -q "project.C+(1,2,3,1)"

	echo "producing projections of resoDz vs pt..."
	root.exe -b -q "project.C+(1,2,4,1)"

	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#300# | sed s#SET_PLOW#-3000# | sed s#SET_PHIGH#3000# | \
	    sed s#SETBINS#50# | sed s#SETPTMIN#1.0# | sed s#SETPTMAX#1.4#   \
	    > project.C
	
	echo "producing projections of resoD0 vs eta..."
	root.exe -b -q "project.C+(1,2,3,2)"

	echo "producing projections of resoDz vs eta..."
	root.exe -b -q "project.C+(1,2,4,2)"
	
	cd -
	;;

    macro3) echo "run macro"
	if [ ! -d outputMacro/$taskName ] 
	then mkdir -p outputMacro/$taskName 
	fi
	
	####### settings of the root macro project.C ########
        # prefix(1-2),category(1-3),type(1-2),projection(1-3)
	#####################################################

	input=$PWD/output/$taskName/merged.$outputFile
	echo "input: " $input
	cd outputMacro/$taskName
	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#200# | sed s#SET_PLOW#-2000# | sed s#SET_PHIGH#2000# | \
	    sed s#SETBINS#60#	> project.C
	echo "producing projections of respD0 vs pt..."
	root.exe -b -q "project.C+(1,3,3,1)"

	echo "producing projections of respDz vs pt..."
	root.exe -b -q "project.C+(1,3,4,1)"

	sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
	    sed s#SET_PBINS#300# | sed s#SET_PLOW#-3000# | sed s#SET_PHIGH#3000# | \
	    sed s#SETBINS#50# | sed s#SETPTMIN#1.0# | sed s#SETPTMAX#1.4#   \
	    > project.C
	
	echo "producing projections of respD0 vs eta..."
	root.exe -b -q "project.C+(1,3,3,2)"

	echo "producing projections of respDz vs eta..."
	root.exe -b -q "project.C+(1,3,4,2)"
	
	cd -
	;;


    *) echo "ERROR: input type of action is not recognized. Choose between \"submit\" and \"copy\" "
	exit 1
	;;
esac


