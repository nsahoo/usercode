#!/bin/bash

### DOCUMENTATION ###
# task: 
#  --> dataRaw
#  --> resid
#  --> resol
#
# action:
#  --> submit
#  --> copy
#  --> newMacro
#  --> allMacroRaw
#  --> allMacroResol
######################

DIR="${CMSSW_BASE}/src/IpResoStudies/EDAnalyzers/"

#echo "DIR: " $DIR

task=$1
action=$2

case $task in
    dataRaw)  ### preparing setup for DATA RawResiduals tasks
	attempt=1
	dataset=/MinimumBias/Commissioning10-May6thReReco-v1/RECO 
	globalTag='GR_R_35X_V8B::All'
	taskName=RawResidualsDATA-v$attempt
	pset=residuals
	outputFile="rawResidualsDATA.root"
	;;

    resid)  ### preparing setup for MC RawResiduals tasks
	attempt=5
	dataset=/MinBias/Spring10-START3X_V26A_356ReReco-v1/GEN-SIM-RECO 
	globalTag='START3X_V26A::All'
	taskName=RawResidualsMC-v$attempt
	pset=residuals
	outputFile="rawResidualsMC.root"
	;;
    resol)  ### preparing setup for MC response and resolution tasks
	attempt=4
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
	cat $DIR/crab/crab.template.cfg |sed "s#SET_DATASET#$dataset#" | \
	    sed "s#SET_TASK_NAME#$taskName#" > tmp.cfg
	cat tmp.cfg| sed "s#SET_PSET#$PWD/pset.$taskName.py#" > crab.$taskName.cfg; rm tmp.cfg
	crab -create -cfg crab.$taskName.cfg
	rm crab.$taskName.cfg; rm pset.$taskName.py;
	crab -submit -c $taskName
	;;
    

    copy) echo "action is copy: copying file from SE and merge them"
	#hadoopPath="/hadoop/cms/store/user/mangano/IPStudies2010/mangano/MinBias/$taskName"
	hadoopPath="/hadoop/cms/store/user/mangano/IPStudies2010/mangano/MinimumBias/$taskName"
	if [ ! -d output/$taskName ] 
	then mkdir -p output/$taskName 
	fi
	scp   bmangano@uaf-6.t2.ucsd.edu:$hadoopPath/*/*.root output/$taskName
	cd output/$taskName
	listFiles=$(ls *.root)
	echo "list: " $listFiles
	hadd merged.$outputFile $listFiles >& /dev/null 
	#rm $listFiles
	cd -
	;;
    

    newMacro) echo "run macro"
	if [ ! -d outputMacro/$taskName ] 
	then mkdir -p outputMacro/$taskName 
	fi
	
	####### settings of the root macro project.C ########
        # prefix(1-2),category(1-3),type(1-6),projection(1-3)
	#####################################################
	prefix=$3
	category=$4
	type=$5
	projection=$6

	input=$PWD/output/$taskName/merged.$outputFile
	echo "input: " $input
	cd outputMacro/$taskName

	if(( $projection == 1 )); then 
	    if(( $category == 1 || $category == 2 )); then
		sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
		    sed s#SET_PBINS#400# | sed s#SET_PLOW#-2000# | sed s#SET_PHIGH#2000# | \
		    sed s#SETBINS#60# |	\
		    sed s/project\(/project_${prefix}_${category}_${type}_${projection}\(/	\
		    > project_${prefix}_${category}_${type}_${projection}.C
			root.exe -b -q "project_${prefix}_${category}_${type}_${projection}.C+($prefix,$category,$type,$projection)"
	    elif(( $category == 3)); then
		sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
		    sed s#SET_PBINS#400# | sed s#SET_PLOW#-200# | sed s#SET_PHIGH#200# | \
		    sed s#SETBINS#60# |	\
		    sed s/project\(/project_${prefix}_${category}_${type}_${projection}\(/	\
		    > project_${prefix}_${category}_${type}_${projection}.C
			root.exe -b -q "project_${prefix}_${category}_${type}_${projection}.C+($prefix,$category,$type,$projection)"
	    fi	
	elif(( $projection == 2 )); then
	    if(( $category == 1 || $category == 2 )); then
		sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
		    sed s#SET_PBINS#600# | sed s#SET_PLOW#-3000# | sed s#SET_PHIGH#3000# | \
		    sed s#SETBINS#50# | sed s#SETPTMIN#1.0# | sed s#SETPTMAX#1.4# |   \
		    sed s/project\(/project_${prefix}_${category}_${type}_${projection}\(/	\
		    > project_${prefix}_${category}_${type}_${projection}.C
			root.exe -b -q "project_${prefix}_${category}_${type}_${projection}.C+($prefix,$category,$type,$projection)"
	    elif(( $category == 3)); then
		sed s#SET_INPUTSIM#$input# $DIR/macro/project.C | \
		    sed s#SET_PBINS#600# | sed s#SET_PLOW#-300# | sed s#SET_PHIGH#300# | \
		    sed s#SETBINS#50# | sed s#SETPTMIN#1.0# | sed s#SETPTMAX#1.4# |   \
		    sed s/project\(/project_${prefix}_${category}_${type}_${projection}\(/	\
		    > project_${prefix}_${category}_${type}_${projection}.C
			root.exe -b -q "project_${prefix}_${category}_${type}_${projection}.C+($prefix,$category,$type,$projection)"
	    fi
	fi

	cd -
	;;


    allMacroRaw)
	##./script.sh resid newMacro 1 1 1 1 >& log1 & #sim raw d0 pt
	##./script.sh resid newMacro 1 1 1 2 >& log2 & #sim raw d0 eta
	#./script.sh resid newMacro 1 1 1 3 >& log3 & #sim raw d0 phi
	##./script.sh resid newMacro 1 1 2 1 >& log4 & #sim raw dz pt
	##./script.sh resid newMacro 1 1 2 2 >& log5 & #sim raw dz eta
	#./script.sh resid newMacro 1 1 2 3 >& log6 & #sim raw dz phi

	./script.sh dataRaw newMacro 2 1 1 1 >& log7 & #data raw d0 pt
	./script.sh dataRaw newMacro 2 1 1 2 >& log8 & #data raw d0 eta
	#./script.sh dataRaw newMacro 2 1 1 3 >& log9 & #data raw d0 phi
	./script.sh dataRaw newMacro 2 1 2 1 >& log10 & #data raw dz pt
	./script.sh dataRaw newMacro 2 1 2 2 >& log11 & #data raw dz eta
	#./script.sh dataRaw newMacro 2 1 2 3 >& log12 & #data raw dz phi

	;;


    allMacroResol)
	./script.sh resol newMacro 1 2 3 1 >& log1 & #reso d0 pt     
	./script.sh resol newMacro 1 2 3 2 >& log2 & #reso d0 eta    
	./script.sh resol newMacro 1 2 4 1 >& log3 & #reso dz pt     
	./script.sh resol newMacro 1 2 4 2 >& log4 & #reso dz eta    

	./script.sh resol newMacro 1 3 5 1 >& log5 & #resp d0 pt      
	./script.sh resol newMacro 1 3 5 2 >& log6 & #resp d0 eta   
	./script.sh resol newMacro 1 3 6 1 >& log7 & #resp dz pt    
	./script.sh resol newMacro 1 3 6 2 >& log8 & #resp dz eta   
	;;



    *) echo "ERROR: input type of action is not recognized. Choose between \"submit\" and \"copy\" "
	exit 1
	;;
esac


