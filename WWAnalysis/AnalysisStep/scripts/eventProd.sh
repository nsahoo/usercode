#!/bin/bash


### DOCUMENTATION ###
#
# action:
#  --> submit
######################

action=$1

DIR="${CMSSW_BASE}/src/WWAnalysis/AnalysisStep"

label=oct21
dataset=boh

case $2 in
    ### Wjets:
    1)
	dataset=/Wmunu/mangano-SkimStep005-bkg1-oct15-2455c73529abf8a64e80894749518608/USER
	taskName=EventProd006-bkg1-$label
	;;
    2)
	dataset=/Wenu/mangano-SkimStep005-bkg2-oct15-2455c73529abf8a64e80894749518608/USER
	taskName=EventProd006-bkg2-$label
	;;
    3)
	dataset=/Wtaunu/mangano-SkimStep005-bkg3-oct15-2455c73529abf8a64e80894749518608/USER
	taskName=EventProd006-bkg3-$label
	;;
    4)
	dataset=/Zgamma_2l_enriched/mangano-SkimStep005-bkg4-v9-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg4-$label
	;;
    5)
	dataset=/WW_2l_7TeV/mangano-SkimStep005-bkg5-v9-2455c73529abf8a64e80894749518608/USER
	taskName=EventProd006-bkg5-$label
	;;
    6)
	dataset=/gg_WW_2l2nu-GG2WW/mangano-SkimStep005-bkg6-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg6-$label
	;;
    7)
	dataset=/WZ/mangano-SkimStep005-bkg7-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg7-$label
	;;
    8)
	dataset=/ZZ/mangano-SkimStep005-bkg8-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg8-$label
	;;
#    9)
#	dataset=Wgammma
#	taskName=EventProd006-bkg9-$label
#	;;
#    10)
#	dataset=Zgamma
#	taskName=EventProd006-bkg10-$label
#	;;
    11)
	dataset=/TTbar_2l/mangano-SkimStep005-bkg11-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg11-$label
	;;
    12)
	dataset=/SingleTop_tChannel-madgraph/mangano-SkimStep005-bkg12-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg12-$label
	;;
    13)
	dataset=/SingleTop_sChannel-madgraph/mangano-SkimStep005-bkg13-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg13-$label
	;;
    14)
	dataset=/SingleTop_tWChannel-madgraph/mangano-SkimStep005-bkg14-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg14-$label
	;;
    15)
	dataset=/WW/mangano-SkimStep005-bkg15-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg15-$label
	;;
    16)
	dataset=/WJets-madgraph/mangano-SkimStep005-bkg16-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg16-$label
	;;
    17)
	dataset=/Zee/mangano-SkimStep005-bkg17-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg17-$label
	;;
    18)
	dataset=/Zmumu/mangano-SkimStep005-bkg18-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg18-$label
	;;
    19)
	dataset=/Ztautau/mangano-SkimStep005-bkg19-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg19-$label
	;;
    20)
	dataset=/TTbarJets-madgraph/mangano-SkimStep005-bkg20-oct15-2455c73529abf8a64e80894749518608/USER 
	taskName=EventProd006-bkg20-$label
	;;
    ### signals
    115)
	dataset=/H160_2W_2lnu_gluonfusion_7TeV/mangano-SkimStep005-sig115-oct15-2455c73529abf8a64e80894749518608/USER
	taskName=EventProd006-sig115-$label
	;;


    *)
	echo "puppa"
	;;
esac
#


pset=eventProducer_ALL_cfg
globalTag='START38_V12::All'


case $action in
    submit) echo "action is submit crab jobs"
	cat $DIR/test/$pset.py |sed "s#SET_GLOBALTAG#$globalTag#" > pset.$taskName.py       

	cat $DIR/crab/crab.template.cfg |sed "s#SET_DATASET#$dataset#" | \
	    sed "s#SET_TASK_NAME#$taskName#" > tmp.cfg
	cat tmp.cfg| sed "s#SET_PSET#$PWD/pset.$taskName.py#" > crab.$taskName.cfg; rm tmp.cfg
	crab -create -cfg crab.$taskName.cfg
	rm crab.$taskName.cfg;
	#rm pset.$taskName.py*;
	crab -submit -c $taskName
	;;

    copy) 
	echo "action = copy is not implemented yet"
	;;


    *) echo "ERROR: input type of action is not recognized. Choose between \"submit\" and \"copy\" "
	exit 1
	;;
esac

