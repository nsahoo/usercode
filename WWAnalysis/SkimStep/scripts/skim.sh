#!/bin/bash


### DOCUMENTATION ###
#
# action:
#  --> submit
#  --> copy
######################

action=$1

DIR="${CMSSW_BASE}/src/WWAnalysis/SkimStep"

label=oct15
dataset=boh

case $2 in
    ### Wjets:
    1)
	dataset=/Wmunu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                       
	taskName=SkimStep005-bkg1-$label
	;;
    2)
	dataset=/Wenu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                        
	taskName=SkimStep005-bkg2-$label
	;;
    3)
	dataset=/Wtaunu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                      
	taskName=SkimStep005-bkg3-$label
	;;

    ### Zjets:
    4)
	dataset=/Zgamma_2l_enriched/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO          
	taskName=SkimStep005-bkg4-$label
	;;

    ### WW:
    5)
	dataset=/WW_2l_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                  
	taskName=SkimStep005-bkg5-$label
	;;

    ### gg->WW:
    6)
	dataset=/gg_WW_2l2nu-GG2WW/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-bkg6-$label
	;;

    ### WZ:
    7)
	dataset=/WZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                          
	taskName=SkimStep005-bkg7-$label
	;;

    ### ZZ:
    8)	
	dataset=/ZZ/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                          
	taskName=SkimStep005-bkg8-$label
	;;

    ### V+gamma:
    9)
	dataset=/Wgamma/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                      
	taskName=SkimStep005-bkg9-$label
	;;
    10)
	dataset=/Zgamma/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                      
	taskName=SkimStep005-bkg10-$label
	;;

    ### ttbar:
    11)
	dataset=/TTbar_2l/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO                    
	taskName=SkimStep005-bkg11-$label
	;;

    ### single-top:
    12)
	dataset=/SingleTop_tChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-bkg12-$label
	;;
    13)
	dataset=/SingleTop_sChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg13-$label
	;;
    14)	
	# num events: 466 437
	dataset=/SingleTop_tWChannel-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-bkg14-$label
	;;

    ### for WW analysis:
    15)
	# num events: 122980
	dataset=/WW/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-bkg15-$label
	;;
    16)
	# num events: 10 068 895
	dataset=/WJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg16-$label
	;;
    17)
	# num events: 2 461 575
	dataset=/Zee/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg17-$label
	;;
    18)
	# num events: 2 111 268
	dataset=/Zmumu/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg18-$label
	;;
    19)
	# num events: 2 195 255
	dataset=/Ztautau/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg19-$label
	;;
    20)
	# num events: 1 483 404
	dataset=/TTbarJets-madgraph/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO 
	taskName=SkimStep005-bkg20-$label
	;;

    ### signals:
    111)
	dataset=/H120_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig111-$label
	;;
    112)
	dataset=/H130_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig112-$label
	;;
    113)
	dataset=/H140_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig113-$label
	;;
    114)
	dataset=/H150_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig115-$label
	;;
    115)
	dataset=/H160_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig115-$label
	;;
    116)
	dataset=/H170_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig116-$label
	;;

    117)
	dataset=/H200_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig117-$label
	;;

    118)
	dataset=/H300_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig118-$label
	;;

    119)
	dataset=/H400_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig119-$label
	;;

    120)
	dataset=/H500_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig120-$label
	;;
    121)
	dataset=/H600_2W_2lnu_gluonfusion_7TeV/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
	taskName=SkimStep005-sig121-$label
	;;
    *)
	echo "puppa"
	;;
esac
#


pset=hwwSkim
globalTag='START3X_V26::All'


case $action in
    submit) echo "action is submit crab jobs"
	#cat $DIR/test/$pset.py |sed "s#SET_GLOBALTAG#$globalTag#" | \
	#    sed "s#SET_OUTPUT#$outputFile#" > pset.$taskName.py
	
	cp $DIR/test/$pset.py  pset.$taskName.py

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

