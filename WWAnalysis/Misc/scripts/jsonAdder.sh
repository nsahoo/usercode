#!/bin/bash

# Need FWCore/PythonUtilities
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile

if [ $1 -eq 1 ]; then
 echo "puppa1"   
 # take the intersection of all v4 PromptReco samples 
 compareJSON.py --and lumiSummary.100.json        lumiSummary.101.json     output1.json
 compareJSON.py --and lumiSummary.102.json        lumiSummary.103.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.104.json     ucsdCrab.PromptRecoV4.json
 rm output*.json

 # take the intersection of all ReReco samples 
 compareJSON.py --and lumiSummary.150.json        lumiSummary.151.json     output1.json
 compareJSON.py --and lumiSummary.152.json        lumiSummary.153.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.154.json     ucsdCrab.May10ReReco.json
 rm output*.json

 # take the intersection of all ReReco samples 
 compareJSON.py --and lumiSummary.110.json        lumiSummary.111.json     output1.json
 compareJSON.py --and lumiSummary.112.json        lumiSummary.113.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.114.json     ucsdCrab.PromptRecoV5a.json
 rm output*.json
fi

if [ $1 -eq 2 ]; then
    May10ReRecoCertification=Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.json
    PromptRecoCertification=Cert_160404-172255_7TeV_PromptReco_Collisions11_JSON_no_dielec.json

    compareJSON.py --and  $May10ReRecoCertification   ucsdCrab.May10ReReco.json    certifiedUCSD.May10ReReco.json
    compareJSON.py --and  $PromptRecoCertification    ucsdCrab.PromptRecoV4.json   certifiedUCSD.PromptRecoV4.json
    compareJSON.py --and  $PromptRecoCertification    ucsdCrab.PromptRecoV5a.json   certifiedUCSD.PromptRecoV5.json

    lumiCalc.py -i certifiedUCSD.May10ReReco.json --nowarning overview >& certifiedUCSD.May10ReReco.summary &
    lumiCalc.py -i certifiedUCSD.PromptRecoV4.json --nowarning overview >& certifiedUCSD.PromptRecoV4.summary &
    lumiCalc.py -i certifiedUCSD.PromptRecoV5.json --nowarning overview >& certifiedUCSD.PromptRecoV5.summary &

fi

if [ $1 -eq 3 ]; then
#filterJSON.py --min 160329 --max 163869 Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.json --output Cert_May10ReReco.json
#filterJSON.py --min 165071 --max 168437 Cert_160404-172255_7TeV_PromptReco_Collisions11_JSON_no_dielec.json --output Cert_PromptRecoV4.json
#filterJSON.py --min 170053 --max 172255 Cert_160404-172255_7TeV_PromptReco_Collisions11_JSON_no_dielec.json --output Cert_PromptRecoV5a.json
    lumiCalc.py -i Cert_May10ReReco.json --nowarning overview >& Cert_May10ReReco.summary & 
    lumiCalc.py -i Cert_PromptRecoV4.json --nowarning overview >& Cert_PromptRecoV4.summary & 
    lumiCalc.py -i Cert_PromptRecoV5a.json --nowarning overview >& Cert_PromptRecoV5a.summary & 
fi



if [ $1 -eq 4 ]; then
    lumiCalc2.py -i Cert_May10ReReco.json overview >& Cert_May10ReReco.summary2 &
    lumiCalc2.py -i Cert_PromptRecoV4.json overview >& Cert_PromptRecoV4.summary2 &
    lumiCalc2.py -i Cert_PromptRecoV5a.json overview >& Cert_PromptRecoV5a.summary2 &

    lumiCalc2.py -i ucsdCrab.May10ReReco.json overview >& ucsdCrab.May10ReReco.summary2 &
    lumiCalc2.py -i ucsdCrab.PromptRecoV4.json overview >& ucsdCrab.PromptRecoV4.summary2 &
    lumiCalc2.py -i ucsdCrab.PromptRecoV5a.json overview >& ucsdCrab.PromptRecoV5a.summary2 &

    lumiCalc2.py -i certifiedUCSD.May10ReReco.json overview >& certifiedUCSD.May10ReReco.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV4.json overview >& certifiedUCSD.PromptRecoV4.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV5a.json overview >& certifiedUCSD.PromptRecoV5a.summary2 &
    

fi



if [ $1 -eq 5 ]; then
    compareJSON.py --or certifiedUCSD.PromptRecoV4.json certifiedUCSD.May10ReReco.json  output1

fi


# take the union of PromptRecoV4 and May10ReReco
#compareJSON.py --or  Jsons/certifiedUCSD.May10ReReco.json  Jsons/certifiedUCSD.PromptRecov4.json  output7
#compareJSON.py --or  output7 Jsons/certifiedUCSD.PromptRecov4b.json Jsons/certifiedLatinos.42X.json

# clean up
# rm *.json

# # create the new PU histogram
# estimatePileup.py -i Jsons/certifiedLatinos.json  --maxPileupBin=24 Scales/certifiedPileUp.root &
# pids="$!"
# estimatePileup.py -i Jsons/certifiedUCSD.json  --maxPileupBin=24 Scales/certifiedPileUpUCSD.root &
# pids="$pids $!"
# wait $pids
# # and update the vector
# python scripts/createPileUpVector.py

# two lines for counting the lumi in the jsons
# for x in Jsons/*.json; do lumiCalc.py -i $x --nowarning recorded > $x.out 2>&1 & done
# for x in Jsons/*.out; do printf "%-30s%20.2f\n" `echo $x | awk -F/ '{print $2}' | awk -F".out" '{print $1}'`  `calc $(grep -A 2 Selected $x | tail -n 1 | c 4)/1000000`; done
