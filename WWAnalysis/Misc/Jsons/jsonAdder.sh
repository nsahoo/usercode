#!/bin/bash

# Need FWCore/PythonUtilities
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile

if [ $1 -eq 1 ]; then
 echo "puppa1"   
 # take the intersection of May10 samples 
 compareJSON.py --and lumiSummary.150.json        lumiSummary.151.json     output1.json
 compareJSON.py --and lumiSummary.152.json        lumiSummary.153.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.154.json     ucsdCrab.May10ReReco.json
 rm output*.json

 # take the intersection of all v4 PromptReco samples 
 compareJSON.py --and lumiSummary.100.json        lumiSummary.101.json     output1.json
 compareJSON.py --and lumiSummary.102.json        lumiSummary.103.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.104.json     ucsdCrab.PromptRecoV4.json
 rm output*.json


 # take the intersection of Aug05 samples 
 compareJSON.py --and lumiSummary.160.json        lumiSummary.161.json     output1.json
 compareJSON.py --and lumiSummary.162.json        lumiSummary.163.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.164.json     ucsdCrab.Aug05ReReco.json
 rm output*.json

 # take the intersection of all v6a PromptReco samples 
 compareJSON.py --and lumiSummary.121.json        lumiSummary.122.json     output1.json 
 compareJSON.py --and lumiSummary.123.json        lumiSummary.124.json     output2.json
 compareJSON.py --and output1.json        output2.json     ucsdCrab.PromptRecoV6a.json
 rm output*.json

 # take the intersection of all v6b PromptReco samples 
 compareJSON.py --and lumiSummary.131.json        lumiSummary.132.json     output1.json 
 compareJSON.py --and lumiSummary.133.json        lumiSummary.134.json     output2.json
 compareJSON.py --and output1.json        output2.json     ucsdCrab.PromptRecoV6b.json
 rm output*.json

 lumiCalc2.py -i ucsdCrab.May10ReReco.json overview >& summary/ucsdCrab.May10ReReco.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV4.json overview >& summary/ucsdCrab.PromptRecoV4.summary2 &
 lumiCalc2.py -i ucsdCrab.Aug05ReReco.json overview >& summary/ucsdCrab.Aug05ReReco.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV6a.json overview >& summary/ucsdCrab.PromptRecoV6a.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV6b.json overview >& summary/ucsdCrab.PromptRecoV6b.summary2 &


fi

if [ $1 -eq 2 ]; then
    folder=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing
    cp $folder/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt officialCert.May10ReReco.json
    cp $folder/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt officialCert.Aug05ReReco.json

    folder=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt
    cp $folder/Cert_160404-173692_7TeV_PromptReco_Collisions11_JSON.txt officialCert.PromptReco.json
    cp $folder/json_DCSONLY_BadEle17Ele8.txt_171050-171578.txt json_DCSONLY_BadEle17Ele8.txt_171050-171578.json


    compareJSON.py --sub officialCert.PromptReco.json json_DCSONLY_BadEle17Ele8.txt_171050-171578.json officialCert.PromptReco_noDiElectronProblem.json
    compareJSON.py --sub officialCert.Aug05ReReco.json json_DCSONLY_BadEle17Ele8.txt_171050-171578.json officialCert.Aug05ReReco_noDiElectronProblem.json
    rm json_DCSONLY_BadEle17Ele8.txt_171050-171578.json

    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=165071 --max=168437 --output=officialCert.PromptRecoV4.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=172620 --max=172802 --output=officialCert.PromptRecoV6a.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=172803 --max=173692 --output=officialCert.PromptRecoV6b.json

    lumiCalc2.py -i officialCert.May10ReReco.json overview >& summary/officialCert.May10ReReco.summary2 &
    lumiCalc2.py -i officialCert.Aug05ReReco.json overview >& summary/officialCert.Aug05ReReco.summary2 &
    lumiCalc2.py -i officialCert.Aug05ReReco_noDiElectronProblem.json overview >& summary/officialCert.Aug05ReReco_noDiElectronProblem.summary2 &

    lumiCalc2.py -i officialCert.PromptRecoV4.json overview >& summary/officialCert.PromptRecoV4.summary2 &
    lumiCalc2.py -i officialCert.PromptRecoV6a.json overview >& summary/officialCert.PromptRecoV6a.summary2 &
    lumiCalc2.py -i officialCert.PromptRecoV6b.json overview >& summary/officialCert.PromptRecoV6b.summary2 &

    rm officialCert.PromptReco_noDiElectronProblem.json
    rm officialCert.PromptReco.json
fi

if [ $1 -eq 3 ]; then
    compareJSON.py --and officialCert.May10ReReco.json    ucsdCrab.May10ReReco.json    certifiedUCSD.May10ReReco.json
    compareJSON.py --and officialCert.Aug05ReReco_noDiElectronProblem.json    ucsdCrab.Aug05ReReco.json    certifiedUCSD.Aug05ReReco.json
    compareJSON.py --and officialCert.PromptRecoV4.json   ucsdCrab.PromptRecoV4.json   certifiedUCSD.PromptRecoV4.json
    compareJSON.py --and officialCert.PromptRecoV6a.json   ucsdCrab.PromptRecoV6a.json   certifiedUCSD.PromptRecoV6a.json
    compareJSON.py --and officialCert.PromptRecoV6b.json   ucsdCrab.PromptRecoV6b.json   certifiedUCSD.PromptRecoV6b.json


    lumiCalc2.py -i certifiedUCSD.May10ReReco.json overview >& summary/certifiedUCSD.May10ReReco.summary2 &
    lumiCalc2.py -i certifiedUCSD.Aug05ReReco.json overview >& summary/certifiedUCSD.Aug05ReReco.summary2 &

    lumiCalc2.py -i certifiedUCSD.PromptRecoV4.json overview >& summary/certifiedUCSD.PromptRecoV4.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV6a.json overview >& summary/certifiedUCSD.PromptRecoV6a.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV6b.json overview >& summary/certifiedUCSD.PromptRecoV6b.summary2 &
fi



if [ $1 -eq 4 ]; then
    compareJSON.py --or certifiedUCSD.May10ReReco.json certifiedUCSD.PromptRecoV4.json   output1.json
    compareJSON.py --or certifiedUCSD.Aug05ReReco.json certifiedUCSD.PromptRecoV6a.json   output2.json
    compareJSON.py --or output1.json output2.json output3.json
    compareJSON.py --or output3.json certifiedUCSD.PromptRecoV6b.json certifiedLatinos.42X.json
    rm output*.json
    
    lumiCalc2.py -i certifiedLatinos.42X.json overview >& summary/certifiedLatinos.42X.summary2 &
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
