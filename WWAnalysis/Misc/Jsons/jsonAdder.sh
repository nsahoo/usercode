#!/bin/bash

# Need FWCore/PythonUtilities
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile

if [ $1 -eq 1 ]; then
 echo "make the intersaction of all the json files from the crab jobs"   
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
 compareJSON.py --and lumiSummary.120.json        lumiSummary.121.json     output1.json 
 compareJSON.py --and lumiSummary.122.json        lumiSummary.123.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.124.json     ucsdCrab.PromptRecoV6a.json
 rm output*.json

 # take the intersection of all v6b PromptReco samples 
 compareJSON.py --and lumiSummary.130.json        lumiSummary.131.json     output1.json 
 compareJSON.py --and lumiSummary.132.json        lumiSummary.133.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.134.json     ucsdCrab.PromptRecoV6b.json
 rm output*.json

 # take the intersection of all 2001Bv1a PromptReco samples 
 compareJSON.py --and lumiSummary.140a.json        lumiSummary.141a.json     output1.json 
 compareJSON.py --and lumiSummary.142a.json        lumiSummary.143a.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.144a.json     ucsdCrab.PromptReco2011BV1a.json
 rm output*.json

 # take the intersection of all 2001Bv1b PromptReco samples 
 compareJSON.py --and lumiSummary.140b.json        lumiSummary.141b.json     output1.json 
 compareJSON.py --and lumiSummary.142b.json        lumiSummary.143b.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.144b.json     ucsdCrab.PromptReco2011BV1b.json
 rm output*.json

 # take the intersection of all 2001Bv1c PromptReco samples
 compareJSON.py --and lumiSummary.140c.json        lumiSummary.141c.json     output1.json
 compareJSON.py --and lumiSummary.142c.json        lumiSummary.143c.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.144c.json     ucsdCrab.PromptReco2011BV1c.json
 rm output*.json

 # take the intersection of all 2001Bv1d PromptReco samples
 compareJSON.py --and lumiSummary.140d.json        lumiSummary.141d.json     output1.json
 compareJSON.py --and lumiSummary.142d.json        lumiSummary.143d.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.144d.json     ucsdCrab.PromptReco2011BV1d.json
 rm output*.json

 # take the intersection of all 2001Bv1e PromptReco samples
 compareJSON.py --and lumiSummary.140e.json        lumiSummary.141e.json     output1.json
 compareJSON.py --and lumiSummary.142e.json        lumiSummary.143e.json     output2.json
 compareJSON.py --and output1.json        output2.json     output3.json
 compareJSON.py --and output3.json        lumiSummary.144e.json     ucsdCrab.PromptReco2011BV1e.json
 rm output*.json

 lumiCalc2.py -i ucsdCrab.May10ReReco.json overview >& summary/ucsdCrab.May10ReReco.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV4.json overview >& summary/ucsdCrab.PromptRecoV4.summary2 &
 lumiCalc2.py -i ucsdCrab.Aug05ReReco.json overview >& summary/ucsdCrab.Aug05ReReco.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV6a.json overview >& summary/ucsdCrab.PromptRecoV6a.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptRecoV6b.json overview >& summary/ucsdCrab.PromptRecoV6b.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptReco2011BV1a.json overview >& summary/ucsdCrab.PromptReco2011BV1a.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptReco2011BV1b.json overview >& summary/ucsdCrab.PromptReco2011BV1b.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptReco2011BV1c.json overview >& summary/ucsdCrab.PromptReco2011BV1c.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptReco2011BV1d.json overview >& summary/ucsdCrab.PromptReco2011BV1d.summary2 &
 lumiCalc2.py -i ucsdCrab.PromptReco2011BV1e.json overview >& summary/ucsdCrab.PromptReco2011BV1e.summary2 &


fi

if [ $1 -eq 2 ]; then
    folder=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing
    cp $folder/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt officialCert.May10ReReco.json
    cp $folder/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt officialCert.Aug05ReReco.json

    folder=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt
    cp $folder/Cert_160404-178078_7TeV_PromptReco_Collisions11_JSON.txt officialCert.PromptReco.json
    cp $folder/json_DCSONLY_BadEle17Ele8.txt_171050-171578.txt json_DCSONLY_BadEle17Ele8.txt_171050-171578.json


    compareJSON.py --sub officialCert.PromptReco.json json_DCSONLY_BadEle17Ele8.txt_171050-171578.json officialCert.PromptReco_noDiElectronProblem.json
    compareJSON.py --sub officialCert.Aug05ReReco.json json_DCSONLY_BadEle17Ele8.txt_171050-171578.json officialCert.Aug05ReReco_noDiElectronProblem.json
    rm json_DCSONLY_BadEle17Ele8.txt_171050-171578.json

    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=165071 --max=168437 --output=officialCert.PromptRecoV4.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=172620 --max=172802 --output=officialCert.PromptRecoV6a.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=172803 --max=173692 --output=officialCert.PromptRecoV6b.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=175832 --max=176023 --output=officialCert.PromptReco2011BV1a.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=176024 --max=176309 --output=officialCert.PromptReco2011BV1b.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=176310 --max=177053 --output=officialCert.PromptReco2011BV1c.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=177054 --max=177515 --output=officialCert.PromptReco2011BV1d.json
    filterJSON.py officialCert.PromptReco_noDiElectronProblem.json --min=177718 --max=178078 --output=officialCert.PromptReco2011BV1e.json

    lumiCalc2.py -i officialCert.May10ReReco.json overview >& summary/officialCert.May10ReReco.summary2 &
    lumiCalc2.py -i officialCert.Aug05ReReco.json overview >& summary/officialCert.Aug05ReReco.summary2 &
    lumiCalc2.py -i officialCert.Aug05ReReco_noDiElectronProblem.json overview >& summary/officialCert.Aug05ReReco_noDiElectronProblem.summary2 &

    #lumiCalc2.py -i officialCert.PromptRecoV4.json overview >& summary/officialCert.PromptRecoV4.summary2 &
    #lumiCalc2.py -i officialCert.PromptRecoV6a.json overview >& summary/officialCert.PromptRecoV6a.summary2 &
    #lumiCalc2.py -i officialCert.PromptRecoV6b.json overview >& summary/officialCert.PromptRecoV6b.summary2 &
    #lumiCalc2.py -i officialCert.PromptReco2011BV1a.json overview >& summary/officialCert.PromptReco2011BV1a.summary2 &
    #lumiCalc2.py -i officialCert.PromptReco2011BV1b.json overview >& summary/officialCert.PromptReco2011BV1b.summary2 &
    #lumiCalc2.py -i officialCert.PromptReco2011BV1c.json overview >& summary/officialCert.PromptReco2011BV1c.summary2 &
    #lumiCalc2.py -i officialCert.PromptReco2011BV1d.json overview >& summary/officialCert.PromptReco2011BV1d.summary2 &
    lumiCalc2.py -i officialCert.PromptReco2011BV1e.json overview >& summary/officialCert.PromptReco2011BV1e.summary2 &

    rm officialCert.PromptReco_noDiElectronProblem.json
    rm officialCert.PromptReco.json
fi

if [ $1 -eq 3 ]; then
    echo "make the intersaction between crab output jsons and official certification jsons"
    compareJSON.py --and officialCert.May10ReReco.json    ucsdCrab.May10ReReco.json    certifiedUCSD.May10ReReco.json
    compareJSON.py --and officialCert.Aug05ReReco_noDiElectronProblem.json    ucsdCrab.Aug05ReReco.json    certifiedUCSD.Aug05ReReco.json
    compareJSON.py --and officialCert.PromptRecoV4.json   ucsdCrab.PromptRecoV4.json   certifiedUCSD.PromptRecoV4.json
    compareJSON.py --and officialCert.PromptRecoV6a.json   ucsdCrab.PromptRecoV6a.json   certifiedUCSD.PromptRecoV6a.json
    compareJSON.py --and officialCert.PromptRecoV6b.json   ucsdCrab.PromptRecoV6b.json   certifiedUCSD.PromptRecoV6b.json
    compareJSON.py --and officialCert.PromptReco2011BV1a.json   ucsdCrab.PromptReco2011BV1a.json   certifiedUCSD.PromptReco2011BV1a.json
    compareJSON.py --and officialCert.PromptReco2011BV1b.json   ucsdCrab.PromptReco2011BV1b.json   certifiedUCSD.PromptReco2011BV1b.json
    compareJSON.py --and officialCert.PromptReco2011BV1c.json   ucsdCrab.PromptReco2011BV1c.json   certifiedUCSD.PromptReco2011BV1c.json
    compareJSON.py --and officialCert.PromptReco2011BV1d.json   ucsdCrab.PromptReco2011BV1d.json   certifiedUCSD.PromptReco2011BV1d.json
    compareJSON.py --and officialCert.PromptReco2011BV1e.json   ucsdCrab.PromptReco2011BV1e.json   certifiedUCSD.PromptReco2011BV1e.json

    lumiCalc2.py -i certifiedUCSD.May10ReReco.json overview >& summary/certifiedUCSD.May10ReReco.summary2 &
    lumiCalc2.py -i certifiedUCSD.Aug05ReReco.json overview >& summary/certifiedUCSD.Aug05ReReco.summary2 &

    lumiCalc2.py -i certifiedUCSD.PromptRecoV4.json overview >& summary/certifiedUCSD.PromptRecoV4.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV6a.json overview >& summary/certifiedUCSD.PromptRecoV6a.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptRecoV6b.json overview >& summary/certifiedUCSD.PromptRecoV6b.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptReco2011BV1a.json overview >& summary/certifiedUCSD.PromptReco2011BV1a.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptReco2011BV1b.json overview >& summary/certifiedUCSD.PromptReco2011BV1b.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptReco2011BV1c.json overview >& summary/certifiedUCSD.PromptReco2011BV1c.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptReco2011BV1d.json overview >& summary/certifiedUCSD.PromptReco2011BV1d.summary2 &
    lumiCalc2.py -i certifiedUCSD.PromptReco2011BV1e.json overview >& summary/certifiedUCSD.PromptReco2011BV1e.summary2 &
fi



if [ $1 -eq 4 ]; then
    echo "make the or of all the relevant json and produce the final certifiedLatinos.42X.json file"
    compareJSON.py --or certifiedUCSD.May10ReReco.json certifiedUCSD.PromptRecoV4.json   output1.json
    compareJSON.py --or certifiedUCSD.Aug05ReReco.json certifiedUCSD.PromptRecoV6a.json   output2.json
    compareJSON.py --or output1.json output2.json output3.json
    compareJSON.py --or certifiedUCSD.PromptRecoV6b.json certifiedUCSD.PromptReco2011BV1a.json  output4.json
    compareJSON.py --or output3.json output4.json output5.json
    compareJSON.py --or output5.json certifiedUCSD.PromptReco2011BV1b.json output6.json
    compareJSON.py --or output6.json certifiedUCSD.PromptReco2011BV1c.json output7.json
    compareJSON.py --or output7.json certifiedUCSD.PromptReco2011BV1d.json output8.json
    compareJSON.py --or output8.json certifiedUCSD.PromptReco2011BV1e.json certifiedLatinos.42X.json
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
