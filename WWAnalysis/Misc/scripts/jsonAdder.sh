#!/bin/bash

# Need FWCore/PythonUtilities

# take the intersection of all ReReco samples 
compareJSON.py --and Jsons/090.SingleMu2011AReRecoMay10.json        Jsons/091.DoubleMu2011AReRecoMay10.json     output1.json
compareJSON.py --and Jsons/092.DoubleElectron2011AReRecoMay10.json  Jsons/093.MuEG2011AReRecoMay10.json         output2.json
compareJSON.py --and output1.json                                   output2.json            Jsons/ucsdCrabOutput.May10ReReco.json

# take the intersection of all v4 PromptReco samples 
compareJSON.py --and Jsons/084.SingleMuon2011Av4.json Jsons/085.DoubleElectron2011Av4.json   output3.json
compareJSON.py --and Jsons/086.DoubleMuon2011Av4.json Jsons/087.MuEG2011Av4.json             output4.json
compareJSON.py --and output3.json                     output4.json                           Jsons/ucsdCrab42Xv4.json

# get the certified UCSD json for MayReReco
compareJSON.py --and  Jsons/ucsdCrabOutput.May10ReReco.json   Jsons/certified2011.ReReco.json    Jsons/certifiedUCSD.May10ReReco.json

# get the certified UCSD json for PromptRecoV4
compareJSON.py --and  Jsons/ucsdCrab42Xv4.json   Jsons/certified2011.PromptReco.json    Jsons/certifiedUCSD.PromptRecov4.json


# take the union of PromptRecoV4 and May10ReReco
compareJSON.py --or  Jsons/certifiedUCSD.May10ReReco.json  Jsons/certifiedUCSD.PromptRecov4.json  Jsons/certifiedLatinos.42X.json

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
