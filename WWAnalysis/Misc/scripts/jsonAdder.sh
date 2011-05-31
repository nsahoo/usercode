#!/bin/bash

# Need FWCore/PythonUtilities

# take the intersection of all v1 samples 
compareJSON.py --and Jsons/074.SingleMuon2011A.json   Jsons/075.DoubleElectron2011A.json     output1.json
compareJSON.py --and Jsons/076.DoubleMuon2011A.json   Jsons/077.MuEG2011A.json               output2.json
compareJSON.py --and output1.json                     output2.json                           Jsons/ucsdCrab41Xv1.json

# take the intersection of all v2 samples 
compareJSON.py --and Jsons/079.SingleMuon2011Av2.json Jsons/080.DoubleElectron2011Av2.json   output3.json
compareJSON.py --and Jsons/081.DoubleMuon2011Av2.json Jsons/082.MuEG2011Av2.json             output4.json
compareJSON.py --and output3.json                     output4.json                           Jsons/ucsdCrab41Xv2.json

# take the union of v1 and v2
compareJSON.py --or  Jsons/ucsdCrab41Xv1.json         Jsons/ucsdCrab41Xv2.json               Jsons/ucsdCrabOutput.json

# get the certified UCSD json
compareJSON.py --and  Jsons/ucsdCrabOutput.json       Jsons/certified2011.json               Jsons/certifiedUCSD.json

# then take the intersection with ema's json
compareJSON.py --and Jsons/ucsdCrabOutput.json        Jsons/romeCrabOutput.json              combinedCrabOutput.json

# and finally the intersection with ema's json:
compareJSON.py --and combinedCrabOutput.json          Jsons/certified2011.json               Jsons/certifiedLatinos.json

# clean up
rm *.json

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
