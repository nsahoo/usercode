#!/bin/bash

# Need FWCore/PythonUtilities

# take the intersection of all v1 samples 
compareJSON.py --and Jsons/074.SingleMuon2011A.json   Jsons/075.DoubleElectron2011A.json     output1.json
compareJSON.py --and Jsons/076.DoubleMuon2011A.json   Jsons/077.MuEG2011A.json               output2.json
compareJSON.py --and output1.json                     output2.json                           output5.json

# take the intersection of all v2 samples 
compareJSON.py --and Jsons/079.SingleMuon2011Av2.json Jsons/080.DoubleElectron2011Av2.json   output3.json
compareJSON.py --and Jsons/081.DoubleMuon2011Av2.json Jsons/082.MuEG2011Av2.json             output4.json
compareJSON.py --and output3.json                     output4.json                           output6.json

# take the union of v1 and v2
compareJSON.py --or  output5.json                     output6.json                           Jsons/ucsdCrabOutput.json

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
# 
# # and update the vector
# python scripts/createPileUpVector.py
