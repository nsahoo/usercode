#!/usr/bin/env python

from ROOT import TH1F,TFile
import os

f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/certifiedPileUp.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]

fout = open( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/python/certifiedPileUp_cfi.py",'w' )


print >> fout, "import FWCore.ParameterSet.Config as cms"
print >> fout
print >> fout, "puVectorCertified = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

fout.close()
