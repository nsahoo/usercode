#!/usr/bin/env python

from ROOT import TH1F,TFile
import os

#Data
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/certifiedPileUp.root" )
h = f.Get("pileup")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()

fout = open( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/python/certifiedPileUp_cfi.py",'w' )

print >> fout, "import FWCore.ParameterSet.Config as cms"
print >> fout
print >> fout, "puVectorCertified = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

# MC
f = TFile( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Scales/s4MCPileUp.root" )
h = f.Get("htemp")
dataCertified = [ h.GetBinContent(i)/h.GetSumOfWeights() for i in range(1,h.GetNbinsX()+1) ]
f.Close()


print >> fout
print >> fout, "puS4fromMC = [",
for data in dataCertified:
    print >> fout, data, ", ",
print >> fout, "]"
print >> fout

fout.close()
