#!/usr/bin/python
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--recompile",
                  action="store_true", dest="recompile", default=False,
                  help="recompile the fake rate estimator")
parser.add_option("-d", "--drawonly",
                  action="store_true", dest="drawonly", default=False,
                  help="do not run. Only draw the results based on existing ROOT files")
parser.add_option("-e", "--electrons",
                  action="store_true", dest="electrons", default=True,
                  help="run the fake rate for electrons")
parser.add_option("-m", "--muons",
                  action="store_false", dest="electrons", default=False,
                  help="run the fake rate for muons")

(options, args) = parser.parse_args()

if(options.drawonly == False):
    if(options.recompile):
        print 'recompiling...'
        os.system('make clean')
        os.system('make fakerate')
        print 'done complining.'

os.system('./fakerate fr')
os.system('hadd -f fakerates_zll1e.root fr_zlle-ElectronMisid*')
os.system('hadd -f fakerates_zll1m.root fr_zllm-MuonMisid*')
os.system('rm -f fr_zlle-ElectronMisid*')
os.system('rm -f fr_zllm-MuonMisid*')        

from ROOT import gROOT
gROOT.LoadMacro('src/drawFR.cc+')
if(options.electrons == True):
    from ROOT import compareElectronIdsUnbiased
    compareElectronIdsUnbiased()
    print "And now fit with thr paramteric function..."
    from ROOT import fitElectronIdsUnbiased
    fitElectronIdsUnbiased(0)
    fitElectronIdsUnbiased(1)
    fitElectronIdsUnbiased(2)
else:
    from ROOT import compareMuonIdsUnbiased
    compareMuonIdsUnbiased()
    print "And now fit with thr paramteric function..."
    from ROOT import fitMuonIdsUnbiased
    fitMuonIdsUnbiased(0)
    fitMuonIdsUnbiased(1)

# make the 2D map and save in a ROOT file
#gROOT.LoadMacro('src/doFRMaps.cc+')
#if(options.electrons == True):
#    from ROOT import doEleFRMapsHzz4l
#    doEleFRMapsHzz4l(0)
#    doEleFRMapsHzz4l(1)
#else:
#    from ROOT import doMuFRMapsHzz4l
#    doMuFRMapsHzz4l(0)
#    doMuFRMapsHzz4l(1)
