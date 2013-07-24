from math import *
import os,sys
from optparse import OptionParser
from array import array
import subprocess

import ROOT
ROOT.gROOT.SetBatch(1)

parser = OptionParser()#usage="usage: %prog [options] workspace min max \nrun with --help to get list of options")
(options, args) = parser.parse_args()
options.rS = True

ROOT.gROOT.LoadMacro('RooLogon.C')
ROOT.gROOT.LoadMacro('fitMassCaliforniaPtEtaBinsEBEv2.cxx+')
from ROOT import fitMassCaliforniaPtEtaBinsEBEv2

resonances = [0,1,2,4]
# 0 = Z->mm
# 1 = J/psi->mm
# 2 = Ups->mm
# 4 = Z->ee

for data in range(2):
    if data==0: print "===> RUNNING ON DATA..."
    else: print "===> RUNNING ON MC..."
    for res in resonances:
        print "   ### running resonance ",res
        ijobs = list()
        if res == 0: ijobs = [41,42,43]
        elif res < 3: ijobs = [40,41]
        else: ijobs = [12]
        for ijob in ijobs:
            print "       %%% running bin ",ijob
            runfile='scripts/fit-data'+str(data)+'-res'+str(res)+'-ibin-'+str(ijob)+'.src'
            f = open(runfile,'w')
            f.write('cd /afs/cern.ch/work/e/emanuele/TnP/CMSSW_5_3_7/src\n')
            f.write('eval `scram ru -sh` \n')
            f.write('root -b RooLogon.C << EOF\n')
            f.write('.L fitMassCaliforniaPtEtaBinsEBEv2.cxx+\n')
            f.write('fitMassCaliforniaPtEtaBinsEBEv2('+str(data)+','+str(res)+',7,'+str(ijob)+')\n')
            f.write('.q \n')
            f.write('EOF\n')
            #fitMassCaliforniaPtEtaBinsEBEv2(data,res,7,ijob)
            bsub = 'bsub -q 1nh -J '+runfile+' -o '+os.getcwd()+'/log/job-ch'+runfile+'.log source '+os.getcwd()+'/'+runfile
            os.system(bsub)
            print "       %%% done bin ",ijob
        print "   ### done resonance ",res
    print "===> DONE FOR DATA/MC..."
print "very done."
