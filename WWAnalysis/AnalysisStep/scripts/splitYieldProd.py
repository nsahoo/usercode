# print sys.argv
# args = sys.argv[1:]
import getopt
import sys

from WWAnalysis.AnalysisStep.data.v03.nov4Samples_cff       import nov4Samples
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples
from WWAnalysis.AnalysisStep.data.v03.wJetsSamples_cff      import wJetsSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff         import wwSamples


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                    "an:vswjh",
                                   ["all",
                                    "files-per-job=",
                                    "nov4Samples",
                                    "sep172010BSamples",
                                    "wwSamples",
                                    "wJetsSamples",
                                    "help"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(sys.argv[0])
        sys.exit(1)

    
    doAll = False
    doNov4 = False
    doSep17 = False
    doWW = True
    doWJet = False
    filesPer = 1

    for o, a in opts:
        if o in ('-a','--all'):
            doAll = True
        elif o in ('-n', '--files-per-job'):
            filesPer = int(a)
        elif o in ('-v', '--nov4Samples'):
            doNov4 = True;
            doSep17 = doWW = doWJet = False;
        elif o in ('-s', '--sep172010BSamples'):
            doSep17 = True;
            doNov4 = doWW = doWJet = False;
        elif o in ('-w', '--wwSamples'):
            doWW = True;
            doNov4 = doSep17 = doWJet = False;
        elif o in ('-j', '--wJetsSamples'):
            doWJet = True;
            doNov4 = doSep17 = doWW = False;
        elif o in ('-h', '--help'):
            usage(sys.argv[0])
            sys.exit(0)  
        else:
            assert False, "unhandled option"

    if doAll:
        doNov4 = doSep17 = doWW = doWJet = True;

    
    
    myDicts = {}
    if doNov4:
        myDicts.update(nov4Samples) 
    if doSep17:
        myDicts.update(sep172010BSamples) 
    if doWJet:
        myDicts.update(wJetsSamples) 
    if doWW:
        myDicts.update(wwSamples) 


    for key in myDicts:
        for i in range(0,len(myDicts[key][3]),filesPer):

            f = open( "job."+key+"."+myDicts[key][0]+"."+str(i)+".py" , 'w')
            f.write(prefix)
            if( myDicts[key][2] ) : f.write(dataLumiString)

            files = myDicts[key][3][i:min(i+filesPer,len(myDicts[key][3]))]
            fileString = ""
            for file in files:
                fileString += "'" + str(file) + "'" + ", "
            thisString =  'process.yieldParams.inputParams = cms.PSet(\n'
            thisString += '    %s = cms.PSet(\n' % myDicts[key][0]
            thisString += '        files = cms.vstring(%s),\n'  % fileString
            thisString += '        scale = cms.double(%e),\n' % myDicts[key][1]
            thisString += '        printEvents = cms.bool(%s),\n' % str(myDicts[key][2])
            thisString += '    )\n'
            thisString += ')\n'

            f.write(thisString)
            f.write("process.yieldParams.fileOutName = cms.string('%s')\n" % ("job."+key+"."+myDicts[key][0]+"."+str(i)+".root"))
            f.close()


prefix = """

import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
process = cms.Process("HistoProducer")
process.load("WWAnalysis.AnalysisStep.yieldProducer_cfi")
process.yieldParams = process.FWLiteParams.clone()

process.yieldParams.selectionParams.wwmumu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu = process.yieldParams.selectionParams.wwelel.clone()
process.yieldParams.selectionParams.wwelmu.mZ   = cms.double(-1.0) 
process.yieldParams.selectionParams.wwelmu.pMet = cms.double(20.0) 

"""

dataLumiString = """
lumis = LumiList(filename = '/nfs/bluearc/group/jsons/Skim_Mu_And_El_byGiovanni/Final.json')
process.yieldParams.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.yieldParams.lumisToProcess = lumis.getCMSSWString().split(',')
process.yieldParams.resFile = 'WWAnalysis/AnalysisStep/data/START38_V13_AK5PF_L2L3Residual.txt'

"""





def usage(a):
    print """
Usage:

%s [OPTIONS] 

Options:

 -h, --help              This message.
 -a, --all               Use all available files (trumps all) (default False)
 -n, --files-per-job     Number of files per job  (default 1)
 -v, --nov4Samples       Use Nov 4th (default False)
 -s, --sep172010BSamples Use Sep 17th (default False)
 -w, --wwSamples         Use WW and Backgrounds (default True)
 -j, --wJetsSamples      Use WJets Samples (default False)
 """ % a


if __name__ == "__main__":
    main()

