from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile
from ROOT import kCyan, kGreen, kMagenta, kBlue

import sys
import math 

def myFormat(v,e):
    if e>0:
        le = int(math.floor(math.log(e,10)))+1
        right = -1*min(0,le-2)
    else:
        right = 1
    tString = '{0:.'+'{0}'.format(right)+'f},{1:.'+'{0}'.format(right)+'f}'
    return tString.format(v,e)


if len(sys.argv) < 2:
    print "I need a file!"
    sys.exit()

lumi=35.49
if(len(sys.argv)>=3):
    lumi=float(sys.argv[2])

from WWAnalysis.AnalysisStep.data.v03.higgsSamples_cff import higgsSamples
from WWAnalysis.AnalysisStep.data.v03.failureSamples_cff import failureSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwSamples
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples

f = TFile(sys.argv[1])
rawDir = f.Get('raw')
keys = [key.GetName() for key in rawDir.GetListOfKeys()]

nbins = 22
channels = ['wwelel','wwelmu','wwmumu','all']

backgrounds = [ wwSamples[k][0] for k in wwSamples ]
signals = [ higgsSamples[k][0] for k in higgsSamples ] 
failures = [ failureSamples[k][0] for k in failureSamples ] 
datasamples = [ sep172010BSamples[k][0] for k in sep172010BSamples ]
allsamples = backgrounds + signals + datasamples + failures

enhancements =      dict(zip(signals, [ higgsSamples[k][4] for k in higgsSamples ]))
enhancements.update(dict(zip(failures,[ failureSamples[k][4] for k in failureSamples ])))

scales =       dict(zip(signals,     [ higgsSamples[k][1]       for k in higgsSamples ]))
scales.update( dict(zip(failures,     [ failureSamples[k][1]     for k in failureSamples     ])))
scales.update( dict(zip(datasamples, [ sep172010BSamples[k][1]  for k in sep172010BSamples  ])))
scales.update( dict(zip(backgrounds, [ wwSamples[k][1]          for k in wwSamples          ])))

bkgYields = {}
sigYields = {}
enYields = {}
dataYields = {}
bkgErr2 = {}
sigErr2 = {}
enErr2 = {}
dataErr2 = {}

countsRaw = {}
countsScaled = {}
errorsSquaredScaled = {}
for channel in channels:
    countsRaw[channel] = {}
    countsScaled[channel] = {}
    errorsSquaredScaled[channel] = {}
    bkgYields[channel] = {}
    sigYields[channel] = {}
    enYields[channel] = {}
    dataYields[channel] = {}
    bkgErr2[channel] = {}
    sigErr2[channel] = {}
    enErr2[channel] = {}
    dataErr2[channel] = {}
    for i in range(1,nbins+1):
        bkgYields[channel][i] = 0
        sigYields[channel][i] = 0
        enYields[channel][i] = 0
        dataYields[channel][i] = 0
        bkgErr2[channel][i] = 0
        sigErr2[channel][i] = 0
        enErr2[channel][i] = 0
        dataErr2[channel][i] = 0

    for k in keys:
        if k.split('_')[0] in channel:
            countsRaw[channel][k] = {}
            countsScaled[channel][k] = {}
            errorsSquaredScaled[channel][k] = {}
            hRaw = rawDir.Get(k)
            sample = k.replace(channel+'_','',1)
            for i in range(1,nbins+1):
                if sample in datasamples:
                    countsRaw[channel][k][i] = hRaw.GetBinContent(i)
                    countsScaled[channel][k][i] = countsRaw[channel][k][i]*scales[sample]
                    errorsSquaredScaled[channel][k][i] = pow( math.sqrt(countsRaw[channel][k][i])*scales[sample], 2)
                    dataYields[channel][i] += countsRaw[channel][k][i]*scales[sample]
                    dataErr2[channel][i] += pow( math.sqrt(countsRaw[channel][k][i])*scales[sample], 2)
                elif sample in signals or sample in failures:
                    countsRaw[channel][k][i] = hRaw.GetBinContent(i)
                    countsScaled[channel][k][i] = countsRaw[channel][k][i]*scales[sample]*lumi
                    errorsSquaredScaled[channel][k][i] = pow( math.sqrt(countsRaw[channel][k][i])*scales[sample]*lumi, 2)
                    sigYields[channel][i] += countsRaw[channel][k][i]*scales[sample]*lumi
                    sigErr2[channel][i] += pow( math.sqrt(countsRaw[channel][k][i])*scales[sample]*lumi, 2)
                    enYields[channel][i] += countsRaw[channel][k][i]*scales[sample]*lumi*enhancements[sample]
                    enErr2[channel][i] += pow( math.sqrt(countsRaw[channel][k][i])*scales[sample]*lumi*enhancements[sample], 2)
                elif sample in backgrounds:
                    countsRaw[channel][k][i] = hRaw.GetBinContent(i)
                    countsScaled[channel][k][i] = countsRaw[channel][k][i]*scales[sample]*lumi
                    errorsSquaredScaled[channel][k][i] = pow( math.sqrt(countsRaw[channel][k][i])*scales[sample]*lumi, 2)
                    bkgYields[channel][i] += countsRaw[channel][k][i]*scales[sample]*lumi
                    bkgErr2[channel][i] += pow( math.sqrt(countsRaw[channel][k][i])*scales[sample]*lumi, 2)

                
linesRaw = {}
linesScaled = {}
for channel in channels:
    linesRaw[channel] = {}
    linesScaled[channel] = {}
    linesRaw[channel][0] = "Cut,"
    linesScaled[channel][0] = "Cut,"
    for i in range(1,nbins+1):
        linesRaw[channel][i] = "{0},".format(i)
        linesScaled[channel][i] = "{0},".format(i)
    
    for k in keys:
        sample = k.replace(channel+'_','',1)
        if channel in k and sample in allsamples:
            linesRaw[channel][0] += "{0},{1},".format(k,k+' Error')
            linesScaled[channel][0] += "{0},{1},".format(k,k+' Error')
            for cut in countsRaw[channel][k]:
                linesRaw[channel][cut] += "{0},".format(myFormat(countsRaw[channel][k][cut],math.sqrt(countsRaw[channel][k][cut])))
            for cut in countsScaled[channel][k]:
                linesScaled[channel][cut] += "{0},".format(myFormat(countsScaled[channel][k][cut],math.sqrt(errorsSquaredScaled[channel][k][cut])))


print "Yield Breakdowns @ {0} pb^-1".format(lumi)
print "Channel,Data,Data Error,SM Signal,SM Error,4th Gen Signal,4th Gen Error,Background,Background Error"
for channel in channels:
    print "{0},{1},{2},{3},{4}".format(
        channel,
        myFormat(dataYields[channel][nbins],math.sqrt(dataErr2[channel][nbins])),
        myFormat(sigYields[channel][nbins], math.sqrt(sigErr2[channel][nbins])),
        myFormat(enYields[channel][nbins],  math.sqrt(enErr2[channel][nbins])),
        myFormat(bkgYields[channel][nbins], math.sqrt(bkgErr2[channel][nbins]))
    )
print

for channel in channels:
    print "Cut Breakdown by number of events for {0}.".format(channel)
    for line in linesRaw[channel]:
        print linesRaw[channel][line]
    print
    
    print "Cut Breakdown by yield to {0} pb^-1 for {1}.".format(lumi,channel)
    for line in linesScaled[channel]:
        print linesScaled[channel][line]
    print
 
