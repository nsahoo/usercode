import sys, os, math, re, copy

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend, THStack
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack
from ROOT import TObject, gInterpreter, Long
from WWAnalysis.AnalysisStep.scaleFactors_cff import *
from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams
from WWAnalysis.AnalysisStep.cutPSets_cfi import *
from WWAnalysis.AnalysisStep.cutMassDependPSets_cfi import *


gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )
from ROOT import LatinoPlot


# loopables and other default that can probably be configured later
masses = [120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]
lumi = 1000./1000.

# the suffix of whatever the name of the SkimEvent branch was, i.e. wwelelIPLHT
suffixes = ['IPLHT']

# four channels plus their corresponding LaTeX labels
channels =                       ['wwelel'  ,'wwmumu'  ,'wwelmu'  ,'wwmuel'  ]
channelLabel = dict(zip(channels,['ee      ','$\mu\mu$','e$\mu$  ','$\mu$e  ']))

# Sample setup
labels      =                 [       'HWW',      'WW',        'WZ/ZZ',  'Z+Jets',        'Top',    'W+Jets', ]
myEnum      = dict(zip(labels,[           0,         1,              4,         2,            3,           5, ]))
tempSamples = dict(zip(labels,[ h160Samples, wwSamples, diBosonSamples, dySamples, ttbarSamples, wJetSamples, ]))

allMassSamples = {}
for mass in masses:
    allMassSamples[mass] = copy.deepcopy(tempSamples)
    allMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))

def getCombinedHistogram(f,sample,channels,prefix,suffix):
    returnPlot = None
    for channel in channels:
        for key in sample:
            plotName = prefix+"/"+channel+suffix+"/"+key+"."+sample[key][0]
            thisPlot = f.Get(plotName)

            if thisPlot == None: 
                thisPlot = f.Get(plotName.replace('.','_'))

            if thisPlot == None: 
                print "WARNING: {0} wasn't found in {1}".format(plotName,f.GetPath())
                continue

            if returnPlot != None:
                returnPlot.Add(thisPlot,lumi*sample[key][1])
            else:
                returnPlot = thisPlot.Clone(key+"stack")
                returnPlot.Scale(lumi*sample[key][1])
    return returnPlot
    
def getStackedPlot(f,samples,channels,prefix,suffix):
    addablePlots = {}
    for label in labels:
        addablePlots[label] = getCombinedHistogram(f,samples[label],channels,prefix,suffix)
    
    hs = LatinoPlot()
    hs.setLumi(lumi*1000)
    re1='(.*?)' # should be the label
    re2='(\\[)' # left bracket
    re3='(.*?)' # units
    re4='(\\])' # right bracket
    rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
    for label in labels:
        hs.setMCHist(Long(myEnum[label]),addablePlots[label])
        # I hate myself for this ...
        m = rg.search(addablePlots[label].GetXaxis().GetTitle())
        if m:
            hs.setLabel(m.group(1))
            hs.setUnits(m.group(3))
        else:
            hs.setBreakdown()
    return hs
    
# def getDataPlot(plotPath):
#     combinedData = None
#     for key in dataSamples:
#         thisPlot = f.Get(plotPath+key+"."+dataSamples[key][0])
#     
#         if thisPlot == None: 
#             continue
#     
#         if combinedData == None:
#             combinedData = thisPlot.Clone("data")
#         else:
#             combinedData.Add(thisPlot)
# 
#     combinedData.SetLineColor(kBlack);
#     combinedData.SetLineWidth(2);
#     return combinedData
    

def getYieldAndError(f,sample,channels,cut=-1):
    suffix = suffixes[0]
    prefix = "eventHists/yields"
    combinedHist = getCombinedHistogram(f,sample,channels,prefix,suffix)
    bin = cut if cut != -1 else combinedHist.GetNbinsX()
    return (combinedHist.GetBinContent(bin),combinedHist.GetBinError(bin))

def getYield(f,sample,channels,cut=-1):
    return getYieldAndError(f,sample,channels,cut)[0]

def getError(f,sample,channels,cut=-1):
    return getYieldAndError(f,sample,channels,cut)[1]

def convertToRealLatex(s):
    if s.find('#') == -1: return s
    else: return "$"+s.replace('#','\\')+"$"


# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                



# currentPlot = "dPhill"
# currentCut = "12"
# plotPath = "eventHists/bycut/"+channels[0]+suffixes[0]+"/"+currentPlot+"/"+currentCut+"/"

# Can loop here, but just using defaults from above for now
# hs = getStackedPlot(prefix,suffixes[0])

# c1 = TCanvas("c1","c1")
# gStyle.SetOptStat(0)

# hs.Draw(1)
# c1.Print("lin.png")

# c1.SetLogy()
# hs.Draw(1)
# c1.Print("log.png")

# User Should pass the root file
# fileName = sys.argv[0]
# f = TFile(fileName)


for mass in masses: 
    f = TFile("data/{0}.root".format(mass))
    outFile = open("tex/yield{0}.tex".format(mass),'w')
    print >> outFile,  "\\begin{table}[p]\\begin{center}"
    print >> outFile,  "\\caption{{Yields for $m_{{H}}$ = {0} \\label{{tab:yield{0} }}}}".format(mass)
    print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline"
    print >> outFile,  " final state ",
    for label in labels:
        print >> outFile,  " & ",label,
    print >> outFile,  "\\\\ \\hline"
    for chan in channels:
        print >> outFile,  channelLabel[chan],
        for label in labels:
            yieldAndError = getYieldAndError(f,allMassSamples[mass][label],[chan])
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format(yieldAndError[0],yieldAndError[1]),
        print >> outFile,  "\\\\ " 
    print >> outFile,  "\\hline",
    print >> outFile,  "all     ",
    for label in labels:
        yieldAndError = getYieldAndError(f,allMassSamples[mass][label],channels)
        print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] ),
    print >> outFile,  "\\\\ \\hline" 
    print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
    print >> outFile
    print >> outFile
    outFile.close()



for mass in masses:

    fullSet = cloneVPSet(defaultWW)
    addMassDependentCuts(fullSet,getattr(sys.modules[__name__],"hThursday{0}".format(mass)))

    f = TFile("data/{0}.root".format(mass))
    outFile = open("tex/breakdown{0}.tex".format(mass),'w')

    print >> outFile,  "\\begin{table}[p]\\begin{center}"
    print >> outFile,  "\\caption{{Breakdown for gg to Higgs to $2\\ell 2\\nu$ for $m_{{H}}$ = {0} \\label{{tab:breakdown{0} }}}}".format(mass)
    print >> outFile,  "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline"
    print >> outFile,  " cut level ",
    for chan in channels:
        print >> outFile,  " & ",channelLabel[chan],
    print >> outFile,  "& all \\\\ \\hline"
    for cut in [4,8,11,15]:
        print >> outFile,  convertToRealLatex(fullSet[-1].label.value() if cut == -1 else fullSet[cut-1].label.value()),
        for chan in channels:
            yieldAndError = getYieldAndError(f,{'101{0}'.format(mass):allMassSamples[mass]["HWW"]["101{0}".format(mass)]},[chan],cut) 
            print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$".format( yieldAndError[0],yieldAndError[1] )
        yieldAndError = getYieldAndError(f,{'101{0}'.format(mass):allMassSamples[mass]["HWW"]["101{0}".format(mass)]},channels,cut) 
        print >> outFile,  " & ${0:6.2f}\pm{1:6.2f}$ \\\\".format( yieldAndError[0],yieldAndError[1] )
    print >> outFile,  "\\hline" 
    print >> outFile,  "\\end{tabular}\\end{center}\\end{table}"
    print >> outFile
    print >> outFile
    outFile.close()

