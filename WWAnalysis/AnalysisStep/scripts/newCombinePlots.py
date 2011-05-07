import sys, os, math, re

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend, THStack
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack
from ROOT import TObject, gInterpreter, Long
from WWAnalysis.AnalysisStep.scaleFactors_cff import *
from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams

gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )
from ROOT import LatinoPlot

# from ROOT import LatinoPlot

# loopables and other default that can probably be configured later
channels = ['wwelel','wwelmu','wwmuel','wwmumu']
suffixes = ['IPLHT']
labels  =                 [       'HWW',      'WW',  'Z+Jets',        'Top',        'WW/WZ',    'W+Jets', ]
myEnum  = dict(zip(labels,[           0,         1,         2,            3,              4,           5, ]))
colors  = dict(zip(labels,[          10, kOrange-9, kViolet-9,   kMagenta-9,       kGreen-9,     kCyan-9, ]))
samples = dict(zip(labels,[ h160Samples, wwSamples, dySamples, ttbarSamples, diBosonSamples, wJetSamples, ]))
plotInfos = FWLiteParams.histParams.clone()
rebinBy = 6
lumi=28.36/1000.
    
sys.argv = sys.argv[1:]
if len(sys.argv) < 1:
    print "I need a file!"
    sys.exit()

# User Should pass the root file
fileName = sys.argv[0]
f = TFile(fileName)

c1 = TCanvas("c1","c1")
gStyle.SetOptStat(0)
tempHist = {}


# def setUpLegend(colors,labels):
#     legend = {}
#     legend[0] = TLegend(0.2,0.8,0.4,0.9)
#     legend[1] = TLegend(0.45,0.8,0.65,0.9)
#     legend[2] = TLegend(0.7,0.8,0.9,0.9)
#     legend[0].SetFillColor(0)
#     legend[0].SetLineColor(0)
#     legend[0].SetLineWidth(0)
#     legend[1].SetFillColor(0)
#     legend[1].SetLineColor(0)
#     legend[1].SetLineWidth(0)
#     legend[2].SetFillColor(0)
#     legend[2].SetLineColor(0)
#     legend[2].SetLineWidth(0)
#     
#     i = 0
#     for label in labels:
#         tempHist[i] = TH1F("asdf"+str(i),"asdf",2,0,2)
#         tempHist[i].SetLineColor(kBlack)
#         tempHist[i].SetLineWidth(2)
#         tempHist[i].SetFillColor(colors[label]);
#         legend[i/2].AddEntry(tempHist[i],label,"f")
#         i = i+1
#     return legend


def getStackedPlot(plotPath):
    addablePlots = {}
    for label in labels:
        for key in samples[label]:
            thisPlot = f.Get(plotPath+key+"."+samples[label][key][0])
    
            if thisPlot == None: 
                continue
    
            if label in addablePlots:
                addablePlots[label].Add(thisPlot,lumi*samples[label][key][1])
    
            else:
                addablePlots[label] = thisPlot.Clone(label+"stack")
                addablePlots[label].Scale(lumi*samples[label][key][1])
    
    hs = LatinoPlot()
    hs.setLumi(lumi*1000)
    re1='(.*?)' # Non-greedy match on filler
    re2='(\\[)' # Any Single Character 1
    re3='(.*?)' # Word 1
    re4='(\\])' # Any Single Character 2
    rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
    for label in labels:
        hs.setMCHist(Long(myEnum[label]),addablePlots[label])
        # I hate myself for this ...
        m = rg.search(addablePlots[label].GetXaxis().GetTitle())
        hs.setLabel(m.group(1))
        hs.setUnits(m.group(3))
    return hs
    
def getDataPlot(plotPath):
    combinedData = None
    for key in dataSamples:
        thisPlot = f.Get(plotPath+key+"."+dataSamples[key][0])
    
        if thisPlot == None: 
            continue
    
        if combinedData == None:
            combinedData = thisPlot.Clone("data")
        else:
            combinedData.Add(thisPlot)

    combinedData.SetLineColor(kBlack);
    combinedData.SetLineWidth(2);
    return combinedData
    
def getYLabel(currentPlot):
    re1='.*?'   # Non-greedy match on filler
    re2='(\\[)' # Any Single Character 1
    re3='((?:[a-z][a-z]+))' # Word 1
    re4='(\\])' # Any Single Character 2
    
    rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
    m = rg.search(getattr(plotInfos,currentPlot).xtitle.value())
    suffix = ""
    if m:
        suffix =  m.group(2)
    
    range = getattr(plotInfos,currentPlot).high.value() - getattr(plotInfos,currentPlot).low.value()
    nbins = range / getattr(plotInfos,currentPlot).nbins.value() * rebinBy
    return ("Events / %d " % nbins)+suffix
    


# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                



currentPlot = "dPhill"
currentCut = "12"
plotPath = "eventHists/bycut/"+channels[0]+suffixes[0]+"/"+currentPlot+"/"+currentCut+"/"

# Can loop here, but just using defaults from above for now
hs = getStackedPlot(plotPath)
hs.Draw(10)
c1.Print("lin.png")
