import sys
import math 
import re

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend, THStack
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack
from ROOT import TObject 
from WWAnalysis.AnalysisStep.scaleFactors_cff import *
from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams

# loopables and other default that can probably be configured later
channels = ['wwelel','wwelmu','wwmuel','wwmumu']
suffixes = ['IPLHT']
labels  =                 [       'HWW',      'WW',  'Z+Jets',        'Top',        'WW/WZ',    'W+Jets', ]
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


def setUpLegend(colors,labels):
    legend = {}
    legend[0] = TLegend(0.2,0.8,0.4,0.9)
    legend[1] = TLegend(0.45,0.8,0.65,0.9)
    legend[2] = TLegend(0.7,0.8,0.9,0.9)
    legend[0].SetFillColor(0)
    legend[0].SetLineColor(0)
    legend[0].SetLineWidth(0)
    legend[1].SetFillColor(0)
    legend[1].SetLineColor(0)
    legend[1].SetLineWidth(0)
    legend[2].SetFillColor(0)
    legend[2].SetLineColor(0)
    legend[2].SetLineWidth(0)
    
    i = 0
    for label in labels:
        tempHist[i] = TH1F("asdf"+str(i),"asdf",2,0,2)
        tempHist[i].SetLineColor(kBlack)
        tempHist[i].SetLineWidth(2)
        tempHist[i].SetFillColor(colors[label]);
        legend[i/2].AddEntry(tempHist[i],label,"f")
        i = i+1
    return legend


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
    
    hs = THStack("hs","hs")
    for label in labels:
        addablePlots[label].SetLineColor(kBlack);
        addablePlots[label].SetLineWidth(2);
        addablePlots[label].SetFillColor(colors[label]);
        addablePlots[label].Rebin(rebinBy)
        hs.Add(addablePlots[label])
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


#Setup the legend
legend = setUpLegend(colors,labels)

currentPlot = "mll"
currentCut = "06"
plotPath = "eventHists/bycut/"+channels[0]+suffixes[0]+"/"+currentPlot+"/"+currentCut+"/"

# Can loop here, but just using defaults from above for now
hs = getStackedPlot(plotPath)

hs.Draw("HIST")
hs.GetXaxis().SetTitle( getattr(plotInfos,currentPlot).xtitle.value() )
hs.GetYaxis().SetTitle( getYLabel(currentPlot) )

#Add the data 
# data = getDataPlot(plotPath)
# data.Draw("SAME")

#Draw the legend
for l in legend: legend[l].Draw()

#Add the data 

#Add a line for preliminary?

# Make room for the stuff at the top
logMax = 1000.0 * hs.GetMaximum()
linMax = 2.0    * hs.GetMaximum()
  
hs.SetMaximum(linMax)
c1.Print("lin.png")
# for fileType in fileTypes:
#     dir = plotPath[:-1].replace("eventHists","plots/"+fileType)+".lin."+fileType)
#     if not os.path.exists(dir):
#         os.makedirs(dir)
#     c1.Print(
  
c1.SetLogy()
hs.SetMaximum(logMax)
hs.SetMinimum(logMax/2e5)
c1.Print("log.png")












