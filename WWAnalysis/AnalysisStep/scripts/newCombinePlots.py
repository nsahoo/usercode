import sys, os, math, re, copy

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend, THStack
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack
from ROOT import TObject, gInterpreter, Long
from WWAnalysis.AnalysisStep.scaleFactors_cff import *
from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams

gROOT.LoadMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoPlot.C+" )
gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )
from ROOT import LatinoPlot

# from ROOT import LatinoPlot

masses = [120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600]

# loopables and other default that can probably be configured later
channels =                       ['wwelel'  ,'wwmumu'  ,'wwelmu'  ,'wwmuel'  ]
channelLabel = dict(zip(channels,['ee      ','$\mu\mu$','e$\mu$  ','$\mu$e  ']))
suffixes = ['IPLHT']
labels  =                 [       'HWW',      'WW',        'WW/WZ',  'Z+Jets',        'Top',    'W+Jets', ]
myEnum  = dict(zip(labels,[           0,         1,              4,         2,            3,           5, ]))
colors  = dict(zip(labels,[          10, kOrange-9,       kGreen-9, kViolet-9,   kMagenta-9,     kCyan-9, ]))
samples = dict(zip(labels,[ h160Samples, wwSamples, diBosonSamples, dySamples, ttbarSamples, wJetSamples, ]))
plotInfos = FWLiteParams.histParams.clone()
rebinBy = 6
lumi=1.0

allMassSamples = {}
for mass in masses:
    allMassSamples[mass] = copy.deepcopy(samples)
    allMassSamples[mass]['HWW'] = getattr(sys.modules[__name__],"h{0}Samples".format(mass))
    
sys.argv = sys.argv[1:]
if len(sys.argv) < 1:
    print "I need a file!"
    sys.exit()

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


# def getStackedPlot(prefix,suffix):
#     addablePlots = {}
#     for label in labels:
#         for channel in channels:
#             for key in samples[label]:
#                 plotName = prefix+"/"+channel+suffix+"/"+key+"."+samples[label][key][0]
#                 thisPlot = f.Get(plotName)
#     
#                 if thisPlot == None: 
#                     thisPlot = f.Get(plotName.replace('.','_'))
#     
#                 if thisPlot == None: 
#                     print "WARNING: This histo wasn't found: " + plotName
#                     continue
#     
#                 if label in addablePlots:
#                     addablePlots[label].Add(thisPlot,lumi*samples[label][key][1])
#     
#                 else:
#                     addablePlots[label] = thisPlot.Clone(label+"stack")
#                     addablePlots[label].Scale(lumi*samples[label][key][1])
#     
#     hs = LatinoPlot()
#     hs.setLumi(lumi*1000)
#     re1='(.*?)' # Non-greedy match on filler
#     re2='(\\[)' # Any Single Character 1
#     re3='(.*?)' # Word 1
#     re4='(\\])' # Any Single Character 2
#     rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
#     for label in labels:
#         hs.setMCHist(Long(myEnum[label]),addablePlots[label])
#         # I hate myself for this ...
#         m = rg.search(addablePlots[label].GetXaxis().GetTitle())
#         if m:
#             hs.setLabel(m.group(1))
#             hs.setUnits(m.group(3))
#         else:
#             hs.setBreakdown()
#     return hs
    
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
    
# def getYLabel(currentPlot):
#     re1='.*?'   # Non-greedy match on filler
#     re2='(\\[)' # Any Single Character 1
#     re3='((?:[a-z][a-z]+))' # Word 1
#     re4='(\\])' # Any Single Character 2
#     
#     rg = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)
#     m = rg.search(getattr(plotInfos,currentPlot).xtitle.value())
#     suffix = ""
#     if m:
#         suffix =  m.group(2)
#     
#     range = getattr(plotInfos,currentPlot).high.value() - getattr(plotInfos,currentPlot).low.value()
#     nbins = range / getattr(plotInfos,currentPlot).nbins.value() * rebinBy
#     return ("Events / %d " % nbins)+suffix
    

def getYield(samples,channel):
    suffix = suffixes[0]
    prefix = "eventHists/yields"
    sum = 0
    for key in samples:
        plotName = prefix+"/"+channel+suffix+"/"+key+"."+samples[key][0]
        thisPlot = f.Get(plotName)
        if thisPlot == None: 
            #print "WARNING: This histo wasn't found: " + plotName
            continue
        sum += thisPlot.GetBinContent(thisPlot.GetNbinsX()) * samples[key][1] * lumi
    return sum
    


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
prefix = "eventHists/yields"

# Can loop here, but just using defaults from above for now
# hs = getStackedPlot(prefix,suffixes[0])
# hs.Draw(1)
# c1.Print("lin.png")

# c1.SetLogy()
# hs.Draw(1)
# c1.Print("log.png")

# User Should pass the root file
# fileName = sys.argv[0]
# f = TFile(fileName)


for mass in masses: 
    f = TFile("{0}.root".format(mass))
    print "\\begin{table}[b]\\begin{center}"
    print "\\caption{{Yields for $m_{{H}}$ = {0} \\label{{tab:yield{0} }}}}".format(mass)
    print "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline"
    print " final state ",
    for label in labels:
        print " & ",label,
    print "\\\\ \\hline"
    for chan in channels:
        print channelLabel[chan],
        for label in labels:
            print " & {0:6.2f}".format(getYield(allMassSamples[mass][label],chan)),
        print "\\\\ " 
    print "\\hline",
    print "all     ",
    for label in labels:
        sum = 0
        for chan in channels:
            sum += getYield(allMassSamples[mass][label],chan)
        print " & {0:6.2f}".format(sum),
    print "\\\\ \\hline" 
    print "\\end{tabular}\\end{center}\\end{table}"
    print 
    print
