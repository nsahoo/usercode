import sys, os, math, re, copy

from WWAnalysis.AnalysisStep.setupSamples_cff import *


# ____             _         __  __       _       
#|  _ \           (_)       |  \/  |     (_)      
#| |_) | ___  __ _ _ _ __   | \  / | __ _ _ _ __  
#|  _ < / _ \/ _` | | '_ \  | |\/| |/ _` | | '_ \ 
#| |_) |  __/ (_| | | | | | | |  | | (_| | | | | |
#|____/ \___|\__, |_|_| |_| |_|  |_|\__,_|_|_| |_|
#             __/ |                               
#            |___/                                


# fileName = sys.argv[0]
# f = TFile(fileName)
f = TFile('data/note120.root')

currentPlot = "mll"
currentCut = "06"
prefix = "eventHists/bycut"
postfix = "CONVLHT/" + currentPlot+"/"+currentCut

# Can loop here, but just using defaults from above for now
# hs = getStackedPlot(f,allMassSamples['120'],channels,prefix,postfix)
# addDataPlot(hs,f,channels,prefix,postfix)
hs = getStackedPlot(f,allMassSamples['120'],['wwmumu'],prefix,postfix)
addDataPlot(hs,f,['wwmumu'],prefix,postfix)
# hs = getStackedPlot(f,allMassSamples['120'],['wwelmu','wwmuel'],prefix,postfix)
# addDataPlot(hs,f,['wwelmu','wwmuel'],prefix,postfix)

c1 = TCanvas("c1","c1")
gStyle.SetOptStat(0)

hs.Draw(5)
c1.Print("lin.png")

# c1.SetLogy()
# hs.Draw(1)
# c1.Print("log.png")

