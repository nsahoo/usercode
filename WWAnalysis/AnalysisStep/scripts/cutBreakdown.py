import sys
import math 

#aaaah fixme
sys.argv = sys.argv[1:]
if len(sys.argv) < 1:
    print "I need a file!"
    sys.exit()

fileName = sys.argv[0]
sys.argv = sys.argv[1:]
lumi=35.49
if(len(sys.argv)>0):
    lumi=float(sys.argv[0])
    sys.argv = sys.argv[1:]

sys.argv.append( '-b-' )

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack


from WWAnalysis.AnalysisStep.data.v03.higgsSamples_cff import higgsSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwGroupings
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples

from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams as yieldParams

c1 = TCanvas("c1","c1")
c1.SetBottomMargin(0.2)
gStyle.SetOptStat(0)
gStyle.SetPadBottomMargin(0.2);

#get file #TODO Error Check
f = TFile(fileName)

#list of the hists in the cut progression folder
cutProgressionFolder = f.Get('raw')
cutProgressionKeys = [key.GetName().split('_',1)[1] for key in cutProgressionFolder.GetListOfKeys()]

# channels = ['wwelel','wwelmu','wwmumu','all']
channels = ['all']
nbins = 16

#names of all the samples
backgrounds = [ wwSamples[k][0] for k in wwSamples ]
signals = [ higgsSamples[k][0] for k in higgsSamples ]
dataSamples = [ sep172010BSamples[k][0] for k in sep172010BSamples ]
allSamples = backgrounds + signals + dataSamples

#enhancement factors
enhancements =       dict(zip(signals,     [ higgsSamples[k][4]       for k in higgsSamples       ]))
enhancements.update( dict(zip(dataSamples, [ sep172010BSamples[k][4]  for k in sep172010BSamples  ])))
enhancements.update( dict(zip(backgrounds, [ wwSamples[k][4]          for k in wwSamples          ])))

#scales in 1pb-1
scales =       dict(zip(signals,     [ higgsSamples[k][1]       for k in higgsSamples       ]))
scales.update( dict(zip(dataSamples, [ sep172010BSamples[k][1]  for k in sep172010BSamples  ])))
scales.update( dict(zip(backgrounds, [ wwSamples[k][1]          for k in wwSamples          ])))


#build groupings for initial adds a la wwGroupings
allGroupings = {}
allGroupings.update( wwGroupings )
allGroupings["hww"] = [val for val in signals if val in cutProgressionKeys]
# allGroupings["data"] = [val for val in dataSamples if val in cutProgressionKeys]

legend = {}
legend[0] = TLegend(0.2,0.8,0.4,0.9)
legend[1] = TLegend(0.45,0.8,0.65,0.9)
legend[2] = TLegend(0.7,0.8,0.9,0.9)
legend[3] = TLegend(0.55,0.7,0.80,0.75)
legend[0].SetFillColor(0)
legend[0].SetLineColor(0)
legend[0].SetLineWidth(0)
legend[1].SetFillColor(0)
legend[1].SetLineColor(0)
legend[1].SetLineWidth(0)
legend[2].SetFillColor(0)
legend[2].SetLineColor(0)
legend[2].SetLineWidth(0)
legend[3].SetFillColor(0)
legend[3].SetLineColor(0)
legend[3].SetLineWidth(0)
sequence = ['hww','ww','zJets','top','diBoson','wJets']
labels = ['4th Gen','WW','Z+Jets','Top','WW/WZ','W+Jets']
colors = [kWhite,kOrange-9,kViolet-9,kMagenta-9,kGreen-9,kCyan-9]
tempHist = {}
for i in range(len(colors)):
    tempHist[i] = TH1F("asdf"+str(i),"asdf",2,0,2)
    tempHist[i].SetLineColor(kBlack)
    tempHist[i].SetLineWidth(2)
    tempHist[i].SetFillColor(colors[i]);
    legend[i/2].AddEntry(tempHist[i],labels[i],"f")
tempHist['data'] = TH1F("asdf234","asdf",2,0,2)
tempHist['data'].SetLineColor(kBlack)
tempHist['data'].SetLineWidth(2)
legend[3].AddEntry(tempHist["data"],"{0:.1f} ".format(lumi)+"pb^{-1}","lp")

binLabels = [ 
    "Initial Skim",
    "20/20",
    "#slash{E}_{T}",
    "M_{ll}",
    "Z Veto",
    "p#slash{E}_{T}",
    "Jet Veto",
    "Soft Muon",
    "Extra Lepton",
    "Top Veto",
    "M_{ll} 2",
    "p_{T} Max",
    "p_{T} Min",
#     "pMet 2",
#     "#slash{E}/p_{T}^{ll}",
    "d#phi_{ll}",
    "d#phi_{ll}",
    "d#phi_{ll}"
]


# sample / channel / [byCut,nMinus1] / [eachPlot]_CUT
for channel in channels:
    allPlots = {}
    plotName = 'breakdown_{0}'.format(channel)
    for group in allGroupings:
        subPlotName = '{0}_{1}'.format(group,channel)
        allPlots[group] = TH1F( subPlotName, subPlotName, nbins, 0, nbins)
        for bin in range(0,nbins):
            allPlots[group].GetXaxis().SetBinLabel(bin+1,binLabels[bin])
        allPlots[group].GetXaxis().LabelsOption("v")
        for sample in allGroupings[group]:
            thisPlot = f.Get('raw/{0}_{1}'.format(channel,sample))
            thisPlot.Scale( enhancements[sample] * scales[sample] * lumi )
            allPlots[group].Add(thisPlot)
    
    for i in range(len(sequence)):
        for j in range(i+1,len(sequence)):
            allPlots[sequence[i]].Add(allPlots[sequence[j]])
        allPlots[sequence[i]].SetLineColor(kBlack);
        allPlots[sequence[i]].SetLineWidth(2);
        allPlots[sequence[i]].SetFillColor(colors[i]);
        allPlots[sequence[i]].GetYaxis().SetTitle("Events / cut".format(allPlots[sequence[i]].GetXaxis().GetBinWidth(1)))
        allPlots[sequence[i]].SetMinimum(0.01)
        allPlots[sequence[i]].SetMaximum(10*allPlots[sequence[i]].GetMaximum())
        allPlots[sequence[i]].GetYaxis().SetRangeUser(0.01,10*allPlots[sequence[i]].GetMaximum())

        allPlots[sequence[i]].SetBinContent(14, allPlots[sequence[i]].GetBinContent(16) )
        allPlots[sequence[i]].SetBinError(14, allPlots[sequence[i]].GetBinError(16) )
        allPlots[sequence[i]].GetXaxis().SetRangeUser(0,13)
            
        if i==0: allPlots[sequence[i]].Draw("HIST")
        else   : allPlots[sequence[i]].Draw("HIST SAME")

    allPlots["data"] = TH1F( "data_all", "data_all", nbins, 0, nbins)
    for sample in dataSamples:
        thisPlot = f.Get('raw/{0}_{1}'.format(channel,sample))
        allPlots["data"].Add(thisPlot)

    allPlots["data"].SetBinContent(14, allPlots["data"].GetBinContent(16) )
    allPlots["data"].SetBinError(14, allPlots["data"].GetBinError(16) )
    allPlots["data"].GetXaxis().SetRangeUser(0,14)
    allPlots["data"].Draw("SAME E")
    for l in legend:
        legend[l].Draw()
    c1.SetLogy(True)
    c1.Print( 'png/'+plotName+'.png' )
    c1.Print( 'eps/'+plotName+'.eps' )
    c1.Print( 'pdf/'+plotName+'.pdf' )

