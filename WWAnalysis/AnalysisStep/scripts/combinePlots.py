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
from ROOT import TObject 


from WWAnalysis.AnalysisStep.data.v03.higgsSamples_cff import higgsSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwSamples
from WWAnalysis.AnalysisStep.data.v03.wwSamples_cff import wwGroupings
from WWAnalysis.AnalysisStep.data.v03.sep172010BSamples_cff import sep172010BSamples

from WWAnalysis.AnalysisStep.yieldProducer_cfi import FWLiteParams as yieldParams

c1 = TCanvas("c1","c1")
gStyle.SetOptStat(0)

#get file #TODO Error Check
f = TFile(fileName)

#List of folders in the file, should be a folder for each sample
#temporary fix, remove the raw folder, as it contains the breakdowns
#need to nest it like the others
sampleKeys = [key.GetName() for key in f.GetListOfKeys()]
sampleKeys.remove('raw')

#list of the hists in the cut progression folder
# cutProgressionFolder = f.Get('raw')
# cutProgressionKeys = [key.GetName() for key in cutProgressionFolder.GetListOfKeys()]

channels = ['wwelel','wwelmu','wwmumu','all']
# channels = ['all']
# channels = ['wwelmu']
nbins = 16

#names of all the samples
backgrounds = [ wwSamples[k][0] for k in wwSamples ]
signals = [ higgsSamples[k][0] for k in higgsSamples ]
dataSamples = [ sep172010BSamples[k][0] for k in sep172010BSamples ]
allSamples = backgrounds + signals + dataSamples

#enhancement factors
enhancements =       dict(zip(signals,     [ higgsSamples[k][4]      for k in higgsSamples      ]))
enhancements.update( dict(zip(dataSamples, [ sep172010BSamples[k][4] for k in sep172010BSamples ])))
enhancements.update( dict(zip(backgrounds, [ wwSamples[k][4]         for k in wwSamples         ])))

#scales in 1pb-1
scales =       dict(zip(signals,     [ higgsSamples[k][1]       for k in higgsSamples       ]))
scales.update( dict(zip(dataSamples, [ sep172010BSamples[k][1]  for k in sep172010BSamples  ])))
scales.update( dict(zip(backgrounds, [ wwSamples[k][1]          for k in wwSamples          ])))

#check to make sure every folder is there
for key in sampleKeys:
    if key not in allSamples:
        print "Shit, we have a key which is not a sample: %s" % key
        sys.exit(1)

#build groupings for initial adds a la wwGroupings
allGroupings = {}
allGroupings.update( wwGroupings )
allGroupings["hww"] = [val for val in signals if val in sampleKeys]
allGroupings["data"] = [val for val in dataSamples if val in sampleKeys]

# plots = [ yieldParams.histParams.getParameter(k).variable for k in yieldParams.histParams.parameterNames_() ]
plots = yieldParams.histParams.parameterNames_() 

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
sequence = ['hww','ww','zJets','top','diBoson','wJets', 'data']
forG = dict(zip(sequence, ['histos', 'histo0', 'histo1', 'histo2', 'histo3', 'histo4', 'histo5' ]))
labels = ['HWW','WW','Z+Jets','Top','WW/WZ','W+Jets','35 pb^{-1}']
colors = [kWhite,kOrange-9,kViolet-9,kMagenta-9,kGreen-9,kCyan-9,kBlack]
tempHist = {}
for i in range(len(colors)-1):
    tempHist[i] = TH1F("asdf"+str(i),"asdf",2,0,2)
    tempHist[i].SetLineColor(kBlack)
    tempHist[i].SetLineWidth(2)
    tempHist[i].SetFillColor(colors[i]);
    legend[i/2].AddEntry(tempHist[i],labels[i],"f")

# sample / channel / [byCut,nMinus1] / [eachPlot]_CUT
for channel in channels:
    for subCat in ['byCut','nMinus1']:
        for plot in plots:
            for cut in range(nbins):
                allPlots = {}
                plotName = '{0}_{1}_{2}_{3:0>2}'.format(channel,subCat,plot,cut)
                for group in allGroupings:
                    subPlotName = '{0}_{1}_{2}_{3}_{4:0>2}'.format(group,channel,subCat,yieldParams.histParams.getParameter(plot).variable.value(),cut)
                    allPlots[group] = TH1F(
                        subPlotName, subPlotName,
                        yieldParams.histParams.getParameter(plot).nbins.value(),
                        yieldParams.histParams.getParameter(plot).low.value(),
                        yieldParams.histParams.getParameter(plot).high.value()
                    )
                    allPlots[group].GetXaxis().SetTitle(yieldParams.histParams.getParameter(plot).xtitle.value())
                    for sample in allGroupings[group]:
                        thisPlot = f.Get('{0}/{1}/{2}/{3}_{4:0>2}'.format(sample,channel,subCat,yieldParams.histParams.getParameter(plot).variable.value(),cut))
                        if group != 'data' : thisPlot.Scale( enhancements[sample] * scales[sample] * lumi * 0.991 * 0.970 * 0.997 )
                        allPlots[group].Add(thisPlot)
                    if cut == 9 and subCat == 'byCut':
                        gf = TFile(plot+".root","UPDATE")
                        tempH = allPlots[group].Clone(forG[group])
                        tempH.Write(forG[group],TObject.kOverwrite)
                        gf.Close()

                for i in range(len(sequence)-1):
                    for j in range(i+1,len(sequence)-1):
                        if sequence[i] != 'data' and sequence[j] != 'data':  
                            allPlots[sequence[i]].Add(allPlots[sequence[j]])
                    allPlots[sequence[i]].Rebin(10)
                    allPlots[sequence[i]].SetLineColor(kBlack);
                    allPlots[sequence[i]].SetLineWidth(2);
                    allPlots[sequence[i]].SetFillColor(colors[i]);
                    if i==0 : 
                        allPlots[sequence[i]].GetYaxis().SetTitle("Entries / {0:.2}".format(allPlots[sequence[i]].GetXaxis().GetBinWidth(1)))
                        allPlots[sequence[i]].SetMinimum(0.0)
                        allPlots[sequence[i]].SetMaximum(1.2*allPlots[sequence[i]].GetMaximum())
                        allPlots[sequence[i]].GetYaxis().SetRangeUser(0.0,1.2*allPlots[sequence[i]].GetMaximum())
                        allPlots[sequence[i]].Draw("HIST")
                    else : 
                        allPlots[sequence[i]].Draw("HIST SAME")

                allPlots['data'].Rebin(10)
                allPlots['data'].SetLineColor(kBlack);
                allPlots['data'].SetLineWidth(2);
                allPlots['data'].Draw("E SAME")
                for l in legend:
                    legend[l].Draw()
                c1.SetLogy(False)
                c1.Print( channel+'/png/lin/'+plotName+'.png' )
                c1.Print( channel+'/eps/lin/'+plotName+'.eps' )
                c1.Print( channel+'/pdf/lin/'+plotName+'.pdf' )
                for i in range(len(sequence)-1):
                    allPlots[sequence[i]].SetMinimum(0.005)
                    allPlots[sequence[i]].SetMaximum(2*allPlots[sequence[i]].GetMaximum())
                    allPlots[sequence[i]].GetYaxis().SetRangeUser(0.005,2*allPlots[sequence[i]].GetMaximum())
                    if i==0: allPlots[sequence[i]].Draw("HIST")
                    else   : allPlots[sequence[i]].Draw("HIST SAME")
                allPlots['data'].Draw("E SAME")
                for l in legend:
                    legend[l].Draw()
                c1.SetLogy(True)
                c1.Print( channel+'/png/log/'+plotName+'.png' )
                c1.Print( channel+'/eps/log/'+plotName+'.eps' )
                c1.Print( channel+'/pdf/log/'+plotName+'.pdf' )




            
