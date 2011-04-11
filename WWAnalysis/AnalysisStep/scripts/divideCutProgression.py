import sys
import math 

#aaaah fixme
# sys.argv = sys.argv[1:]
# if len(sys.argv) < 1:
#     print "I need a file!"
#     sys.exit()

# fileName = sys.argv[0]
# sys.argv = sys.argv[1:]

# lumi=1000
# if(len(sys.argv)>0):
#     lumi=float(sys.argv[0])
#     sys.argv = sys.argv[1:]
# 
# sys.argv.append( '-b-' )

from ROOT import gROOT, gStyle, TCanvas, TH1F, TFile, TDirectoryFile, TStyle, TLegend
from ROOT import kCyan, kGreen, kMagenta, kBlue, kWhite, kOrange, kViolet, kBlack, kRed

from WWAnalysis.AnalysisStep.data.v03.failureSamples_cff import failureSamples

c1 = TCanvas("c1","c1")
c1.SetBottomMargin(0.2)
gStyle.SetOptStat(0)
gStyle.SetPadBottomMargin(0.2);

#get file #TODO Error Check
# f0 = TFile(fileName)

# channels = ['wwelel','wwelmu','wwmumu','all']
channels = ['wwelel0','wwmumu0']
labels = ['H #rightarrow WW #rightarroww 2e2#nu','H #rightarrow WW #rightarroww 2#mu2#nu']
# channels = ['wwelel']
# channels = ['all']
nbins = 26

denomSamples = [ 
    'ID101160_V80'
]
allSamples = [
    'ID101160_V8A',
    'ID101160_V8B',
    'ID101160_V8C',
    'ID101160_V8D',
    'ID101160_V8E',
    'ID101160_V8F',
    'ID101160_V8G',
    'ID101160_V8H',
    'ID101160_V8I'
]


tempFiles = [ TFile('data/{0}.root'.format(k)) for k in allSamples ] 
allFiles = dict(zip(allSamples,tempFiles))

tempFiles2 = [ TFile('data/{0}.root'.format(k)) for k in denomSamples ] 
denomFiles = dict(zip(denomSamples,tempFiles2))

colors = [kRed,kBlue,kGreen+2,kBlack]

for denom in denomSamples:
    for sample in allSamples:
        plotName = 'breakdown_{0}_{1}'.format(sample[-3:],denom[-3:])
        allPlots = []
        for ich in range(0,len(channels)):

            allFiles[sample].cd()
            thisPlot = allFiles[sample].Get('yields_{0}'.format(channels[ich]))
            tempThisPlot = TH1F(thisPlot)

            denomFiles[denom].cd()
            denomPlot = denomFiles[denom].Get('yields_{0}'.format(channels[ich]))
            tempDenomPlot = TH1F(denomPlot)

            allFiles[sample].cd()
            allPlots.append( tempThisPlot.Clone() )
#             allPlots[-1].Add(tempThisPlot)
            denomFiles[denom].cd()
            allPlots[-1].Divide(tempDenomPlot)
            allPlots[-1].GetYaxis().SetTitle( sample[-3:]+' / V8' )


        legend = TLegend(0.2,0.73,0.5,0.93)
        legend.SetFillColor(0)
        legend.SetLineColor(0)
        legend.SetLineWidth(0)

        for ip in range(0,len(allPlots)):
            allPlots[ip].SetLineColor(0)
            allPlots[ip].SetMarkerColor(colors[ip])
            allPlots[ip].SetFillColor(0)
            allPlots[ip].SetFillStyle(0)
            legend.AddEntry(allPlots[ip],labels[ip],"p")
            if ip == 0:
                allPlots[ip].Draw()
#                 for i in range(nbins):
#                     allPlots[ip].SetBinError(i+1,0.001)
#                     allPlots[ip].GetXaxis().SetBinLabel(i+1,binLabels[i])
                allPlots[ip].LabelsOption('v')
                allPlots[ip].SetMinimum(0.9)
                allPlots[ip].SetMaximum(1.015)
                allPlots[ip].GetYaxis().SetRangeUser(0.8,1.2)
                allPlots[ip].Draw("E")
            else:
                for i in range(nbins):
                    allPlots[ip].SetBinError(i+1,0.001)
                allPlots[ip].Draw("E SAME")

        legend.Draw()
        c1.Print( 'png/'+plotName+'.png' )
        c1.Print( 'eps/'+plotName+'.eps' )
        c1.Print( 'pdf/'+plotName+'.pdf' )
            

