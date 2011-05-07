
from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile, TLegend
from ROOT import kCyan, kGreen, kMagenta, kBlue, kBlack, kRed, gStyle

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


f = TFile(sys.argv[1])
yieldDir = f.Get('eventHists/yields')

# channels = ['wwelel']
# types = ['0']
channels = ['wwmumu','wwmuel','wwelmu','wwelel']
types = ['TopB'  ,'TopMu'  ,'TopTag',
         'TopBMC','TopMuMC','TopTagMC']
# types = ['TopB','TopMu','TopTag']

hist = "023.TTJetsMad"
nbins = 10

num = {}
den = {}
for chan in channels:
    for type in types:
        localDir = yieldDir.Get(chan+type)
        h = localDir.Get(hist)
        if type in num:
            num[type] += h.GetBinContent(h.GetNbinsX())
        else:
            num[type] = h.GetBinContent(h.GetNbinsX())

        if type in den:
            den[type] += h.GetBinContent(h.GetNbinsX()-1)
        else:
            den[type] = h.GetBinContent(h.GetNbinsX()-1)


for type in types:
    print type,num[type],den[type],num[type]/den[type]


num = {}
den = {}
for chan in channels:
    for type in types:
        if 'MC' not in type: continue
        thisNum = "eventHists/bycut/"+chan+type+"/nJets/08/"+hist
        thisDen = "eventHists/bycut/"+chan+type+"/nJets/07/"+hist

        h = f.Get(thisNum)
        if type in num:
            num[type].Add(h)
        else: 
            num[type] = h.Clone()

        h = f.Get(thisDen)
        if type in den:
            den[type].Add(h)
        else: 
            den[type] = h.Clone()


c1 = TCanvas()
gStyle.SetOptStat(0)
# for type in types:
#     if 'MC' not in type: continue
#     num[type].Divide(den[type])
#     num[type].Draw()

    
num['TopBMC'  ].Divide(den['TopBMC'  ])
num['TopMuMC' ].Divide(den['TopMuMC' ])
num['TopTagMC'].Divide(den['TopTagMC'])
num['TopBMC'  ].SetLineColor(kBlack)
num['TopMuMC' ].SetLineColor(kRed)
num['TopTagMC'].SetLineColor(kBlue)
num['TopBMC'  ].SetLineWidth(2)
num['TopMuMC' ].SetLineWidth(2)
num['TopTagMC'].SetLineWidth(2)
num['TopBMC'  ].Draw()
num['TopBMC'  ].GetXaxis().SetRangeUser(0,5)
num['TopBMC'  ].SetMinimum(0.0)
num['TopBMC'  ].SetMaximum(1.5)
num['TopBMC'  ].GetYaxis().SetTitle("Efficiency")
num['TopMuMC' ].Draw("SAME")
num['TopTagMC'].Draw("SAME")


l = TLegend(0.4,0.7,0.9,0.9)
l.SetHeader("TTJets MC by the Thundercats")
l.AddEntry(num['TopBMC'  ],"low pt jet tagging","l")
l.AddEntry(num['TopMuMC' ],"soft muon tagging","l")
l.AddEntry(num['TopTagMC'],"low pt jet and soft muon tagging","l")

print 'TopBMC',
for i in range(1, num['TopBMC'  ].GetNbinsX()+1):
    print num['TopBMC'].GetBinContent(i),
print

print 'TopMuMC',
for i in range(1, num['TopMuMC'  ].GetNbinsX()+1):
    print num['TopMuMC'].GetBinContent(i),
print

print 'TopTagMC',
for i in range(1, num['TopTagMC'  ].GetNbinsX()+1):
    print num['TopTagMC'].GetBinContent(i),
print

l.SetFillColor(0)
l.Draw()


c1.Print("all.png")






