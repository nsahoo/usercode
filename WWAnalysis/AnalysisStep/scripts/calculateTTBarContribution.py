
from ROOT import gROOT, TCanvas, TH1F, TFile, TDirectoryFile, TLegend, gInterpreter
from ROOT import kCyan, kGreen, kMagenta, kBlue, kBlack, kRed, gStyle
from ROOT import TEfficiency, Math

import sys, math, os

gInterpreter.ExecuteMacro( os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/scripts/LatinoStyle.C" )

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

channels = ['wwmumu','wwmuel','wwelmu','wwelel']
types =                 ['TopBMC','TopMuMC', 'TopTagMC']
labels = dict(zip(types,['b-tag', 'soft-#mu','combined']))
colors =                [ kBlue,  kRed      , kBlack]

# hist = "023.TTJetsMad"
hist = "101160.ggToH160toWWto2L2Nu"

c1 = TCanvas()
gStyle.SetOptStat(0)

# First make the SMURF plot
numHists = {}
denHists = {}
divHists = {}
drawString = ""
l = TLegend(0.65,0.55,0.9,0.8)
l.SetBorderSize(     0);
l.SetFillColor (     0);
l.SetTextAlign (    12);
l.SetTextFont  (    42);
l.SetTextSize  (  0.03);
l.SetHeader("t#bar{t} MC")
l.SetFillColor(0)

for type in types:
    for chan in channels:
        thisNum = "eventHists/bycut/"+chan+type+"/nJets/08/"+hist
        thisDen = "eventHists/bycut/"+chan+type+"/nJets/07/"+hist

        h = f.Get(thisNum)
        if type in numHists:
            numHists[type].Add(h)
        else: 
            numHists[type] = h.Clone()

        h = f.Get(thisDen)
        if type in denHists:
            denHists[type].Add(h)
        else: 
            denHists[type] = h.Clone()
    
    divHists[type] = numHists[type].Clone( "div"+type )
    divHists[type].Divide(numHists[type],denHists[type])
    divHists[type].SetLineColor(colors[-1])
    divHists[type].SetLineWidth(2)
    divHists[type].Draw(drawString)
    divHists[type].GetXaxis().SetRangeUser(0,5)
    divHists[type].SetMinimum(0.0)
    divHists[type].SetMaximum(1.2)
    divHists[type].GetYaxis().SetTitle("Efficiency")
    divHists[type].GetXaxis().SetNdivisions(106)
    divHists[type].SetTitle("CMS preliminary")
    l.AddEntry(divHists[type],labels[type],"l")
    drawString = "SAME"
    colors.pop()

l.Draw()
c1.Print("topTagging.smurfFig20.png")
c1.Print("topTagging.smurfFig20.eps")

sigma = 0.68540158589942957
alpha = (1.0 - sigma)/2.

def getEffandError(num,den):
    eff = 0 if num == 0 else num/den
#     hi = TEfficiency.ClopperPearson(int(den),int(num),sigma,True)
#     lo  = TEfficiency.ClopperPearson(int(den),int(num),sigma,False)
#     return (eff, lo, hi)
    return (100*eff, 100*math.sqrt( eff*(1-eff)/den ))

def getEffFromHist(type):
    num = 0
    den = 0
    for chan in channels:
        h = yieldDir.Get(chan+type).Get(hist)
        num += h.GetBinContent(h.GetNbinsX())
        den += h.GetBinContent(h.GetNbinsX()-1)
    return getEffandError(num,den) 
    
print "Actual MC b-tag:    & $%.2f\pm%.2f$ \\\\" % getEffandError(numHists['TopBMC'].GetBinContent(1),denHists['TopBMC'].GetBinContent(1))
print "Actual MC soft-\mu: & $%.2f\pm%.2f$ \\\\" % getEffandError(numHists['TopMuMC'].GetBinContent(1),denHists['TopMuMC'].GetBinContent(1))
print "Actual MC top-tag:  & $%.2f\pm%.2f$ \\\\" % getEffandError(numHists['TopTagMC'].GetBinContent(1),denHists['TopTagMC'].GetBinContent(1)),
print "\hline"

print "Simulated Data b-tag:    & $%.2f\pm%.2f$ \\\\" % getEffFromHist('TopB')
print "Simulated Data soft-\mu: & $%.2f\pm%.2f$ \\\\" % getEffFromHist('TopMu')
print "Simulated Data top-tag:  & $%.2f\pm%.2f$ \\\\" % getEffFromHist('TopTag'),
print "\hline"

# Actual Data b-tag efficiency (n-1/n-2 of TopB from Data):

# Actual Data soft muon efficiency (n-1/n-2 of TopMu from Data):

# Actual Data top-tag efficiency (n-1/n-2 of TopTag from Data):

# Do some calculations

# Also want to eventually make the purity plot

# 
# 
# for type in types:
#     print type,num[type],den[type],num[type]/den[type]


# print 'TopBMC',
# for i in range(1, num['TopBMC'  ].GetNbinsX()+1):
#     print num['TopBMC'].GetBinContent(i),
# print
# 
# print 'TopMuMC',
# for i in range(1, num['TopMuMC'  ].GetNbinsX()+1):
#     print num['TopMuMC'].GetBinContent(i),
# print
# 
# print 'TopTagMC',
# for i in range(1, num['TopTagMC'  ].GetNbinsX()+1):
#     print num['TopTagMC'].GetBinContent(i),
# print

