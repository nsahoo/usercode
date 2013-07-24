from math import *
import os,sys
from optparse import OptionParser
from array import array
import subprocess

import ROOT
ROOT.gROOT.SetBatch(1)

parser = OptionParser()#usage="usage: %prog [options] workspace min max \nrun with --help to get list of options")
(options, args) = parser.parse_args()
options.rS = True


muonPullsFile53X  = ROOT.TFile.Open("mcPull_2012_ptEta_ALLFine.root")
muonPullsFile42X  = ROOT.TFile.Open("mcPull_2011_ptEta_ALLFine.root")
muonPulls = {
    'reco53x':muonPullsFile53X.Get("ALLFine_rms"),
    'reco42x':muonPullsFile42X.Get("ALLFine_rms"),
    'mc53x':muonPullsFile53X.Get("ALLFine_rms"),
    'mc42x':muonPullsFile42X.Get("ALLFine_rms"),
}

def makeFrame(name):
    xbin = array('d',[ 5, 7.5, 10, 13, 16, 20, 25, 30, 35, 40, 50, 60, 80, 100, 200 ])
    ybin = array('d',[0.05*i for i in range(51)])
    ret = ROOT.TH2F(name,name+";p_{T} (GeV);|#eta|",len(xbin)-1,xbin,len(ybin)-1,ybin)
    for bx in xrange(1,len(xbin)+1):
        for by in xrange(1,len(ybin)+1):
            ret.SetBinContent(bx,by,1.0)
    return ret

def applyCorrList(frame,bins,ptMin,ptMax):
    nx, ny = frame.GetNbinsX(), frame.GetNbinsY()
    ax, ay = frame.GetXaxis(), frame.GetYaxis()
    for bx in xrange(1,nx+1):
        for by in xrange(1,ny+1):
            (x,y) = (ax.GetBinCenter(bx), ay.GetBinCenter(by))
            if x < ptMin or x > ptMax: continue
            for (ptLo,ptHi,etaMin,etaMax,res,resE) in bins:
                if etaMin <= y and y <= etaMax:
                    frame.SetBinContent(bx,by, res*frame.GetBinContent(bx,by))
    return frame


def applyCorrHist(frame,hist):
    nx, ny = frame.GetNbinsX(), frame.GetNbinsY()
    ax, ay = frame.GetXaxis(), frame.GetYaxis()
    ahx, ahy = hist.GetXaxis(), hist.GetYaxis()
    nhx, nhy = hist.GetNbinsX(), hist.GetNbinsY()
    for bx in xrange(1,nx+1):
        for by in xrange(1,ny+1):
            (x,y) = (ax.GetBinCenter(bx), ay.GetBinCenter(by))
            (bhx, bhy) = (ahx.FindBin(x), ahy.FindBin(y))
            bhx = max(1,min(nhx,bhx))        
            bhy = max(1,min(nhy,bhy))
            frame.SetBinContent(bx,by, hist.GetBinContent(bhx,bhy)*frame.GetBinContent(bx,by))
    return frame


### https://gist.github.com/1027906
def check_output(*popenargs, **kwargs):
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

def glob(X):
    return check_output('ls -1 '+X,shell=True).strip().split('\n')
files = []
files += glob("plots/2011_v7_all/fit2D_ptEta_all_mcPull_jPsi_ptAvg*5.0_020.0*eta*{0.0_0.8,0.8_1.6,1.6_2.4}*txt")
files += glob("plots/2011_v7_all/fit2D_ptEta_all_mcPull_zMuMu_ptAvg*20*99*eta*{0.0_0.8,0.8_1.6,1.6_2.4}*txt")
#files += glob("plots/2012_v7_all/fit2D_ptEta_all_zEE_ptAvg*30*80*eta*{0.0_0.8,0.8_1.5,1.5_2.0,2.0_2.5}*53x*txt")
files += glob("plots/2011_v7_all/fit2D_ptEta_all_zEE_ptAvg*30*80*eta*{0.0_0.8,0.8_1.5,1.5_2.0,2.0_2.5}*42x*txt")
    
db = {}
for filename in files:
    file = open(filename, "r")
    for line in file:
        cols = line.split()
        (dm,dmE,res,resE) = [float(x) for x in cols[2:6]]
        hcols = cols[0].split("_")
        (ptMin,ptMax,etaMin,etaMax) = [float(hcols[i]) for i in (2,3,5,6)]
        sample = hcols[0]+"_"+hcols[-1]
        if sample not in db: db[sample] = []
        db[sample].append( (ptMin,ptMax,etaMin,etaMax,res,resE) )
print db.keys()

ROOT.gROOT.ProcessLine(".x ~/RooLogon.C");
ROOT.gStyle.SetCanvasDefW(800)
ROOT.gStyle.SetPaperSize(40,40.)
c1 = ROOT.TCanvas("c1")
(w,h) = (800,800)
c1.SetWindowSize(w + (w - c1.GetWw()), h + (h - c1.GetWh()));
c1.SetRightMargin(0.2)
c1.SetLogx(1)

def nicePlot(plot,zmin=0.601,zmax=1.59,ztitle="#sigma(real)/#sigma(pred.)",drawOpt="COLZ",drawFmt=".2f"):
        plot.GetXaxis().SetTitle("p_{T} (GeV)");
        plot.GetXaxis().SetMoreLogLabels(1);
        plot.GetYaxis().SetTitle("|#eta|");
        plot.GetYaxis().SetDecimals(1);
        plot.GetXaxis().SetNoExponent(1);
        plot.GetYaxis().SetNdivisions(505);
        plot.GetXaxis().SetLabelOffset(0.00);
        plot.GetZaxis().SetTitle(ztitle);
        plot.GetZaxis().SetRangeUser(zmin,zmax);
        plot.GetZaxis().SetDecimals(1);
        plot.SetContour(100);
        plot.SetMarkerSize(1.7);
        ROOT.gStyle.SetPaintTextFormat(drawFmt);
        ROOT.gStyle.SetTextFont(42);
        ROOT.gStyle.SetOptStat(0);
        plot.Draw(drawOpt);
        Y = 2012 if '53x' in plot.GetName() else 2011
        for E in [ 'pdf', 'eps', 'png' ]:
            c1.Print("plots/%d_v7_all/000_Overall_%s.%s" % (Y,plot.GetName(),E))
    

finals = {}
#for X in 'reco53x', 'mc53x': #, 'reco42x', 'mc42x':
for X in 'reco42x', 'mc42x': 
    finals['mu_'+X] = applyCorrList(applyCorrList(applyCorrHist(makeFrame("mu_"+X),muonPulls[X]),db['jPsi_'+X],5.,20),db['zMuMu_'+X],20.,999)
    finals['el_'+X] = applyCorrList(makeFrame("el_"+X),db['zEE_'+X],0,9999)
    finals['mu_noPull_'  +X] = applyCorrList(applyCorrList(makeFrame("mu_noPull_"+X),db['jPsi_'+X],5.,20),db['zMuMu_'+X],20.,999)
    finals['mu_justPull_'+X] = applyCorrHist(makeFrame("mu_justPull_"+X),muonPulls[X])
fOut = ROOT.TFile("finalCorrections.root","RECREATE")
for X in finals.itervalues(): 
    nicePlot(X)
    fOut.WriteTObject(X)
fOut.Close()


