#!/usr/bin/env python
from math import *
import re
import os
import ROOT
from copy import *
ROOT.gROOT.SetBatch(True)

from WWAnalysis.AnalysisStep.mcAnalysis import *
from WWAnalysis.AnalysisStep.fomCalculator import *

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] analysis.txt cuts.txt cut expr min max")
addMCAnalysisOptions(parser)
parser.add_option("-0", "--0jetOnly",    dest="zorro", action="store_true", help="Don't try to create also 1 jet selection out of this cuts file.");
parser.add_option("--asimov",     dest="asimov",   action="store_true", default=False, help="Replace observation with expected outcome (always on if lumi != refLumi)")
parser.add_option("--categorize", dest="categories", action="append", default=[], help="categorize")
parser.add_option("--variable",   dest="variable",   type="string", default="2*gammaMRStar:50,0.,200.", help="Expression and binning")
parser.add_option("--cutbased",   dest="cic",   action="store_true", default=False, help="Just do cut based, in categories")
parser.add_option("--name",       dest="name",      type="string", default="shape")
parser.add_option("--mva1j",   dest="mva1j",    help="Attach this MVA for 1-jet events (e.g. BDT_5ch_160)");
parser.add_option("--dataDriven", dest="ddb", action="append", default=[], help="Load the following data-driven backgrounds")
(options, args) = parser.parse_args()

options.out = "hww-{lumi}fb-{name}.mH{mass}.root".format(lumi=options.lumi, name=options.name, mass=options.mass)
options.comboppo = False
channels = ['mumu','muel','elmu','elel']
if options.inclusive: channels = [ 'all' ]
if 'all' in options.ddb: options.ddb = [ 'WJet', 'WW', 'Top', 'DY' ]

mca  = MCAnalysis(args[0],options)
cf0j = CutsFile(args[1],options)
cf1j = CutsFile(cf0j).replace('jet veto', 'one jet', 'njet == 1').replace('soft b-tag veto', 'b-tag veto', 'bveto_ip && nbjet==0')

cats = { '0j':cf0j, '1j':cf1j }
if options.zorro: cats = { '0j':cf0j }
for extracat in options.categories:
    newcats = {}
    for cn, cf in cats.items():
        if extracat=="pt2":
            newcats[cn+"_hpt2"] = CutsFile(cf).add("pt2 > 25","pt2 > 25")
            newcats[cn+"_lpt2"] = CutsFile(cf).add("pt2 < 25","pt2 < 25")
        elif extracat=="dphi":
            newcats[cn+"_hdphi"] = CutsFile(cf).add("dphi > 120", "dphill > 3.14159265/180.*120.")
            newcats[cn+"_ldphi"] = CutsFile(cf).add("dphi < 120", "dphill < 3.14159265/180.*120.")
        elif extracat=="mll":
            newcats[cn+"_hmll"] = CutsFile(cf).add("mll > 100", "mll > 100")
            newcats[cn+"_lmll"] = CutsFile(cf).add("mll < 100", "mll < 100")
        elif extracat=="top":
            newcats[cn]        = cf
            newcats[cn+"_top"] = CutsFile(cf).invert("top veto")
        elif extracat=="Z":
            newcats[cn]        = cf
            newcats[cn+"_Z"] = CutsFile(cf).invert("Z veto")
    cats = newcats

expr, bins = options.variable.split(":")
if options.cic: 
    (expr, bins) = ("0.5", "1,0,1")
    options.keysHist = False

plots = {}
for name, cuts in cats.iteritems():
    if options.mva and options.mva1j: 
        mca.attachMVA(options.mva1j if "1j" in name else options.mva)
    plots[name] = mca.getPlotsForCut("x", expr, bins, cuts.allCuts(), write=False, nodata=True)

if options.asimov:
    for (name,cuts) in cats.iteritems():
        plots[name]['data'] = mergePlots([p for n,p in plots[name].items() if n not in 'ggH vbfH data'], clone=True)
        for (k,h) in  plots[name]['data']: 
            for b in range(1,h.GetNbinsX()+1):
                h.SetBinContent(b, floor(0.5+h.GetBinContent(b)))
else:
    savKeys = options.keysHist; options.keysHist = False
    for (name,cuts) in cats.iteritems():
        dataPlots = mca.getPlotsForCut("x", expr, bins, cuts.allCuts(), write=False, process="data")
        plots[name]['data'] = dataPlots['data']
    options.keysHist = savKeys

## now serialize
shapesFile = options.out if not options.cic else None
if not options.cic:
    for catname, plotmap in plots.iteritems():
        for proc, plotlist in plotmap.iteritems():
            for (k,h) in plotlist: 
                h.SetName(k+catname)
                mca._fOut(proc).WriteTObject(h)

from WWAnalysis.AnalysisStep.datacardWriter import *
builder = DatacardWriter(options)

combcmd  ="combineCards.py "
combcmd0j="combineCards.py "
combcmd1j="combineCards.py "
for catname, plotmap in plots.iteritems():
    jets = 0 if "0j" in catname else 1
    for channel in channels:
        outName = "hww-%sfb-%s.mH%d.%s.txt" % (options.lumi, options.name, options.mass, channel+catname)
        print "Assembling card for mH = %d, channel %s %s --> %s" % (options.mass, channel, catname, outName)
        ## Get yields
        yields = builder.yieldsFromPlots(plotmap,channel,alwaysKeep=(['data']+mca.listSignals()))
        for y in yields:
            print y
        #print yields

        ## Get all nuisances
        nuisanceMap = getCommonSysts(options.mass, channel, jets, options.mass <= 200)
        if len(options.ddb):
            for proc in options.ddb:
                builder.loadDataDrivenYieldsDefault(yields, options.mass, channel+catname, proc)
        else:
            errWW=(0.2 if "WW"   not in mca.listSignals() else 0.0)
            errWJ=(0.5 if "WJet" not in mca.listSignals() else 0.0)
            errDY=(1.0 if "DY"   not in mca.listSignals() else 0.0)
            errTop0j=(0.4 if "Top" not in mca.listSignals() else 0.0)
            errTop1j=(0.2 if "Top" not in mca.listSignals() else 0.0)
            addFakeBackgroundSysts(nuisanceMap, options.mass, channel, jets,
                                   errWW=errWW, errDY=errDY, errTop0j=errTop0j, errTop1j=errTop1j, errWJ=errWJ)
        #print nuisanceMap

        for p,y in yields.iteritems(): y.fillNuisances(nuisanceMap, p, channel, jets)
        #print nuisanceMap

        ## Write datacard
        builder.writeFromYields(yields, nuisanceMap, outName, options.mass, channel+catname, False,
                                shapesFile=shapesFile, signals=mca.listSignals())
        combcmd += "%s=%s " % (channel+catname, outName)
        if jets == 0: combcmd0j += "%s=%s " % (channel+catname, outName)
        else:         combcmd1j += "%s=%s " % (channel+catname, outName)

combcmd   += " > hww-%sfb-%s.mH%d.comb.txt" % (options.lumi, options.name, options.mass)
combcmd0j += " > hww-%sfb-%s.mH%d.comb_0j.txt" % (options.lumi, options.name, options.mass)
combcmd1j += " > hww-%sfb-%s.mH%d.comb_1j.txt" % (options.lumi, options.name, options.mass)
print combcmd;   os.system(combcmd)
print combcmd0j; os.system(combcmd0j)
print combcmd1j; os.system(combcmd1j)
