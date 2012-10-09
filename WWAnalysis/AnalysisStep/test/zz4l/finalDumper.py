#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit, modules
from optparse import OptionParser
from math import *

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

parser = OptionParser(usage="usage: %prog [options] rootfile [what] \nrun with --help to get list of options")
parser.add_option("-b", "--blind",   dest="blind",  default=False, action="store_true",  help="blind OS")
parser.add_option("-m", "--mass-range", dest="massrange", default=(0,9999), type="float", nargs=2, help="Mass range")
parser.add_option("-r", "--run-range",  dest="runrange", default=(0,99999999), type="float", nargs=2, help="Run range")
parser.add_option("--dr", "--min-dr-cut",  dest="minDR", default=0, type="float", help="Run range")
parser.add_option("-t", "--tree",  dest="tree", default=None, type="string", help="Tree to use")
parser.add_option("-c", "--cut",  dest="cut", default=None, type="string", help="Cut to apply")
parser.add_option("-T", "--type",  dest="type", default=None, type="string", help="Type of events to select")
parser.add_option("-N", "--nicola",   dest="nicola",  default=False, action="store_true",  help="print Nicola's header row")
parser.add_option("-F", "--fudge",   dest="fudge",  default=False, action="store_true",  help="print Nicola's header row")
parser.add_option("--id", "--event-id", dest="eventid", default="False", action="store_true", help="print run:lumi:event for skimming")

(options, args) = parser.parse_args()
what = args[1] if len(args) > 1 else "signal"
if what not in [ "signal", "CRss", "CRos" ]: raise RuntimeError, "Unknown what"

def deltaR(eta1,phi1,eta2,phi2):
    dphi=phi1-phi2;
    while dphi >  pi: dphi -= pi
    while dphi < -pi: dphi += pi
    return hypot(dphi, eta1-eta2)
def mllCut4(ev):
    if (ev.z1mll <= 4 or ev.z2mll <= 4): return False
    if ev.qll13 == 0 and ev.mll13 <= 4: return False
    if ev.qll14 == 0 and ev.mll14 <= 4: return False
    if ev.qll23 == 0 and ev.mll23 <= 4: return False
    if ev.qll24 == 0 and ev.mll24 <= 4: return False
    return True
def minDeltaRZZ(ev):
    return min([deltaR(ev.l1eta,ev.l1phi,ev.l3eta,ev.l3phi),
                deltaR(ev.l2eta,ev.l2phi,ev.l3eta,ev.l3phi),
                deltaR(ev.l2eta,ev.l2phi,ev.l4eta,ev.l4phi),
                deltaR(ev.l1eta,ev.l1phi,ev.l4eta,ev.l4phi)])
class Event:
    def __init__(self,tree,index):
        self.tree = tree
        self.index = index
        self.tree.GetEntry(index); self.tree.index = index
        self.id = (tree.run, tree.lumi, tree.event)
        self.type = "4l"
        if abs(self.l1pdgId) == 13 and abs(self.l3pdgId) == 13: self.type = "4mu"
        if abs(self.l1pdgId) == 13 and abs(self.l3pdgId) == 11: self.type = "2mu2e"
        if abs(self.l1pdgId) == 11 and abs(self.l3pdgId) == 13: self.type = "2e2mu"
        if abs(self.l1pdgId) == 11 and abs(self.l3pdgId) == 11: self.type = "4e"
        self.leptype = self.type
        if self.pho1pt > 0 and self.pho2pt > 0: self.type += "+2g"
        elif self.pho1pt + self.pho2pt > 0:     self.type += "+g"
    def __getattr__(self,attr):
        global options
        if self.tree.index != self.index: 
            self.tree.GetEntry(self.index)
            self.tree.index = self.index
        if options.fudge:
            try:
                return getattr(self.tree, attr)
            except AttributeError:
                return -9.99
        else:
            return getattr(self.tree, attr)
    def __lt__(self,other): return self.id <  other.id
    def __le__(self,other): return self.id <= other.id
    def __gt__(self,other): return self.id >  other.id
    def __ge__(self,other): return self.id >= other.id

class SignalDumper:
    def __init__(self,options=None):
        self.options = options
    def __call__(self,ev):
        if ev.mass < self.options.massrange[0] or ev.mass > self.options.massrange[1]: return
        if ev.run < self.options.runrange[0] or ev.run > self.options.runrange[1]: return
        if options.type and (options.type not in ev.type): return
        if options.cut and not eval(options.cut, globals(), {'ev':ev}): return
        if options.nicola: 
            ltype = ev.leptype
            if ltype == "2mu2e": ltype = "2e2mu"
            print " Run= %d evt= %d ls= %d m%s= %.3f mZ1= %.3f mZ2= %.3f massError= %.3f LD=%.3f" % (ev.run, ev.event, ev.lumi, ltype, ev.mass, ev.z1mass, ev.z2mass, ev.massErr, ev.melaLD)
        print "run %06d lumi %4d event %11d : %-8s  mass %6.2f  mz1 %6.2f  mz2 %6.2f  rapidity %+5.3f  m4l %6.2f  massError %4.2f  MELA LD %5.3f " % (ev.run, ev.lumi, ev.event, ev.type, ev.mass, ev.z1mass, ev.z2mass, ev.rap, ev.m4l, ev.massErr, ev.melaLD)
        if options.eventid: print "\t\"%d:%d:%d\"," % (ev.run, ev.lumi, ev.event)
        print "  Z(1,2): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z1mass,ev.z1pt,ev.z1eta,ev.z1rap,ev.z1mll)
        print "  Z(3,4): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z2mass,ev.z2pt,ev.z2eta,ev.z2rap,ev.z2mll)
        #print "  lep  pdgId    pt     eta      phi   id  relIso     sip3d    dxy      dz "
        #print "   1    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1idNew, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1ip2d, ev.l1dz)
        #print "   2    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2idNew, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2ip2d, ev.l2dz)
        #print "   3    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3idNew, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3ip2d, ev.l3dz)
        #print "   4    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4idNew, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4ip2d, ev.l4dz)
        print "  lepton 1: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1idNew, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1ip2d, ev.l1dz)
        print "  lepton 2: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2idNew, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2ip2d, ev.l2dz)
        print "  lepton 3: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3idNew, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3ip2d, ev.l3dz)
        print "  lepton 4: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4idNew, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4ip2d, ev.l4dz)
        if ev.pho1pt > 0:
            print "  photon of Z1:    pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f deltaR %6.4f" % (ev.pho1pt, ev.pho1eta, ev.pho1phi, ev.pho1iso, ev.pho1dr)
        if ev.pho2pt > 0:
            print "  photon of Z2:    pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f deltaR %6.4f" % (ev.pho2pt, ev.pho2eta, ev.pho2phi, ev.pho2iso, ev.pho2dr)
        print "  other pairs: Z(1,3):  mll %7.3f  charge %+1d  Z(1,4):  mll %7.3f  charge %+1d " % ( ev.mll13, ev.qll13, ev.mll14, ev.qll14)
        print "               Z(2,3):  mll %7.3f  charge %+1d  Z(2,4):  mll %7.3f  charge %+1d " % ( ev.mll23, ev.qll23, ev.mll24, ev.qll24)
        print "               3/4 SF m(ll) > 12 %4s   4/4 m(ll) > 4 %4s   6/6 m(ll) > 4 %4s " % ( "pass" if ev.threeMassCut12SF else "fail", "pass" if ev.fourMassCut4Any else "fail", "pass" if ev.fourMassCut4AnyS else "fail" )
        print "  ZZ angles:   cos(theta1) %+6.4f  cos(theta2) %+6.4f  cos(theta*) %+6.4f   phi %+6.4f  phi1 %+6.4f  phi2 %+6.4f " % (ev.melaCosTheta1, ev.melaCosTheta2, ev.melaCosThetaStar, ev.melaPhi, ev.melaPhi1, ev.melaPhi2)
        print "  other info:  vertices: %2d   pfMet: %6.2f" % (ev.recoVertices, ev.pfmet)
        print ""
        print ""
        print ""

class ControlDumper:
    def __init__(self,what,options=None):
        self.what    = what
        self.options = options
    def __call__(self,ev):
        if ev.mass < self.options.massrange[0] or ev.mass > self.options.massrange[1]: return
        if ev.run < self.options.runrange[0] or ev.run > self.options.runrange[1]: return
        if abs(ev.l3pdgId) != abs(ev.l4pdgId) or ev.z2mass <= 12: return
        if not mllCut4(ev): return
        if "CRss" in self.what and ev.l3pdgId != ev.l4pdgId: return
        if "CRos" in self.what and ev.l3pdgId == ev.l4pdgId: return
        minDRZZ = minDeltaRZZ(ev)
        if options.minDR > 0 and minDRZZ < options.minDR: return
        if options.blind and self.what == "CRos": 
            if (110 <= ev.mass and ev.mass <= 140) or (ev.mass > 300): return
        if options.type and (options.type not in ev.type): return
        if options.cut and not eval(options.cut, globals(), {'ev':ev}): return
        print "run %06d lumi %4d event %11d : %-8s  mass %6.2f  mz1 %6.2f  mz2 %6.2f  rapidity %+5.3f  m4l %6.2f  massError %4.2f  MELA LD %5.3f " % (ev.run, ev.lumi, ev.event, ev.type, ev.mass, ev.z1mass, ev.z2mass, ev.rap, ev.m4l, ev.massErr, ev.melaLD)
        if options.eventid: print "\t\"%d:%d:%d\"," % (ev.run, ev.lumi, ev.event)

        print "  Z(1,2): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z1mass,ev.z1pt,ev.z1eta,ev.z1rap,ev.z1mll)
        print "  Z(3,4): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z2mass,ev.z2pt,ev.z2eta,ev.z2rap,ev.z2mll)
        #print "  lep  pdgId    pt     eta      phi   id  relIso     sip3d    dxy      dz "
        #print "   1    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1idNew, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1ip2d, ev.l1dz)
        #print "   2    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2idNew, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2ip2d, ev.l2dz)
        #print "   3    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3idNew, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3ip2d, ev.l3dz)
        #print "   4    %+2d  %7.3f  %+5.3f  %+5.3f   %d   %5.3f    %5.3f  %+6.4f  %+6.4f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4idNew, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4ip2d, ev.l4dz)
        print "  lepton 1: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1idNew, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1ip2d, ev.l1dz)
        print "  lepton 2: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2idNew, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2ip2d, ev.l2dz)
        print "  lepton 3: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3idNew, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3ip2d, ev.l3dz)
        print "  lepton 4: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  id %d relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4idNew, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4ip2d, ev.l4dz)
        if ev.pho1pt > 0:
            print "  photon of Z1:    pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f deltaR %6.4f" % (ev.pho1pt, ev.pho1eta, ev.pho1phi, ev.pho1iso, ev.pho1dr)
        if ev.pho2pt > 0:
            print "  photon of Z2:    pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f deltaR %6.4f" % (ev.pho2pt, ev.pho2eta, ev.pho2phi, ev.pho2iso, ev.pho2dr)
        print "  other pairs: Z(1,3):  mll %7.3f  charge %+1d  Z(1,4):  mll %7.3f  charge %+1d " % ( ev.mll13, ev.qll13, ev.mll14, ev.qll14)
        print "               Z(2,3):  mll %7.3f  charge %+1d  Z(2,4):  mll %7.3f  charge %+1d " % ( ev.mll23, ev.qll23, ev.mll24, ev.qll24)
        print "               3/4 SF m(ll) > 12 %4s   4/4 m(ll) > 4 %4s   6/6 m(ll) > 4 %4s " % ( "pass" if ev.threeMassCut12SF else "fail", "pass" if ev.fourMassCut4Any else "fail", "pass" if ev.fourMassCut4AnyS else "fail" )
        print "               min deltaR between leptons of different Z's %.6f" % minDRZZ
        print "  ZZ angles:   cos(theta1) %+6.4f  cos(theta2) %+6.4f  cos(theta*) %+6.4f   phi %+6.4f  phi1 %+6.4f  phi2 %+6.4f " % (ev.melaCosTheta1, ev.melaCosTheta2, ev.melaCosThetaStar, ev.melaPhi, ev.melaPhi1, ev.melaPhi2)
        print "  other info:  vertices: %2d   pfMet: %6.2f" % (ev.recoVertices, ev.pfmet)
        print ""
        print ""
        print ""
        
file = ROOT.TFile(args[0])
treename = "zz4lTree/probe_tree" if what == "signal" else "zxTree/probe_tree"
if options.tree: treename = options.tree
tree = file.Get(treename)
events = {}
for i in xrange(tree.GetEntries()):
    ev = Event(tree,i)
    if ev.id not in events: events[ev.id] = ev
evsorted = events.values(); evsorted.sort()
dump = SignalDumper(options) if what == "signal" else ControlDumper(what,options)
for ev in evsorted: dump(ev)
