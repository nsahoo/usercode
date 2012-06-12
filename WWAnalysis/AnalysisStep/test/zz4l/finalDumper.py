#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit, modules
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

parser = OptionParser(usage="usage: %prog [options] rootfile \nrun with --help to get list of options")
(options, args) = parser.parse_args()

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
        if self.pho1pt > 0 and self.pho2pt > 0: self.type += "+2g"
        elif self.pho1pt + self.pho2pt > 0:     self.type += "+g"
    def __getattr__(self,attr):
        if self.tree.index != self.index: 
            self.tree.GetEntry(self.index)
            self.tree.index = self.index
        return getattr(self.tree, attr)
    def __lt__(self,other): return self.id <  other.id
    def __le__(self,other): return self.id <= other.id
    def __gt__(self,other): return self.id >  other.id
    def __ge__(self,other): return self.id >= other.id
class Dumper:
    def __init__(self,options=None):
        self.options = options
    def __call__(self,ev):
        print "run %06d lumi %4d event %11d : %-8s  mass %6.2f  mz1 %6.2f  mz2 %6.2f  rapidity %+5.3f  m4l %6.2f  massError %4.2f  MELA LD %5.3f " % (ev.run, ev.lumi, ev.event, ev.type, ev.mass, ev.z1mass, ev.z2mass, ev.rap, ev.m4l, ev.massErr, ev.melaLD)
        print "  Z(1,2): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z1mass,ev.z1pt,ev.z1eta,ev.z1rap,ev.z1mll)
        print "  Z(3,4): mass %7.3f  pt %7.3f eta %+5.3f rapidity %+5.3f  mll %7.3f" % (ev.z2mass,ev.z2pt,ev.z2eta,ev.z2rap,ev.z2mll)
        print "  lepton 1: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1ip2d, ev.l1dz)
        print "  lepton 2: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2ip2d, ev.l2dz)
        print "  lepton 3: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3ip2d, ev.l3dz)
        print "  lepton 4: id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dxy %+6.4f dz %+6.4f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4ip2d, ev.l4dz)
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
        
file = ROOT.TFile(args[0])
tree = file.Get("zz4lTree/probe_tree")
events = {}
for i in xrange(tree.GetEntries()):
    ev = Event(tree,i)
    events[ev.id] = ev
evsorted = events.values(); evsorted.sort()
dump = Dumper()
for ev in evsorted: dump(ev)
