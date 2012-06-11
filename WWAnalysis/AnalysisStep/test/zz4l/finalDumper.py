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
        print "run %06d lumi %4d event %11d : %-8s  mass %6.2f  mz1 %6.2f  mz2 %6.2f  m4l %6.2f  massError %4.2f LD %5.3f " % (ev.run, ev.lumi, ev.event, ev.type, ev.mass, ev.z1mass, ev.z2mass, ev.m4l, ev.massErr, ev.melaLD)
        print "  lepton 1 id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dz %+5.3f" % (ev.l1pdgId,ev.l1pt,ev.l1eta,ev.l1phi, ev.l1pfIsoComb04EACorr/ev.l1pt, ev.l1sip3d, ev.l1dz)
        print "  lepton 2 id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dz %+5.3f" % (ev.l2pdgId,ev.l2pt,ev.l2eta,ev.l2phi, ev.l2pfIsoComb04EACorr/ev.l2pt, ev.l2sip3d, ev.l2dz)
        print "  lepton 3 id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dz %+5.3f" % (ev.l3pdgId,ev.l3pt,ev.l3eta,ev.l3phi, ev.l3pfIsoComb04EACorr/ev.l3pt, ev.l3sip3d, ev.l3dz)
        print "  lepton 4 id %+2d pt %7.3f eta %+5.3f phi %+5.3f  relIso %5.3f sip3d %5.3f dz %+5.3f" % (ev.l4pdgId,ev.l4pt,ev.l4eta,ev.l4phi, ev.l4pfIsoComb04EACorr/ev.l4pt, ev.l4sip3d, ev.l4dz)
        if ev.pho1pt > 0:
            print "  photon of Z1: pt %7.3f, eta %+5.3f" % (ev.pho1pt, ev.pho1eta)
        if ev.pho2pt > 0:
            print "  photon of Z2: pt %7.3f, eta %+5.3f" % (ev.pho2pt, ev.pho2eta)
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
