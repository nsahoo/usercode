#!/usr/bin/env python
from math import *
import re
import ROOT
from copy import *
ROOT.gROOT.SetBatch(True)

class CutsFile:
    def __init__(self,txtfileOrCuts,options=None):
        if type(txtfileOrCuts) == list:
            self._cuts = deepcopy(txtfileOrCuts[:])
        elif isinstance(txtfileOrCuts,CutsFile):
            self._cuts = deepcopy(txtfileOrCuts.cuts())
        else:
            self._cuts = []
            file = open(txtfileOrCuts, "r")
            if not file: raise RuntimeError, "Cannot open "+txtfileOrCuts+"\n"
            for cr,cn,cv in options.cutsToAdd:
                if re.match(cr,"entry point"): self._cuts.append((cn,cv))
            for line in file:
                (name,cut) = [x.strip() for x in line.split(":")]
                if name == "entry point" and cut == "1": continue
                self._cuts.append((name,cut))
                for cr,cn,cv in options.cutsToAdd:
                    if re.match(cr,name): self._cuts.append((cn,cv))
                if options.upToCut and re.search(options.upToCut,name):
                    break
            for ci in options.cutsToInvert:  self.invert(ci)
            for ci in options.cutsToExclude: self.remove(ci)
            for cr,cn,cv in options.cutsToReplace: self.replace(cr,cn,cv)
    def remove(self,cut):
        self._cuts = [(cn,cv) for (cn,cv) in self._cuts if not re.search(cut,cn)]
        return self
    def invert(self,cut):
        for i,(cn,cv) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                if cn.startswith("not ") and re.match(r"!\(.*\)", cv):
                    self._cuts[i] = (cn[4:], cv[2:-1])
                else:
                    self._cuts[i] = ("not "+cn, "!("+cv+")")
        return self
    def replace(self,cut,newname,newcut):       
        for i,(cn,cv) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                self._cuts[i] = (newname, newcut)
        return self
    def cuts(self):
        return self._cuts[:]
    def nMinusOne(self):
        return CutsFile(self.nMinusOneCuts())
    def nMinusOneCuts(self):
        ret = []
        for cn,cv in self._cuts:
            nm1 = " && ".join("(%s)" % cv1 for cn1,cv1 in self._cuts if cn1 != cn)
            ret.append(("all but "+cn, nm1))
        return ret
    def allCuts(self):
        return " && ".join("(%s)" % x[1] for x in self._cuts)

class PlotsFile:
    def __init__(self,txtfileOrPlots,options=None):
        if type(txtfileOrPlots) == list:
            self._plots = txtfileOrPlots[:]
        else:
            self._plots = []
            file = open(txtfileOrPlots, "r")
            if not file: raise RuntimeError, "Cannot open "+txtfileOrPlots+"\n"
            for line in file:
                (name,expr,bins) = [x.strip() for x in line.split(":")]
                self._plots.append((name,expr,bins))
    def plots(self):
        return self._plots[:]

class TreeToYield:
    def __init__(self,root,options,report=None):
        self._fname = root
        self._tfile = ROOT.TFile.Open(root)
        if not self._tfile: raise RuntimeError, "Cannot open %s\n" % root
        self._options = options
        self._trees = []
        for h in ("mumu","muel","elmu","elel",):
            t = self._tfile.Get((options.tree % h)+"/probe_tree")
            if not t: raise RuntimeError, "Cannot find tree %s/probe_tree in file %s\n" % (options.tree % h, root)
            self._trees.append((h,t))
        self._weight  = (options.weight and self._trees[0][1].GetBranch("weight") != None)
    def attachMVA(self,name):
        self._fnameMVA = self._fname.replace(".root","."+name+".root")
        self._tfileMVA = ROOT.TFile.Open(self._fnameMVA)
        if not self._tfileMVA: raise RuntimeError, "Cannot open %s\n" % self._fnameMVA
        self._treesMVA = []
        for h,t0 in self._trees:
            t = self._tfile.Get((options.tree % h)+"/"+name)
            if not t: raise RuntimeError, "Cannot find tree %s/%s in file %s\n" % (options.tree % h, name, self._fnameMVA)
            self._treesMVA.append((h,t))
            t0.AddFriend(t)
    def getYields(self,cuts):
        report = []; cut = ""
        cutseq = [ ['entry point','1'] ]
        sequential = False
        if self._options.nMinusOne: 
            cutseq = cuts.nMinusOneCuts()
            cutseq += [ ['all',cuts.allCuts()] ]
            sequential = False
        elif self._options.final:
            cutseq += [ ['all', cuts.allCuts()] ]
        else:
            cutseq += cuts.cuts();
            sequential = True
        for cn,cv in cutseq:
            if sequential:
                if cut: cut += " && "
                cut += "(%s)" % cv
            else:
                cut = cv
            report.append((cn,self._getYields(cut)))
        return report
    def prettyPrint(self,report):
        clen = max([len(cut) for cut,yields in report])
        nch  = len(report[0][1])
        cfmt = "%%-%ds   " % clen;
        print cfmt % "   cut",
        for (hypo,nev) in report[0][1]:
            print "      %4s        " % hypo,
        print ""
        print "-"*(18*(nch+1)+clen+3)
        for i,(cut,yields) in enumerate(report):
            print cfmt % cut,
            for j,(hypo,nev) in enumerate(yields):
                den = report[i-1][1][j][1] if i>0 else 0
                fraction = nev/float(den) if den > 0 else 1
                if self._options.nMinusOne: 
                    fraction = report[-1][1][j][1]/nev if nev > 0 else 1
                if self._weight and nev < 1000:
                    print "%7.2f  %6.2f%%   " % (nev, fraction * 100),
                else:
                    print "%7d  %6.2f%%   " % (nev, fraction * 100),
            print ""
    def getPlots(self,plots,cut):
        ret = [ [name, self.getPlots(expr,name,bins,cut)] for (expr,name,bins) in plots.plots()]
        return ret
    def getPlots(self,expr,name,bins,cut):
        plots = [ [k,self._getPlot(t,expr,name+"_"+k,bins,cut)] for (k,t) in self._trees ]
        hall  = plots[0][1].Clone(name+"_all"); hall.Reset()
        for k,h in plots: hall.Add(h)
        all   = [ ['all', hall] ]
        if self._options.inclusive:
            plots = all
        else:
            plots += all
        return plots
    def dumpEvents(self,cut,vars=['run','lumi','event']):
        for (k,t) in self._trees:
            print "Dump for channel ",k
            t.Scan(":".join(vars), cut)
            print
    def _getYields(self,cut):
        yields = [ [k,self._getYield(t,cut)] for (k,t) in self._trees ]
        all    = [ ['all', sum(x for h,x in yields)] ]
        if self._options.inclusive:
            yields = all
        else:
            yields += all
        return yields
    def _getYield(self,tree,cut):
        if self._weight:
            histo = self._getPlot(tree,"0.5","dummy","1,0.,1.",cut)
            return histo.GetBinContent(1) 
        else: 
            npass = tree.Draw("1",cut,"goff");
            return npass
    def _getPlot(self,tree,expr,name,bins,cut):
            if self._weight: cut = "weight*"+str(self._options.lumi)+"*("+cut+")"
            nev = tree.Draw("%s>>%s(%s)" % (expr,"htemp",bins), cut ,"goff")
            if nev == 0:
                (nb,xmin,xmax) = bins.split(",")
                histo = ROOT.TH1F(name,name,int(nb),float(xmin),float(xmax))
            else:
                histo = ROOT.gROOT.FindObject("htemp").Clone(name)
                ROOT.gROOT.FindObject("htemp").Delete()
            return histo

def mergeReports(reports):
    one = reports[0]
    for two in reports[1:]:
        for i,(c,x) in enumerate(two):
            for j,xj in enumerate(x):
                one[i][1][j][1] += xj[1]
    return one

