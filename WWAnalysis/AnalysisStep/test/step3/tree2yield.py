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

class MCAnalysis:
    def __init__(self,samples,options):
        self._foutName = options.out if options.out else samples.replace(".txt","")+".root"
        self._fout     = None
        self._options = options
        self._data        = []
        self._signals     = []
        self._backgrounds = [] 
        for line in open(samples,'r'):
            field = line.split()
            if len(field) == 2 and field[0] in ["signal", "background", "data"]:
                rootfile = "tree_%s.root" % field[1]
                if field[0] == "data":
                    self._data.append(TreeToYield(rootfile,options))
                elif field[0] == "signal":
                    self._signals.append(TreeToYield(rootfile % options.mass,options))
                elif field[0] == "background":
                    self._backgrounds.append(TreeToYield(rootfile,options))
        if len(self._signals) == 0: raise RuntimeError, "No signals!"
        if len(self._backgrounds) == 0: raise RuntimeError, "No backgrounds!"
    def getYields(self,cuts):
        ret = { 'signal':self._getYields(self._signals, cuts),
                'background':self._getYields(self._backgrounds, cuts) }
        if self._data:
            ret['data'] = self._getYields(self._data, cuts)
        return ret
    def prettyPrint(self,reports):
        print " ==== SIGNAL ==== "
        self._signals[0].prettyPrint(reports['signal'])
        print " ==== BACKGROUND ==== "
        self._backgrounds[0].prettyPrint(reports['background'])
        if self._data:
            print " ==== DATA ==== "
            self._data[0].prettyPrint(reports['data'])
    def getPlots(self,plots,cut):
        for (expr,name,bins) in plots.plots():
            self.getPlotsForCut(expr,name,bins,cut)
    def getPlotsForCut(self,expr,name,bins,cut):
        for (k,h) in self._getPlots(expr,name,bins,cut,self._signals):
            self._fOut("signal").WriteTObject(h)
        for (k,h) in self._getPlots(expr,name,bins,cut,self._backgrounds):
            self._fOut("background").WriteTObject(h)
        if self._data:
            for (k,h) in self._getPlots(expr,name,bins,cut,self._data):
                self._fOut("data").WriteTObject(h)
    def dumpEvents(self,cut,vars=['run','lumi','event']):
        for tty in self._data: tty.dumpEvents(cut,vars)
    def _getYields(self,ttylist,cuts):
        return mergeReports([tty.getYields(cuts) for tty in ttylist])
    def _getPlots(self,expr,name,bins,cut,ttylist):
        plots = ttylist[0].getPlots(expr,name,bins,cut)
        for tty in ttylist[1:]: 
            self._mergePlots(plots, tty.getPlots(expr,name,bins,cut))
        return plots
    def _mergePlots(self,one,two):
        for i,(k,h) in enumerate(two):
            one[i][1].Add(h)
    def _fOut(self,dir=None):
        if dir == None:
            if not self._fout: self._fout = ROOT.TFile.Open(self._foutName, "RECREATE")
            return self._fout
        else:
            self._fOut()
            tdir = self._fout.GetDirectory(dir)
            if tdir: return tdir
            else:    return self._fout.mkdir(dir)
        
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
    parser.add_option("-t", "--tree",   dest="tree", default='%sTree', help="Pattern for tree name");
    parser.add_option("-o", "--out",    dest="out",  help="Output file name. by default equal to input -'.txt' +'.root'");
    parser.add_option("-D", "--dump",   dest="dump", action="store_true", help="Dump events passing selection");
    parser.add_option("-p", "--plots",  dest="plots", type="string", metavar="FILE", help="Make the plots defined in plot file");
    parser.add_option("-w", "--weight", dest="weight", action="store_true", help="Use weight (in MC events)");
    parser.add_option("-l", "--lumi",   dest="lumi", type="float", default="1.0", help="Luminosity (in 1/fb)");
    parser.add_option("-m", "--mass",   dest="mass", type="int", default="160", help="Higgs boson mass");
    parser.add_option("-f", "--final",  dest="final", action="store_true", help="Just compute final yield after all cuts");
    parser.add_option("-i", "--inclusive",  dest="inclusive", action="store_true", help="Only show totals, not each final state separately");
    parser.add_option("-N", "--n-minus-one", dest="nMinusOne", action="store_true", help="Compute n-minus-one yields and plots")
    parser.add_option("-U", "--up-to-cut",   dest="upToCut",   type="string", help="Run selection only up to the cut matched by this regexp, included.") 
    parser.add_option("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
    parser.add_option("-I", "--invert-cut",  dest="cutsToInvert",  action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
    parser.add_option("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
    parser.add_option("-A", "--add-cut",     dest="cutsToAdd",     action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
    (options, args) = parser.parse_args()
    tty = TreeToYield(args[0],options) if ".root" in args[0] else MCAnalysis(args[0],options)
    cf  = CutsFile(args[1],options)
    if options.plots:
        pf = PlotsFile(options.plots, options)
        tty.getPlots(pf, cf.allCuts())
    else:
        report = tty.getYields(cf)
        tty.prettyPrint(report)
    #tty.getPlots("gammaMRStar","gammaMRStar","200,0.,200.",cf.allCuts())
    #if options.dump: tty.dumpEvents(cf.allCuts())
