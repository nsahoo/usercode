from math import *
import re
import ROOT
ROOT.gROOT.SetBatch(True)

class CutsFile:
    def __init__(self,txtfile):
        self._cuts = []
        file = open(txtfile, "r")
        if not file: raise RuntimeError, "Cannot open "+txtfile+"\n"
        for line in file:
            (name,cut) = line.split(":")
            self._cuts.append((name.strip(),cut.strip()))
    def cuts(self):
        return self._cuts[:]
    def nMinusOneCuts(self):
        ret = []
        for c in self._cuts:
            nm1 = " && ".join("(%s)" % x[1] for x in self._cuts if x[0] != c[0])
            ret.append(("all but "+name, nm1))
        return ret
    def allCuts(self):
        return " && ".join("(%s)" % x[1] for x in self._cuts)

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
        for cn,cv in cuts:
            if cut: cut += " && "
            cut += "(%s)" % cv
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
                if self._weight:
                    print "%7.2f  %6.2f%%   " % (nev, fraction * 100),
                else:
                    print "%7d  %6.2f%%   " % (nev, fraction * 100),
            print ""
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
            tree.Draw("%s>>%s(%s)" % (expr,name,bins), cut ,"goff")
            histo = ROOT.gROOT.FindObject(name).Clone(name)
            return histo

class MCAnalysis:
    def __init__(self,samples,options):
        self._foutName = fout = samples.replace(".txt","")+".root"
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
    def getPlots(self,expr,name,bins,cut):
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
        report = ttylist[0].getYields(cuts)
        for tty in ttylist[1:]: 
            self._mergeReport(report, tty.getYields(cuts))
        return report
    def _getPlots(self,expr,name,bins,ttylist,cut):
        plots = ttylist[0].getPlots(expr,name,bins,cut)
        for tty in ttylist[1:]: 
            self._mergePlots(plots, tty.getPlots(expr,name,bins,cut))
    def _mergeReport(self,one,two):
        for i,(c,x) in enumerate(two):
            for j,xj in enumerate(x):
                one[i][1][j][1] += xj[1]
    def _mergePlots(self,one,two):
        for i,(k,h) in enumerate(two):
            one[i][1].Add(h)
    def _fOut(self,dir=None):
        if dir == None:
            if not self._fout: self._fout = ROOT.TFile.Open(self._foutName, "RECREATE")
            return self._fout
        else:
            self._fOut()
            dir = self._fout.GetDirectory(dir)
            if dir: return dir
            else:   return self._fout.mkdir(dir)
        
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
    parser.add_option("-t", "--tree",   dest="tree", default='%sTree', help="Pattern for tree name");
    parser.add_option("-D", "--dump",   dest="dump", action="store_true", help="Dump events passing selection");
    parser.add_option("-w", "--weight", dest="weight", action="store_true", help="Use weight (in MC events)");
    parser.add_option("-l", "--lumi",   dest="lumi", type="float", default="1.0", help="Luminosity (in 1/fb)");
    parser.add_option("-m", "--mass",   dest="mass", type="int", default="160", help="Higgs boson mass");
    parser.add_option("-f", "--final",  dest="final", action="store_true", help="Just compute final yield after all cuts");
    parser.add_option("-I", "--inclusive",  dest="inclusive", action="store_true", help="Only print out yields for all channels combined");
    (options, args) = parser.parse_args()
    tty = TreeToYield(args[0],options) if ".root" in args[0] else MCAnalysis(args[0],options)
    cf  = CutsFile(args[1])
    report = tty.getYields(cf.cuts())
    tty.prettyPrint(report)
    if options.dump: tty.dumpEvents(cf.allCuts())
