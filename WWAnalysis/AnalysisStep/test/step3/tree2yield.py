from math import *
import re
import ROOT

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
        cfmt = "%%-%ds   " % clen;
        print cfmt % "   cut",
        for (hypo,nev) in report[0][1]:
            print "      %4s        " % hypo,
        print ""
        print "-"*(18*5+clen+3)
        for i,(cut,yields) in enumerate(report):
            print cfmt % cut,
            for j,(hypo,nev) in enumerate(yields):
                den = report[i-1][1][j][1] if i>0 else 0
                fraction = nev/float(den) if den > 0 else 1
                if self._options.weight:
                    print "%7.2f  %6.2f%%   " % (nev, fraction * 100),
                else:
                    print "%7d  %6.2f%%   " % (nev, fraction * 100),
            print ""
    def getPlots(self,expr,name,bins,cut):
        return [ (k,self._getPlot(t,expr,name,bins,cut)) for (k,t) in self._trees ]
    def dumpEvents(self,cut,vars=['run','lumi','event']):
        for (k,t) in self._trees:
            print "Dump for channel ",k
            t.Scan(":".join(vars), cut)
            print
    def _getYields(self,cut):
        yields = [ (k,self._getYield(t,cut)) for (k,t) in self._trees ]
        yields += [ ('all', sum(x for h,x in yields)) ]
        return yields
    def _getYield(self,tree,cut):
        if self._options.weight:
            histo = self._getPlot(tree,"0.5","dummy","(1,0.,1.)",cut)
            return histo.GetBinContent(1) 
        else: 
            npass = tree.Draw("1",cut,"goff");
            return npass
    def _getPlot(self,tree,expr,name,bins,cut):
            if self._options.weight: cut = "weight*("+cut+")"
            tree.Draw("%s>>%s(%s)" % (expr,name,bins), cut ,"goff")
            histo = ROOT.gROOT.FindObject(name).Clone(name)
            return histo

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
    parser.add_option("-t", "--tree",   dest="tree", default='%sTree', help="Pattern for tree name");
    parser.add_option("-D", "--dump",   dest="dump", action="store_true", help="Dump events passing selection");
    parser.add_option("-w", "--weight", dest="weight", action="store_true", help="Use weight (in MC events)");
    (options, args) = parser.parse_args()
    tty = TreeToYield(args[0],options)
    cf  = CutsFile(args[1])
    report = tty.getYields(cf.cuts())
    tty.prettyPrint(report)
    if options.dump: tty.dumpEvents(cf.allCuts())
