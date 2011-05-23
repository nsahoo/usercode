#!/usr/bin/env python
from WWAnalysis.AnalysisStep.tree2yield import *

class MCAnalysis:
    def __init__(self,samples,options):
        self._foutName = options.out if options.out else samples.replace(".txt","")+".root"
        self._fout     = None
        self._options = options
        self._allData     = {}
        self._data        = []
        self._signals     = []
        self._backgrounds = [] 
        for line in open(samples,'r'):
            field = line.split(':')
            rootfile = "tree_%s.root" % field[1].strip()
            signal = ("%d" in rootfile)
            if signal: rootfile = rootfile % options.mass
            tty = TreeToYield(rootfile, options)
            if len(field) == 3: tty.setScaleFactor(float(field[2]))
            if signal: 
                self._signals.append(tty)
            elif field[0] == "data":
                self._data.append(tty)
            else:
                self._backgrounds.append(tty)
            if field[0] in self._allData: self._allData[field[0]].append(TreeToYield(rootfile,options))
            else                        : self._allData[field[0]] =     [TreeToYield(rootfile,options)]
        if len(self._signals) == 0: raise RuntimeError, "No signals!"
        if len(self._backgrounds) == 0: raise RuntimeError, "No backgrounds!"
    def getYields(self,cuts,nodata=False):
        ret = { }
        for key in self._allData:
            if key == 'data' and nodata: continue
            ret[key] = self._getYields(self._allData[key],cuts)
        return ret
    def prettyPrint(self,reports):
        for key in self._allData:
            print " ==== {0} ====".format(key)
            self._allData[key][0].prettyPrint(reports[key])
    def getPlots(self,plots,cut):
        for (name,expr,bins) in plots.plots():
            self.getPlotsForCut(expr,name,bins,cut)
    def getPlotsForCut(self,expr,name,bins,cut,write=True,nodata=False):
        report = {}
        for key in self._allData:
            if key == 'data' and nodata: continue
            report[key] = self._getPlots(expr,name,bins,cut,self._allData[key])
            for (k,h) in report[key]:
                if write: self._fOut(key).WriteTObject(h)
        return report;
    def dumpEvents(self,cut,vars=['run','lumi','event']):
        for tty in self._data: tty.dumpEvents(cut,vars)
    def _getYields(self,ttylist,cuts):
        return mergeReports([tty.getYields(cuts) for tty in ttylist])
    def _getPlots(self,expr,name,bins,cut,ttylist):
        return mergePlots([tty.getPlots(expr,name,bins,cut) for tty in ttylist])
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
    addTreeToYieldOptions(parser)
    parser.add_option("-o", "--out",    dest="out",  help="Output file name. by default equal to input -'.txt' +'.root'");
    parser.add_option("-D", "--dump",   dest="dump", action="store_true", help="Dump events passing selection");
    parser.add_option("-p", "--plots",  dest="plots", type="string", metavar="FILE", help="Make the plots defined in plot file");
    parser.add_option("-m", "--mass",   dest="mass", type="int", default="160", help="Higgs boson mass");
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
