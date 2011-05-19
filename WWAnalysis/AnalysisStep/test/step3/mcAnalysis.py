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
            field = line.split()
            rootfile = "tree_%s.root" % field[1]
            if len(field) == 2 and field[0] in ["signal", "background", "data"]:
                if field[0] == "data":
                    self._data.append(TreeToYield(rootfile,options))
                elif field[0] == "signal":
                    self._signals.append(TreeToYield(rootfile % options.mass,options))
                elif field[0] == "background":
                    self._backgrounds.append(TreeToYield(rootfile,options))
            else:
                if field[0] in self._allData: self._allData[field[0]].append(TreeToYield(rootfile,options))
                else                        : self._allData[field[0]] =     [TreeToYield(rootfile,options)]
        if len(self._signals) == 0: raise RuntimeError, "No signals!"
        if len(self._backgrounds) == 0: raise RuntimeError, "No backgrounds!"
    def getYields(self,cuts):
        ret = { 'signal':self._getYields(self._signals, cuts),
                'background':self._getYields(self._backgrounds, cuts) }
        if self._data:
            ret['data'] = self._getYields(self._data, cuts)
        for key in self._allData:
            ret[key] = self._getYields(self._allData[key],cuts)
        return ret
    def prettyPrint(self,reports):
        print " ==== SIGNAL ==== "
        self._signals[0].prettyPrint(reports['signal'])
        print " ==== BACKGROUND ==== "
        self._backgrounds[0].prettyPrint(reports['background'])
        if self._data:
            print " ==== DATA ==== "
            self._data[0].prettyPrint(reports['data'])
        for key in self._allData:
            print " ==== {0} ====".format(key)
            self._allData[key][0].prettyPrint(reports[key])
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
        for key in self._allData:
            for (k,h) in self._getPlots(expr,name,bins,cut,self._allData[key]):
                self._fOut(key).WriteTObject(h)
    def dumpEvents(self,cut,vars=['run','lumi','event']):
        for tty in self._data: tty.dumpEvents(cut,vars)
    def _getYields(self,ttylist,cuts):
        return mergeReports([tty.getYields(cuts) for tty in ttylist])
    def _getPlots(self,expr,name,bins,cut,ttylist):
        return mergePlots([[tty.getPlots(expr,name,bins,cut) for tty in ttylist])
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
