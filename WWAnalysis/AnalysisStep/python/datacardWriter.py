from WWAnalysis.AnalysisStep.systematicUncertainties import *

class Yield:
    def __init__(self,*args,**kwargs):
        """Things that work:
                Yield(number)  -> create a simple fixed yield
                Yield(number, error,   type="lnN" (default) or "gmN", name="...")
                Yield(N,alpha,         type="gamma",nameN="...")
                Yield(N,alpha,alphaErr,type="gamma", typeAlpha="lnN" (default) or "gmN", nameN="..", nameAlpha="...")"""
        if len(args) == 0: raise RuntimeError, "You must specify at least one value"
        self._args = args
        self._kwargs = kwargs
        if len(args) == 1:
            self.val   = float(args[0])
            self._type = "fixed"
        elif len(args) >= 2:
            self._type = kwargs['type'] if kwargs.has_key('type') else 'lnN'
            if self._type == 'lnN':
                self.val   = float(args[0])
                self._eff = 1+float(args[1])/self.val if self.val else 1.0
                self._name  = kwargs['name']
            if self._type == 'gmN':
                self.val = float(args[0])
                self._eff = float(args[1])/self.val if self.val else 0.0
                self._name  = kwargs['name']
            if self._type == 'gamma':
                self._N     = float(args[0])
                self._alpha = float(args[1])
                self.val   = self._N * self._alpha
                self._nameN  = kwargs['nameN']
                if len(args) >= 3:
                    self._hasErrAlpha = True
                    self._typeAlpha = kwargs['typeAlpha'] if kwargs.has_key('typeAlpha') else 'lnN'
                    self._nameAlpha  = kwargs['nameAlpha']
                    self._effAlpha = 1+float(args[2])/float(args[1]) if self._typeAlpha == 'lnN' else float(args[2])/float(args[1])
                else:
                    self._hasErrAlpha = False
        self.zero = ((self.val == 0) if self._type != "gamma" else (self._alpha == 0))
    def fillNuisances(self,nuisanceMap, process, channel, jets):
        if self._type == 'fixed': return
        elif self._type != 'gamma':
            name = self._name.format(process=process, channel=channel, jets=jets)
            if not nuisanceMap.has_key(name): nuisanceMap[name] = [[self._type],{}]
            if nuisanceMap[name][0][0] != self._type: raise RuntimeError, "Type mismatch for "+name
            nuisanceMap[name][1][process] = self._eff
        else: 
            nameN = self._nameN.format(process=process, channel=channel, jets=jets)
            if not nuisanceMap.has_key(nameN): nuisanceMap[nameN] = [["gmN", self._N],{}]
            if nuisanceMap[nameN][0][0] != "gmN":   raise RuntimeError, "Type mismatch for "+nameN
            if nuisanceMap[nameN][0][1] != self._N: raise RuntimeError, "N mismatch for "+nameN
            nuisanceMap[nameN][1][process] = self._alpha
            if self._hasErrAlpha:
                nameAlpha = self._nameAlpha.format(process=process, channel=channel, jets=jets)
                if not nuisanceMap.has_key(nameAlpha): nuisanceMap[nameAlpha] = [[self._typeAlpha],{}]
                if nuisanceMap[nameAlpha][0][0] != self._typeAlpha: raise RuntimeError, "Type mismatch for "+nameAlpha
                nuisanceMap[nameAlpha][1][process] = self._effAlpha

class DatacardWriter:
    def __init__(self,options):
        self.options = options
    def yieldsFromPlots(self,plots,channel,alwaysKeep=['data']):
        yields = {}
        for proc, ps in plots.iteritems():
            for cn, cp in ps:
                if cn != channel: continue
                val = cp.Integral()
                if val == 0 and proc not in alwaysKeep: continue
                yields[proc] = Yield(val)
        return yields
    def yieldsFromReport(self,report,channel,alwaysKeep=['data']):
        yields = {}
        for proc, ys in report.iteritems():
            for cn, cy, cye in ys[-1][1]:
                if cn != channel: continue
                if cy == 0 and proc not in alwaysKeep: continue
                yields[proc] = Yield(cy)
        return yields
    def loadDataDrivenYields(self, yields, process, file, mass, type, name, name2=None):
        map = file2map(file)
        if type == "gamma-lnN":
            (N, alpha, alphaErr) = map[mass]
            yields[process] = Yield(N, alpha, alphaErr, type="gamma", nameN=name, nameAlpha=name2)
        elif type == "gamma-gmM":
            (N, alpha, alphaErr) = map[mass]
            yields[process] = Yield(N, alpha, alphaErr, type="gamma", nameN=name, nameAlpha=name2, typeAlpha="gmM")
        elif type == "lnN":
            (y, dy) =  map[mass]
            yields[process] = Yield(y, dy, type="lnN", name=name)
    def loadDataDrivenYieldsDefault(self, yields, mass, channel, process):
        if mass < 200:
            if process == "WW": self.loadDataDrivenYieldsDefaultWW(yields, mass, channel)
        if "mumu" in channel or "elel" in channel:      
            if process == "DY": self.loadDataDrivenYieldsDefaultDY(yields, mass, channel)
        if process == "Top":    self.loadDataDrivenYieldsDefaultTop(yields, mass, channel)
        if process == "WJet":   self.loadDataDrivenYieldsDefaultWJet(yields, mass, channel)
    def loadDataDrivenYieldsDefaultWW(self, yields, mass, channel):
        file = "%s/WWCard_%s%s_%dj.txt" % (SYST_PATH, channel[0], channel[2], 0 if "0j" in channel else 1)
        self.loadDataDrivenYields(yields, "WW", file, mass, "gamma-lnN", name="CMS_ww_WW_stat", name2="CMS_ww_WW_extr")
    def loadDataDrivenYieldsDefaultDY(self, yields, mass, channel):
        DYproc = "DYMM" if "mumu" in channel else "DYEE"
        j = "0j" if "0j" in channel else "1j"
        file = "%s/DYCard_%s%s_%dj.txt" % (SYST_PATH, channel[0], channel[2], 0 if "0j" in channel else 1)
        self.loadDataDrivenYields(yields, DYproc, file, mass, "gamma-gmM", name="CMS_ww_"+DYproc+j+"_stat", name2="CMS_ww_"+DYproc+j+"_extr")
    def loadDataDrivenYieldsDefaultTop(self, yields, mass, channel):
        j = "0j" if "0j" in channel else "1j"
        file = "%s/TopCard_%s%s_%s.txt" % (SYST_PATH, channel[0], channel[2], j)
        self.loadDataDrivenYields(yields, "Top", file, mass, "gamma-lnN", name="CMS_ww_Top%s_stat" % j, name2="CMS_ww_Top%s_extr" % j)
    def loadDataDrivenYieldsDefaultWJet(self, yields, mass, channel):
        file = "%s/WJet_%s%s_%dj.txt" % (SYST_PATH, channel[0], channel[2], 0 if "0j" in channel else 1)
        self.loadDataDrivenYields(yields, "WJet", file, mass, "lnN", name="CMS_fake_%s" % channel[2])
    def writeFromYields(self, yields, nuisanceMap, fname, mass, channel, qqWWfromData, title="", shapesFile=None, signals=['ggH', 'vbfH']):
        """Yields must be in the form map (process -> Yield)"""
        ## Write datacard
        card = open(fname.format(mass=mass,channel=channel), "w")
        card.write("## H->WW for m(H) = %d, channel %s.\n" % (mass, channel))
        if title: card.write(" ##%s\n" % title)
        card.write("imax 1\njmax *\nkmax *\n");
        card.write(("-"*100) + "\n")
        if shapesFile:
            card.write("shapes  *         %s  %s $PROCESS/%s\n" % (channel, shapesFile.format(mass=mass,channel=channel), channel))
            card.write("shapes  data_obs  %s  %s     data/%s\n" % (channel, shapesFile.format(mass=mass,channel=channel), channel))
            card.write(("-"*100) + "\n")
        if self.options.asimov: card.write("# we put as 'observation' the expected background-only outcome\n")
        card.write("bin %s\n" % channel);
        if self.options.asimov and "data" not in yields:
            card.write("observation  %4d\n" % floor(0.5+sum([y.val for p,y in yields.items() if p != 'data' and p not in signals])))
        else:
            card.write("observation  %4d\n" % yields['data'].val)
        card.write(("-"*100) + "\n")
        keyline = [ ((-i if p in signals else i+1), p, y.val) for i,(p,y) in enumerate(yields.items()) if p != 'data']; 
        flen = min(len(channel), 6); fpad = " "*(flen-6)
        keyline.sort()
        card.write("bin                               " + "   ".join(   "%6s" % channel for i,p,y in keyline) + "\n");
        card.write("process                           " + "   ".join(fpad+"%6s"     % p for i,p,y in keyline) + "\n");
        card.write("process                           " + "   ".join(fpad+"   % 3d" % i for i,p,y in keyline) + "\n");
        card.write("rate                              " + "   ".join(fpad+"%6.3f"   % y for i,p,y in keyline) + "\n");
        card.write(("-"*100) + "\n")
        nuisSorted = nuisanceMap.keys(); nuisSorted.sort()
        for name in nuisSorted:
            (pdf,effect) = nuisanceMap[name]
            if len(pdf) == 1: card.write("%-25s %3s     " % (name,pdf[0]))
            else:             card.write("%-25s %3s %3d " % (name,pdf[0],pdf[1]))
            for i,p,y in keyline:
                if effect.has_key(p): 
                    if pdf[0] == 'gmN': card.write(fpad+"%6.4f" % effect[p])
                    else:               card.write(fpad+"  %4.2f" % effect[p])
                else:                   card.write(fpad+"     -")
                card.write("   ")
            card.write("\n");
    
