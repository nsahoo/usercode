from WWAnalysis.AnalysisStep.tree2yield import *
from optparse import OptionParser
parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
parser.add_option("-t", "--tree",  dest="tree", default='%sTree', help="Pattern for tree name");
parser.add_option("-l", "--lumi",   dest="lumi", type="float", default="1.0", help="Luminosity (in 1/fb)");
parser.add_option("-f", "--final",  dest="final", action="store_true", help="Just compute final yield after all cuts");
parser.add_option("--kmm", dest="kmm", type="float", nargs=2, default=(0.592288, 0.592288), help="kappa factor and uncertainty, for mu mu");
parser.add_option("--kee", dest="kee", type="float", nargs=2, default=(0.592288, 0.592288), help="kappa factor and uncertainty, for el el");
parser.add_option("-U", "--up-to-cut",   dest="upToCut",   type="string", help="Run selection only up to the cut matched by this regexp, included.") 
parser.add_option("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
parser.add_option("-I", "--invert-cut",  dest="cutsToInvert",  action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
parser.add_option("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
parser.add_option("-A", "--add-cut",     dest="cutsToAdd",     action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
(options, args) = parser.parse_args()
options.inclusive = False
options.nMinusOne = False
options.weight    = True
datasets = [ "SingleMuon2011A", "SingleMuon2011Av2", "DoubleMuon2011A", "DoubleMuon2011Av2",
             "DoubleElectron2011A", "DoubleElectron2011Av2", "MuEG2011A", "MuEG2011Av2" ]
mcsets = [ "DYtoMuMu", "DYtoElEl" ]
ttysData = [ TreeToYield("tree_%s.root" % D, options) for D in datasets ]
ttysMC   = [ TreeToYield("tree_%s.root" % D, options) for D in mcsets   ] 
cf0j = CutsFile(args[0],options)
cf1j = CutsFile(cf0j).replace('Jet veto', 'One jet', 'njet == 1 && dphilljet*sameflav < 165./180.*3.1415926').replace('B veto', 'B veto', 'bveto && nbjet == 0')
cfZ0j = CutsFile(cf0j).remove('m\(?ll\)?').remove('gMR').remove('dphi').replace('Z veto','Z window','abs(mll-91.1876)<15').remove('B veto').remove('m\(?ll\)? >|pT (min|max)|dphi|gMR').replace('pMET','pMET > 35', 'pmmet > 35')
cfZ1j = CutsFile(cf1j).remove('m\(?ll\)?').remove('gMR').remove('dphi').replace('Z veto','Z window','abs(mll-91.1876)<15').remove('B veto').remove('m\(?ll\)? >|pT (min|max)|dphi|gMR').replace('pMET','pMET > 35', 'pmmet > 35')
cfZ0jNoMET = CutsFile(cfZ0j).remove('pMET').remove('met')
## get yields in peak, with MET cut
reportDataZ0j = mergeReports([tty.getYields(cfZ0j) for tty in ttysData])
reportMCZ0j   = mergeReports([tty.getYields(cfZ0j) for tty in ttysMC])
reportDataZ1j = mergeReports([tty.getYields(cfZ1j) for tty in ttysData])
reportMCZ1j   = mergeReports([tty.getYields(cfZ1j) for tty in ttysMC])
print "\nZ peak in Data, WW level no btag, with MET cut, 0 jets"; ttysData[0].prettyPrint(reportDataZ0j)
print "\nZ peak in Sim., WW level no btag, with MET cut, 0 jets";   ttysMC[0].prettyPrint(reportMCZ0j)
print "\nZ peak in Data, WW level no btag, with MET cut, 1 jets"; ttysData[0].prettyPrint(reportDataZ1j)
print "\nZ peak in Sim., WW level no btag, with MET cut, 1 jets";   ttysMC[0].prettyPrint(reportMCZ1j)
Z0jmmData = reportDataZ0j[-1][1][MM][1];  Z0jeeData = reportDataZ0j[-1][1][EE][1]; Z0jxxData = (reportDataZ0j[-1][1][EM][1]+reportDataZ0j[-1][1][ME][1])
Z1jmmData = reportDataZ1j[-1][1][MM][1];  Z1jeeData = reportDataZ1j[-1][1][EE][1]; Z1jxxData = (reportDataZ1j[-1][1][EM][1]+reportDataZ1j[-1][1][ME][1])
Z0jmmMC = reportMCZ0j[-1][1][MM][1];  Z0jeeMC = reportMCZ0j[-1][1][EE][1]; Z0jxxMC = (reportMCZ0j[-1][1][EM][1]+reportMCZ0j[-1][1][ME][1])
Z1jmmMC = reportMCZ1j[-1][1][MM][1];  Z1jeeMC = reportMCZ1j[-1][1][EE][1]; Z1jxxMC = (reportMCZ1j[-1][1][EM][1]+reportMCZ1j[-1][1][ME][1])

reportDataZNoMET = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysData])
reportMCZNoMET   = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysMC])
print "\nZ peak in Data, WW level no btag, no MET, 0 jets"; ttysData[0].prettyPrint(reportDataZNoMET)
print "\nZ peak in Sim., WW level no btag, no MET, 0 jets";   ttysMC[0].prettyPrint(reportMCZNoMET)
keeData = sqrt(reportDataZNoMET[-1][1][EE][1]/float(reportDataZNoMET[-1][1][MM][1]));
kmmData = 1.0/keeData;
keeMC = sqrt(reportMCZNoMET[-1][1][EE][1]/float(reportMCZNoMET[-1][1][MM][1]));
kmmMC = 1.0/keeMC;
print "kee(Data) = ",keeData,"\tkee(MC) = ",keeMC
print "kmm(Data) = ",kmmData,"\tkmm(MC) = ",kmmMC

print "N(Zmm 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z0jmmData, kmmData, Z0jxxData, Z0jmmData - 0.5*kmmData*Z0jxxData, sqrt(Z0jxxData+pow(0.5*kmmData,2)*Z0jxxData))
print "N(Zmm 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z1jmmData, kmmData, Z1jxxData, Z1jmmData - 0.5*kmmData*Z1jxxData, sqrt(Z1jxxData+pow(0.5*kmmData,2)*Z1jxxData))
print "N(Zmm 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z0jmmMC, kmmMC, Z0jxxMC, Z0jmmMC - 0.5*kmmMC*Z0jxxMC, sqrt(Z0jxxMC+pow(0.5*kmmMC,2)*Z0jxxMC))
print "N(Zmm 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z1jmmMC, kmmMC, Z1jxxMC, Z1jmmMC - 0.5*kmmMC*Z1jxxMC, sqrt(Z1jxxMC+pow(0.5*kmmMC,2)*Z1jxxMC))
print "N(Zee 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z0jeeData, keeData, Z0jxxData, Z0jeeData - 0.5*keeData*Z0jxxData, sqrt(Z0jxxData+pow(0.5*keeData,2)*Z0jxxData))
print "N(Zee 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z1jeeData, keeData, Z1jxxData, Z1jeeData - 0.5*keeData*Z1jxxData, sqrt(Z1jxxData+pow(0.5*keeData,2)*Z1jxxData))
print "N(Zee 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z0jeeMC, keeMC, Z0jxxMC, Z0jeeMC - 0.5*keeMC*Z0jxxMC, sqrt(Z0jxxMC+pow(0.5*keeMC,2)*Z0jxxMC))
print "N(Zee 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %4.2f +/- %4.2f" % (Z1jeeMC, keeMC, Z1jxxMC, Z1jeeMC - 0.5*keeMC*Z1jxxMC, sqrt(Z1jxxMC+pow(0.5*keeMC,2)*Z1jxxMC))

#cfZ   = CutsFile(cf ).remove('m\(?ll\)?').remove('gMR').remove('dphi').replace('Z veto','Z window','abs(mll-91.1876)<15')
#cfZWW = CutsFile(cfWW).replace('Z veto','Z window','abs(mll-91.1876)<15')
#cfZNoMet = CutsFile(cfZ).remove("met").remove("pMET")
#cfNoMet   = CutsFile(cf ).remove("met").remove("pMET")
#cfWWNoMet = CutsFile(cfWW).remove("met").remove("pMET")
#reportDataZ    = mergeReports([tty.getYields(cfZ) for tty in ttysData])
#reportMCZ      = mergeReports([tty.getYields(cfZ) for tty in ttysMC])
#reportMCZNoMET = mergeReports([tty.getYields(cfZNoMet) for tty in ttysMC])
#reportMCNoMET  = mergeReports([tty.getYields(cfNoMet)  for tty in ttysMC])
#reportMC        = mergeReports([tty.getYields(cf)  for tty in ttysMC])
#reportMCWW      = mergeReports([tty.getYields(cfWW)  for tty in ttysMC])
#reportMCWWNoMET = mergeReports([tty.getYields(cfWWNoMet)  for tty in ttysMC])
#print "\nZ peak in data, with MET cut";    ttysData[0].prettyPrint(reportDataZ)
#print "\nZ peak in MC, with MET cut";    ttysMC[0].prettyPrint(reportMCZ)
#print "\nZ peak in MC, with no MET cut"; ttysMC[0].prettyPrint(reportMCZNoMET)
#print "\nStandard sel. in MC, with MET cut"; ttysMC[0].prettyPrint(reportMC)
#print "\nStandard sel. in MC, with no MET cut"; ttysMC[0].prettyPrint(reportMCNoMET)
#print "\nWW sel. in MC, with MET cut"; ttysMC[0].prettyPrint(reportMCWW)
#print "\nWW sel. in MC, with no MET cut"; ttysMC[0].prettyPrint(reportMCWWNoMET)
#Rmm = reportMCNoMET[-1][1][0][1]/reportMCZNoMET[-1][1][0][1]
#Ree = reportMCNoMET[-1][1][3][1]/reportMCZNoMET[-1][1][3][1]
#print "Rmm = ",Rmm
#print "Ree = ",Ree
#print "kmm = ",options.kmm[0],"+/-",options.kmm[1]
#print "kee = ",options.kee[0],"+/-",options.kee[1]
#print ""
#print "Estimate from data, mm = ",Rmm*(reportDataZ[-1][1][0][1] - options.kmm[0]*(reportDataZ[-1][1][1][1]+reportDataZ[-1][1][2][1]))
#print "    stat. uncertainty  = ",Rmm*sqrt(reportDataZ[-1][1][0][1] + (options.kmm[0])**2 * (reportDataZ[-1][1][1][1]+reportDataZ[-1][1][2][1]) )
#print "Estimate from data, ee = ",Ree*(reportDataZ[-1][1][2][1] - options.kee[0]*(reportDataZ[-1][1][1][1]+reportDataZ[-1][1][2][1]))
#print "    stat. uncertainty  = ",Ree*sqrt(reportDataZ[-1][1][2][1] + (options.kee[0])**2 * (reportDataZ[-1][1][1][1]+reportDataZ[-1][1][2][1]) )
