from WWAnalysis.AnalysisStep.tree2yield import *
from optparse import OptionParser
parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
addTreeToYieldOptions(parser)
parser.add_option("-0", "--0jetOnly",    dest="zorro", action="store_true", help="Don't try to create also 1 jet selection out of this cuts file.");
parser.add_option("--noVV", dest="doVV", action="store_false", default=True, help="Don't include WZ/ZZ in the MC estimate");
(options, args) = parser.parse_args()
options.inclusive = False
options.nMinusOne = False
options.weight    = True
do1 = False if options.zorro else True
datasets = [ "SingleMuon2011A", "SingleMuon2011Av2", "DoubleMuon2011A", "DoubleMuon2011Av2",
             "DoubleElectron2011A", "DoubleElectron2011Av2", "MuEG2011A", "MuEG2011Av2" ]
mcsets = [ "DYtoMuMu", "DYtoElEl", "DY10toMuMuZ2", "DY10toElElZ2" ];
if options.doVV: 
    mcsets += [ "WZtoAny", "ZZtoAny" ]
    print "Including WZ->any and ZZ->any in the MC samples."
ttysData = [ TreeToYield("tree_%s.root" % D, options) for D in datasets ]
ttysMC   = [ TreeToYield("tree_%s.root" % D, options) for D in mcsets   ] 
cf0j = CutsFile(args[0],options)
cf1j = CutsFile(cf0j).replace('jet veto', 'one jet', 'njet == 1 && dphilljet*sameflav < 165./180.*3.1415926').replace('top veto', 'top veto', 'bveto && nbjet == 0')
cf0jNoMET = CutsFile(cf0j).remove('met')
cf1jNoMET = CutsFile(cf1j).remove('met')
cfZ0j = CutsFile(cf0j).replace('Z veto','Z window','abs(mll-91.1876)<15').remove('gamma|Delta.phi|m_.*ell|p_T.*(min|max)').replace('35/20','pMET > 35', 'pmmet > 35')
cfZ1j = CutsFile(cf1j).replace('Z veto','Z window','abs(mll-91.1876)<15').remove('gamma|Delta.phi|m_.*ell|p_T.*(min|max)').replace('35/20','pMET > 35', 'pmmet > 35')
cfZ0jNoMET = CutsFile(cfZ0j).remove('met').remove('pMET')
cfZ1jNoMET = CutsFile(cfZ1j).remove('met').remove('pMET')
## get yields in peak, with MET cut
reportDataZ0j = mergeReports([tty.getYields(cfZ0j) for tty in ttysData])
reportMCZ0j   = mergeReports([tty.getYields(cfZ0j) for tty in ttysMC])
if do1: reportDataZ1j = mergeReports([tty.getYields(cfZ1j) for tty in ttysData])
if do1: reportMCZ1j   = mergeReports([tty.getYields(cfZ1j) for tty in ttysMC])
print "\nZ peak in Data, WW level no btag, with MET cut, 0 jets"; ttysData[0].prettyPrint(reportDataZ0j)
print "\nZ peak in Sim., WW level no btag, with MET cut, 0 jets";   ttysMC[0].prettyPrint(reportMCZ0j)
if do1:
    print "\nZ peak in Data, WW level no btag, with MET cut, 1 jets"; ttysData[0].prettyPrint(reportDataZ1j)
    print "\nZ peak in Sim., WW level no btag, with MET cut, 1 jets";   ttysMC[0].prettyPrint(reportMCZ1j)
Z0jmmData  = reportDataZ0j[-1][1][MM][1]; Z0jeeData  = reportDataZ0j[-1][1][EE][1]; Z0jxxData  = (reportDataZ0j[-1][1][EM][1]+reportDataZ0j[-1][1][ME][1])
Z0jmmMC    = reportMCZ0j[-1][1][MM][1];   Z0jeeMC    = reportMCZ0j[-1][1][EE][1];   Z0jxxMC    = (reportMCZ0j[-1][1][EM][1]+reportMCZ0j[-1][1][ME][1])
Z0jmmMCErr = reportMCZ0j[-1][1][MM][2];   Z0jeeMCErr = reportMCZ0j[-1][1][EE][2];   Z0jxxMCErr = hypot(reportMCZ0j[-1][1][EM][2], reportMCZ0j[-1][1][ME][2])
if do1:
    Z1jmmData  = reportDataZ1j[-1][1][MM][1]; Z1jeeData  = reportDataZ1j[-1][1][EE][1]; Z1jxxData  = (reportDataZ1j[-1][1][EM][1]+reportDataZ1j[-1][1][ME][1])
    Z1jmmMC    = reportMCZ1j[-1][1][MM][1];   Z1jeeMC    = reportMCZ1j[-1][1][EE][1];   Z1jxxMC    = (reportMCZ1j[-1][1][EM][1]+reportMCZ1j[-1][1][ME][1])
    Z1jmmMCErr = reportMCZ1j[-1][1][MM][2];   Z1jeeMCErr = reportMCZ1j[-1][1][EE][2];   Z1jxxMCErr = hypot(reportMCZ1j[-1][1][EM][2], reportMCZ1j[-1][1][ME][2])

## Now we compute events in the peak without MET cut, to get efficiencies
reportDataZNoMET = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysData])
reportMCZNoMET   = mergeReports([tty.getYields(cfZ0jNoMET) for tty in ttysMC])
print "\nZ peak in Data, WW level no btag, no MET, 0 jets"; ttysData[0].prettyPrint(reportDataZNoMET)
print "\nZ peak in Sim., WW level no btag, no MET, 0 jets";   ttysMC[0].prettyPrint(reportMCZNoMET)
keeData = sqrt(reportDataZNoMET[-1][1][EE][1]/float(reportDataZNoMET[-1][1][MM][1]));
kmmData = 1.0/keeData;
keeMC = sqrt(reportMCZNoMET[-1][1][EE][1]/float(reportMCZNoMET[-1][1][MM][1]));
kmmMC = 1.0/keeMC;

print ""
## Now we compute events in the signal region. we need the weight regularize the yields
avgweightMM  = ttysMC[0].getAverageWeight(cfZ0jNoMET.allCuts());
avgweightEE  = ttysMC[1].getAverageWeight(cfZ0jNoMET.allCuts());
avgweightSMM = ttysMC[0].getAverageWeight(cf0j.allCuts());
avgweightSEE = ttysMC[1].getAverageWeight(cf0j.allCuts());
print "Average event weight in Z peak w/o MET: ",avgweightMM, "(mm),\t",avgweightEE,"(ee)"
print "Average event weight in signal region : ",avgweightSMM,"(mm),\t",avgweightSEE,"(ee)"

## Do flavour subtraction to get DY-only yield
Z0jeeSubData = Z0jeeData - 0.5*keeData*Z0jxxData; Z0jeeSubDataErr =  sqrt(Z0jeeData+pow(0.5*keeData,2)*Z0jxxData)
Z0jeeSubMC   = Z0jeeMC   - 0.5*keeMC  *Z0jxxMC;   Z0jeeSubMCErr   = hypot(Z0jeeMCErr,0.5*keeMC*Z0jxxMCErr)
Z0jmmSubData = Z0jmmData - 0.5*kmmData*Z0jxxData; Z0jmmSubDataErr =  sqrt(Z0jmmData+pow(0.5*kmmData,2)*Z0jxxData)
Z0jmmSubMC   = Z0jmmMC   - 0.5*kmmMC  *Z0jxxMC;   Z0jmmSubMCErr   = hypot(Z0jmmMCErr,0.5*kmmMC*Z0jxxMCErr)
if do1:
    Z1jeeSubData = Z1jeeData - 0.5*keeData*Z1jxxData; Z1jeeSubDataErr =  sqrt(Z1jeeData+pow(0.5*keeData,2)*Z1jxxData)
    Z1jeeSubMC   = Z1jeeMC   - 0.5*keeMC  *Z1jxxMC;   Z1jeeSubMCErr   = hypot(Z1jeeMCErr,0.5*keeMC*Z1jxxMCErr)
    Z1jmmSubData = Z1jmmData - 0.5*kmmData*Z1jxxData; Z1jmmSubDataErr =  sqrt(Z1jmmData+pow(0.5*kmmData,2)*Z1jxxData)
    Z1jmmSubMC   = Z1jmmMC   - 0.5*kmmMC  *Z1jxxMC;   Z1jmmSubMCErr   = hypot(Z1jmmMCErr,0.5*kmmMC*Z1jxxMCErr)

print ""
## Now we compute events in the signal region. we need the weight regularize the yields
avgweightMM  = ttysMC[0].getAverageWeight(cfZ0jNoMET.allCuts());
avgweightEE  = ttysMC[1].getAverageWeight(cfZ0jNoMET.allCuts());
#avgweightSMM = ttysMC[0].getAverageWeight(cf0j.allCuts()); # Gives problems at H->WW level
#avgweightSEE = ttysMC[1].getAverageWeight(cf0j.allCuts()); #
avgweightSMM = avgweightMM
avgweightSEE = avgweightEE
print "Average event weight in Z peak w/o MET: ",avgweightMM, "(mm),\t",avgweightEE,"(ee)"
print "Average event weight in signal region : ",avgweightSMM,"(mm),\t",avgweightSEE,"(ee)"

report0jMC = mergeReports([tty.getYields(cf0j) for tty in ttysMC])
print "\nFull selection in Sim., 0 jets";   ttysMC[0].prettyPrint(report0jMC)
if do1:
    report1jMC = mergeReports([tty.getYields(cf1j) for tty in ttysMC])
    print "\nFull selection in Sim., 1 jets";   ttysMC[0].prettyPrint(report1jMC)

## Then we try to do the same but as a function of the MET cut
plot0j  = mergePlots([tty.getPlots("pmmet_0j","pmmet", "50,0.,50.", cf0jNoMET.allCuts()) for tty in ttysMC])
plotZ0j = mergePlots([tty.getPlots("pmmetZ_0j","pmmet","50,0.,50.",cfZ0jNoMET.allCuts()) for tty in ttysMC])
totZmm0j = plotZ0j[MM][1].Integral(1,51); totSmm0j = plot0j[MM][1].Integral(1,51);
totZee0j = plotZ0j[EE][1].Integral(1,51); totSee0j = plot0j[EE][1].Integral(1,51);
totZme0j = plotZ0j[ME][1].Integral(1,51); totSme0j = plot0j[ME][1].Integral(1,51);
totZem0j = plotZ0j[EM][1].Integral(1,51); totSem0j = plot0j[ME][1].Integral(1,51);
plotData0j  = mergePlots([tty.getPlots("pmmetD_0j","pmmet", "50,0.,50.", cf0jNoMET.allCuts()) for tty in ttysData])
plotDataZ0j = mergePlots([tty.getPlots("pmmetDZ_0j","pmmet","50,0.,50.",cfZ0jNoMET.allCuts()) for tty in ttysData])
totDataZee0j = plotDataZ0j[EE][1].Integral(1,51); totDataSee0j = plotData0j[EE][1].Integral(1,51);
totDataZem0j = plotDataZ0j[EM][1].Integral(1,51); totDataSem0j = plotData0j[EM][1].Integral(1,51);
totDataZme0j = plotDataZ0j[ME][1].Integral(1,51); totDataSme0j = plotData0j[ME][1].Integral(1,51); 
totDataZmm0j = plotDataZ0j[MM][1].Integral(1,51); totDataSmm0j = plotData0j[MM][1].Integral(1,51); 
if do1:
    plot1j  = mergePlots([tty.getPlots("pmmet_1j","pmmet", "50,0.,50.", cf1jNoMET.allCuts()) for tty in ttysMC])
    plotZ1j = mergePlots([tty.getPlots("pmmetZ_1j","pmmet","50,0.,50.",cfZ1jNoMET.allCuts()) for tty in ttysMC])
    totZmm1j = plotZ1j[MM][1].Integral(1,51); totSmm1j = plot1j[MM][1].Integral(1,51);
    totZme1j = plotZ1j[ME][1].Integral(1,51); totSme1j = plot1j[ME][1].Integral(1,51);
    totZem1j = plotZ1j[EM][1].Integral(1,51); totSem1j = plot1j[EM][1].Integral(1,51);
    totZee1j = plotZ1j[EE][1].Integral(1,51); totSee1j = plot1j[EE][1].Integral(1,51);
    plotData1j  = mergePlots([tty.getPlots("pmmetD_1j","pmmet", "50,0.,50.", cf1jNoMET.allCuts()) for tty in ttysData])
    plotDataZ1j = mergePlots([tty.getPlots("pmmetDZ_1j","pmmet","50,0.,50.",cfZ1jNoMET.allCuts()) for tty in ttysData])
    totDataZee1j = plotDataZ1j[EE][1].Integral(1,51); totDataSee1j = plotData1j[EE][1].Integral(1,51);
    totDataZem1j = plotDataZ1j[EM][1].Integral(1,51); totDataSem1j = plotData1j[EM][1].Integral(1,51);
    totDataZme1j = plotDataZ1j[ME][1].Integral(1,51); totDataSme1j = plotData1j[ME][1].Integral(1,51);
    totDataZmm1j = plotDataZ1j[MM][1].Integral(1,51); totDataSmm1j = plotData1j[MM][1].Integral(1,51);

## Compute extrapolation factors for nominal pMET cut
alpha0jee = (report0jMC[-1][1][EE][1]+avgweightSEE)/(reportMCZ0j[-1][1][EE][1]+avgweightEE);
alpha0jeeErr = min(1.0, sqrt(1.0/(report0jMC[-1][1][EE][1]/avgweightSEE+1) + 1.0/(reportMCZ0j[-1][1][EE][1]/avgweightEE+1)));
alpha0jmm = (report0jMC[-1][1][MM][1]+avgweightSMM)/(reportMCZ0j[-1][1][MM][1]+avgweightMM);
alpha0jmmErr = min(1.0, sqrt(1.0/(report0jMC[-1][1][MM][1]/avgweightSMM+1) + 1.0/(reportMCZ0j[-1][1][MM][1]/avgweightMM+1)));
if do1:
    alpha1jee = (report1jMC[-1][1][EE][1]+avgweightSEE)/(reportMCZ1j[-1][1][EE][1]+avgweightEE);
    alpha1jeeErr = min(1.0, sqrt(1.0/(report1jMC[-1][1][EE][1]/avgweightSEE+1) + 1.0/(reportMCZ1j[-1][1][EE][1]/avgweightEE+1)));
    alpha1jmm = (report1jMC[-1][1][MM][1]+avgweightSMM)/(reportMCZ1j[-1][1][MM][1]+avgweightMM);
    alpha1jmmErr = min(1.0, sqrt(1.0/(report1jMC[-1][1][MM][1]/avgweightSMM+1) + 1.0/(reportMCZ1j[-1][1][MM][1]/avgweightMM+1)));

## Let's do a scan of alpha vs pMET cut
tfout = ROOT.TFile("dyEstimate.root", "RECREATE")
print "\nEstimate of alpha vs pMET, regularized, 0 jets. Data is flavour subtracted."
print "\tpMET  ------------ MC -------------  |  ------------- DATA ----------  |  ------------ MC -------------  |  ------------- DATA ----------"
print "\tpMET    Zmm0j    Smm0j   alpha  err  |    Zmm0j    Smm0j   alpha  err  |    Zee0j    See0j   alpha  err  |    Zee0j    See0j   alpha  err"
gZmm0jMC   = ROOT.TGraphAsymmErrors(50); gZee0jMC   = ROOT.TGraphAsymmErrors(50); gSmm0jMC   = ROOT.TGraphAsymmErrors(50); gSee0jMC   = ROOT.TGraphAsymmErrors(50);
gZmm0jData = ROOT.TGraphAsymmErrors(50); gZee0jData = ROOT.TGraphAsymmErrors(50); gSmm0jData = ROOT.TGraphAsymmErrors(50); gSee0jData = ROOT.TGraphAsymmErrors(50);
gAmm0jMC   = ROOT.TGraphAsymmErrors(50); gAee0jMC   = ROOT.TGraphAsymmErrors(50); gAmm0jData = ROOT.TGraphAsymmErrors(50); gAee0jData = ROOT.TGraphAsymmErrors(50)
for b in range(50):
    totDataSmm0jSub = max(0, totDataSmm0j - 0.5*kmmData * (totDataSem0j + totDataSme0j))
    totDataSee0jSub = max(0, totDataSee0j - 0.5*keeData * (totDataSem0j + totDataSme0j))
    totDataSmm0jErr = sqrt(totDataSmm0j + pow(0.5*kmmData,2) * (totDataSem0j + totDataSme0j))/totDataSmm0jSub if totDataSmm0jSub > 0 else 1.0
    totDataSee0jErr = sqrt(totDataSee0j + pow(0.5*keeData,2) * (totDataSem0j + totDataSme0j))/totDataSee0jSub if totDataSee0jSub > 0 else 1.0
    totDataZmm0jSub = max(0, totDataZmm0j - 0.5*kmmData * (totDataZem0j + totDataZme0j))
    totDataZee0jSub = max(0, totDataZee0j - 0.5*keeData * (totDataZem0j + totDataZme0j))
    totDataZmm0jErr = sqrt(totDataZmm0j + pow(0.5*kmmData,2) * (totDataZem0j + totDataZme0j))/totDataZmm0jSub if totDataZmm0jSub > 0 else 1.0
    totDataZee0jErr = sqrt(totDataZee0j + pow(0.5*keeData,2) * (totDataZem0j + totDataZme0j))/totDataZee0jSub if totDataZee0jSub > 0 else 1.0
    Aee0jData = (totDataSee0jSub+1)/(totDataZee0jSub+1);   Aee0jDataErr = sqrt(1.0/(totDataZee0jSub+1) + 1.0/(totDataSee0jErr+1));
    Amm0jData = (totDataSmm0jSub+1)/(totDataZmm0jSub+1);   Amm0jDataErr = sqrt(1.0/(totDataZmm0jSub+1) + 1.0/(totDataSmm0jErr+1));
    totSmm0jSub = max(0, totSmm0j - 0.5*kmmMC * (totSem0j + totSme0j))
    totSee0jSub = max(0, totSee0j - 0.5*keeMC * (totSem0j + totSme0j))
    totSmm0jErr = sqrt(totSmm0j + pow(0.5*kmmMC,2) * (totSem0j + totSme0j))/totSmm0jSub if totSmm0jSub > 0 else 1.0
    totSee0jErr = sqrt(totSee0j + pow(0.5*keeMC,2) * (totSem0j + totSme0j))/totSee0jSub if totSee0jSub > 0 else 1.0
    totZmm0jSub = max(0, totZmm0j - 0.5*kmmMC * (totZem0j + totZme0j))
    totZee0jSub = max(0, totZee0j - 0.5*keeMC * (totZem0j + totZme0j))
    totZmm0jErr = sqrt(totZmm0j + pow(0.5*kmmMC,2) * (totZem0j + totZme0j))/totZmm0jSub if totZmm0jSub > 0 else 1.0
    totZee0jErr = sqrt(totZee0j + pow(0.5*keeMC,2) * (totZem0j + totZme0j))/totZee0jSub if totZee0jSub > 0 else 1.0
    Aee0j = (totSee0jSub+avgweightSEE)/(totZee0jSub+avgweightEE);   Aee0jErr = sqrt(1.0/(totZee0jSub/avgweightEE+1) + 1.0/(totSee0jErr/avgweightSEE+1));
    Amm0j = (totSmm0jSub+avgweightSMM)/(totZmm0jSub+avgweightMM);   Amm0jErr = sqrt(1.0/(totZmm0jSub/avgweightEE+1) + 1.0/(totSmm0jErr/avgweightSMM+1));
    gZmm0jMC.SetPoint(b, 1.0*b, totZmm0jSub); gZmm0jMC.SetPointError(b, 0, 0, totZmm0jErr*totZmm0jSub, totZmm0jErr*totZmm0jSub);
    gZee0jMC.SetPoint(b, 1.0*b, totZee0jSub); gZee0jMC.SetPointError(b, 0, 0, totZee0jErr*totZee0jSub, totZee0jErr*totZee0jSub);
    gSmm0jMC.SetPoint(b, 1.0*b, totSmm0jSub); gSmm0jMC.SetPointError(b, 0, 0, totSmm0jErr*totSmm0jSub, totSmm0jErr*totSmm0jSub);
    gSee0jMC.SetPoint(b, 1.0*b, totSee0jSub); gSee0jMC.SetPointError(b, 0, 0, totSee0jErr*totSee0jSub, totSee0jErr*totSee0jSub);
    gAmm0jMC.SetPoint(b, 1.0*b,    Amm0j);    gAmm0jMC.SetPointError(b, 0, 0, Amm0j*Amm0jErr, Amm0j*Amm0jErr);
    gAee0jMC.SetPoint(b, 1.0*b,    Aee0j);    gAee0jMC.SetPointError(b, 0, 0, Aee0j*Aee0jErr, Aee0j*Aee0jErr);
    gZmm0jData.SetPoint(b, 1.0*b, totDataZmm0jSub); gZmm0jData.SetPointError(b, 0, 0, totDataZmm0jErr, totDataZmm0jErr);
    gZee0jData.SetPoint(b, 1.0*b, totDataZee0jSub); gZee0jData.SetPointError(b, 0, 0, totDataZee0jErr, totDataZee0jErr);
    gSmm0jData.SetPoint(b, 1.0*b, totDataSmm0jSub); gSmm0jData.SetPointError(b, 0, 0, totDataSmm0jErr, totDataSmm0jErr);
    gSee0jData.SetPoint(b, 1.0*b, totDataSee0jSub); gSee0jData.SetPointError(b, 0, 0, totDataSee0jErr, totDataSee0jErr);
    gAmm0jData.SetPoint(b, 1.0*b,       Amm0jData); gAmm0jData.SetPointError(b, 0, 0, Amm0jData*Amm0jDataErr, Amm0jData*Amm0jDataErr);
    gAee0jData.SetPoint(b, 1.0*b,       Aee0jData); gAee0jData.SetPointError(b, 0, 0, Aee0jData*Aee0jDataErr, Aee0jData*Aee0jDataErr);
    print "\t %2d   %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f" % (
            b,
            totZmm0jSub, totSmm0jSub, Amm0j, min(1.0, Amm0jErr),
            totDataZmm0jSub, totDataSmm0jSub, Amm0jData, min(1.0, Amm0jDataErr),
            totZee0jSub, totSee0jSub, Aee0j, min(1.0, Aee0jErr),
            totDataZee0jSub, totDataSee0jSub, Aee0jData, min(1.0, Aee0jDataErr),
          )
    totZmm0j -= plotZ0j[MM][1].GetBinContent(b+1); totSmm0j -= plot0j[MM][1].GetBinContent(b+1)
    totZee0j -= plotZ0j[EE][1].GetBinContent(b+1); totSee0j -= plot0j[EE][1].GetBinContent(b+1)
    totZem0j -= plotZ0j[EM][1].GetBinContent(b+1); totSem0j -= plot0j[EM][1].GetBinContent(b+1)
    totZme0j -= plotZ0j[ME][1].GetBinContent(b+1); totSme0j -= plot0j[ME][1].GetBinContent(b+1)
    totDataZmm0j -= plotDataZ0j[MM][1].GetBinContent(b+1); totDataSmm0j -= plotData0j[MM][1].GetBinContent(b+1)
    totDataZee0j -= plotDataZ0j[EE][1].GetBinContent(b+1); totDataSee0j -= plotData0j[EE][1].GetBinContent(b+1)
    totDataZem0j -= plotDataZ0j[EM][1].GetBinContent(b+1); totDataSem0j -= plotData0j[EM][1].GetBinContent(b+1)
    totDataZme0j -= plotDataZ0j[ME][1].GetBinContent(b+1); totDataSme0j -= plotData0j[ME][1].GetBinContent(b+1)
gZmm0jMC.SetName("gZmm0jMC"); gSmm0jMC.SetName("gSmm0jMC"); gAmm0jMC.SetName("gAmm0jMC"); 
gZee0jMC.SetName("gZee0jMC"); gSee0jMC.SetName("gSee0jMC"); gAee0jMC.SetName("gAee0jMC");
gZmm0jData.SetName("gZmm0jData"); gSmm0jData.SetName("gSmm0jData"); gAmm0jData.SetName("gAmm0jData"); 
gZee0jData.SetName("gZee0jData"); gSee0jData.SetName("gSee0jData"); gAee0jData.SetName("gAee0jData");
gZmm0jMC.Write(); gZee0jMC.Write(); gSmm0jMC.Write(); gSee0jMC.Write(); gAmm0jMC.Write(); gAee0jMC.Write();
gZmm0jData.Write(); gZee0jData.Write(); gSmm0jData.Write(); gSee0jData.Write(); gAmm0jData.Write(); gAee0jData.Write();

if do1:
    gZmm1jMC   = ROOT.TGraphAsymmErrors(50); gZee1jMC   = ROOT.TGraphAsymmErrors(50); gSmm1jMC   = ROOT.TGraphAsymmErrors(50); gSee1jMC   = ROOT.TGraphAsymmErrors(50);
    gZmm1jData = ROOT.TGraphAsymmErrors(50); gZee1jData = ROOT.TGraphAsymmErrors(50); gSmm1jData = ROOT.TGraphAsymmErrors(50); gSee1jData = ROOT.TGraphAsymmErrors(50);
    gAmm1jMC   = ROOT.TGraphAsymmErrors(50); gAee1jMC   = ROOT.TGraphAsymmErrors(50); gAmm1jData = ROOT.TGraphAsymmErrors(50); gAee1jData = ROOT.TGraphAsymmErrors(50)
    print "\nEstimate of alpha vs pMET, regularized, 1 jet. Data is flavour subtracted."
    print "\tpMET  ------------ MC -------------  |  ------------- DATA ----------  |  ------------ MC -------------  |  ------------- DATA ----------"
    print "\tpMET    Zmm1j    Smm1j   alpha  err  |    Zmm1j    Smm1j   alpha  err  |    Zee1j    See1j   alpha  err  |    Zee1j    See1j   alpha  err"
    for b in range(50):
        totDataSmm1jSub = max(0, totDataSmm1j - 0.5*kmmData * (totDataSem1j + totDataSme1j))
        totDataSee1jSub = max(0, totDataSee1j - 0.5*keeData * (totDataSem1j + totDataSme1j))
        totDataSmm1jErr = sqrt(totDataSmm1j + pow(0.5*kmmData,2) * (totDataSem1j + totDataSme1j))/totDataSmm1jSub if totDataSmm1jSub > 0 else 1.0
        totDataSee1jErr = sqrt(totDataSee1j + pow(0.5*keeData,2) * (totDataSem1j + totDataSme1j))/totDataSee1jSub if totDataSee1jSub > 0 else 1.0
        totDataZmm1jSub = max(0, totDataZmm1j - 0.5*kmmData * (totDataZem1j + totDataZme1j))
        totDataZee1jSub = max(0, totDataZee1j - 0.5*keeData * (totDataZem1j + totDataZme1j))
        totDataZmm1jErr = sqrt(totDataZmm1j + pow(0.5*kmmData,2) * (totDataZem1j + totDataZme1j))/totDataZmm1jSub if totDataZmm1jSub > 0 else 1.0
        totDataZee1jErr = sqrt(totDataZee1j + pow(0.5*keeData,2) * (totDataZem1j + totDataZme1j))/totDataZee1jSub if totDataZee1jSub > 0 else 1.0
        Aee1jData = (totDataSee1jSub+1)/(totDataZee1jSub+1);   Aee1jDataErr = sqrt(1.0/(totDataZee1jSub+1) + 1.0/(totDataSee1jErr+1));
        Amm1jData = (totDataSmm1jSub+1)/(totDataZmm1jSub+1);   Amm1jDataErr = sqrt(1.0/(totDataZmm1jSub+1) + 1.0/(totDataSmm1jErr+1));
        totSmm1jSub = max(0, totSmm1j - 0.5*kmmMC * (totSem1j + totSme1j))
        totSee1jSub = max(0, totSee1j - 0.5*keeMC * (totSem1j + totSme1j))
        totSmm1jErr = sqrt(totSmm1j + pow(0.5*kmmMC,2) * (totSem1j + totSme1j))/totSmm1jSub if totSmm1jSub > 0 else 1.0
        totSee1jErr = sqrt(totSee1j + pow(0.5*keeMC,2) * (totSem1j + totSme1j))/totSee1jSub if totSee1jSub > 0 else 1.0
        totZmm1jSub = max(0, totZmm1j - 0.5*kmmMC * (totZem1j + totZme1j))
        totZee1jSub = max(0, totZee1j - 0.5*keeMC * (totZem1j + totZme1j))
        totZmm1jErr = sqrt(totZmm1j + pow(0.5*kmmMC,2) * (totZem1j + totZme1j))/totZmm1jSub if totZmm1jSub > 0 else 1.0
        totZee1jErr = sqrt(totZee1j + pow(0.5*keeMC,2) * (totZem1j + totZme1j))/totZee1jSub if totZee1jSub > 0 else 1.0
        Aee1j = (totSee1jSub+avgweightSEE)/(totZee1jSub+avgweightEE);   Aee1jErr = sqrt(1.0/(totZee1jSub/avgweightEE+1) + 1.0/(totSee1jErr/avgweightSEE+1));
        Amm1j = (totSmm1jSub+avgweightSMM)/(totZmm1jSub+avgweightMM);   Amm1jErr = sqrt(1.0/(totZmm1jSub/avgweightEE+1) + 1.0/(totSmm1jErr/avgweightSMM+1));
        gZmm1jMC.SetPoint(b, 1.0*b, totZmm1jSub); gZmm1jMC.SetPointError(b, 0, 0, totZmm1jErr*totZmm1jSub, totZmm1jErr*totZmm1jSub);
        gZee1jMC.SetPoint(b, 1.0*b, totZee1jSub); gZee1jMC.SetPointError(b, 0, 0, totZee1jErr*totZee1jSub, totZee1jErr*totZee1jSub);
        gSmm1jMC.SetPoint(b, 1.0*b, totSmm1jSub); gSmm1jMC.SetPointError(b, 0, 0, totSmm1jErr*totSmm1jSub, totSmm1jErr*totSmm1jSub);
        gSee1jMC.SetPoint(b, 1.0*b, totSee1jSub); gSee1jMC.SetPointError(b, 0, 0, totSee1jErr*totSee1jSub, totSee1jErr*totSee1jSub);
        gAmm1jMC.SetPoint(b, 1.0*b,    Amm1j);    gAmm1jMC.SetPointError(b, 0, 0, Amm1j*Amm1jErr, Amm1j*Amm1jErr);
        gAee1jMC.SetPoint(b, 1.0*b,    Aee1j);    gAee1jMC.SetPointError(b, 0, 0, Aee1j*Aee1jErr, Aee1j*Aee1jErr);
        gZmm1jData.SetPoint(b, 1.0*b, totDataZmm1jSub); gZmm1jData.SetPointError(b, 0, 0, totDataZmm1jErr, totDataZmm1jErr);
        gZee1jData.SetPoint(b, 1.0*b, totDataZee1jSub); gZee1jData.SetPointError(b, 0, 0, totDataZee1jErr, totDataZee1jErr);
        gSmm1jData.SetPoint(b, 1.0*b, totDataSmm1jSub); gSmm1jData.SetPointError(b, 0, 0, totDataSmm1jErr, totDataSmm1jErr);
        gSee1jData.SetPoint(b, 1.0*b, totDataSee1jSub); gSee1jData.SetPointError(b, 0, 0, totDataSee1jErr, totDataSee1jErr);
        gAmm1jData.SetPoint(b, 1.0*b,       Amm1jData); gAmm1jData.SetPointError(b, 0, 0, Amm1jData*Amm1jDataErr, Amm1jData*Amm1jDataErr);
        gAee1jData.SetPoint(b, 1.0*b,       Aee1jData); gAee1jData.SetPointError(b, 0, 0, Aee1jData*Aee1jDataErr, Aee1jData*Aee1jDataErr);
        print "\t %2d   %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f |  %9.1f %8.2f  %5.3f %4.2f" % (
                b,
                totZmm1jSub, totSmm1jSub, Amm1j, min(1.0, Amm1jErr),
                totDataZmm1jSub, totDataSmm1jSub, Amm1jData, min(1.0, Amm1jDataErr),
                totZee1jSub, totSee1jSub, Aee1j, min(1.0, Aee1jErr),
                totDataZee1jSub, totDataSee1jSub, Aee1jData, min(1.0, Aee1jDataErr),
              )
        totZmm1j -= plotZ1j[MM][1].GetBinContent(b+1); totSmm1j -= plot1j[MM][1].GetBinContent(b+1)
        totZee1j -= plotZ1j[EE][1].GetBinContent(b+1); totSee1j -= plot1j[EE][1].GetBinContent(b+1)
        totZem1j -= plotZ1j[EM][1].GetBinContent(b+1); totSem1j -= plot1j[EM][1].GetBinContent(b+1)
        totZme1j -= plotZ1j[ME][1].GetBinContent(b+1); totSme1j -= plot1j[ME][1].GetBinContent(b+1)
        totDataZmm1j -= plotDataZ1j[MM][1].GetBinContent(b+1); totDataSmm1j -= plotData1j[MM][1].GetBinContent(b+1)
        totDataZee1j -= plotDataZ1j[EE][1].GetBinContent(b+1); totDataSee1j -= plotData1j[EE][1].GetBinContent(b+1)
        totDataZem1j -= plotDataZ1j[EM][1].GetBinContent(b+1); totDataSem1j -= plotData1j[EM][1].GetBinContent(b+1)
        totDataZme1j -= plotDataZ1j[ME][1].GetBinContent(b+1); totDataSme1j -= plotData1j[ME][1].GetBinContent(b+1)
    gZmm1jMC.SetName("gZmm1jMC"); gSmm1jMC.SetName("gSmm1jMC"); gAmm1jMC.SetName("gAmm1jMC"); 
    gZee1jMC.SetName("gZee1jMC"); gSee1jMC.SetName("gSee1jMC"); gAee1jMC.SetName("gAee1jMC");
    gZmm1jData.SetName("gZmm1jData"); gSmm1jData.SetName("gSmm1jData"); gAmm1jData.SetName("gAmm1jData"); 
    gZee1jData.SetName("gZee1jData"); gSee1jData.SetName("gSee1jData"); gAee1jData.SetName("gAee1jData");
    gZmm1jMC.Write(); gZee1jMC.Write(); gSmm1jMC.Write(); gSee1jMC.Write(); gAmm1jMC.Write(); gAee1jMC.Write();
    gZmm1jData.Write(); gZee1jData.Write(); gSmm1jData.Write(); gSee1jData.Write(); gAmm1jData.Write(); gAee1jData.Write();

bg0jmm = alpha0jmm*Z0jmmSubData; bg0jee = alpha0jee*Z0jeeSubData;
if do1:
    bg1jee = alpha1jee*Z1jeeSubData; bg1jmm = alpha1jmm*Z1jmmSubData;

print "\n"
print "\n --- Efficiency factors, from Z peak with no MET cut.  --- "
print "kee(Data) = sqrt(N(ee)/N(mm)) = %6.4f\tkee(MC) = %6.4f" % (keeData ,keeMC)
print "kmm(Data) = sqrt(N(mm)/N(ee)) = %6.4f\tkmm(MC) = %6.4f" % (kmmData ,kmmMC)

print "\n --- Yields in the peak, with MET cut                  --- "
print "N(Zmm 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jmmData, kmmData, Z0jxxData, Z0jmmSubData, Z0jmmSubDataErr)
print "N(Zmm 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jmmMC  , kmmMC  , Z0jxxMC  , Z0jmmSubMC  , Z0jmmSubMCErr  )
print "N(Zee 0j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jeeData, keeData, Z0jxxData, Z0jeeSubData, Z0jeeSubDataErr)
print "N(Zee 0j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z0jeeMC  , keeMC  , Z0jxxMC  , Z0jeeSubMC  , Z0jeeSubMCErr  )

if do1:
    print "N(Zmm 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jmmData, kmmData, Z1jxxData, Z1jmmSubData, Z1jmmSubDataErr)
    print "N(Zmm 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jmmMC  , kmmMC  , Z1jxxMC  , Z1jmmSubMC  , Z1jmmSubMCErr  )
    print "N(Zee 1j, Data) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jeeData, keeData, Z1jxxData, Z1jeeSubData, Z1jeeSubDataErr)
    print "N(Zee 1j, Sim.) = (%5.1f - %5.3f/2 * %5.1f) = %5.2f +/- %4.2f" % (Z1jeeMC  , keeMC  , Z1jxxMC  , Z1jeeSubMC  , Z1jeeSubMCErr  )


print "\n --- Extrapolation factors to final selection, from MC --- "
print "alpha mm, 0j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
         report0jMC[-1][1][MM][1], Z0jmmMC, report0jMC[-1][1][MM][1]/avgweightMM, Z0jmmMC/avgweightMM)
print "alpha ee, 0j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
         report0jMC[-1][1][EE][1], Z0jeeMC, report0jMC[-1][1][EE][1]/avgweightEE, Z0jeeMC/avgweightEE)
if do1:
    print "alpha mm, 1j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
             report1jMC[-1][1][MM][1], Z1jmmMC, report1jMC[-1][1][MM][1]/avgweightMM, Z1jmmMC/avgweightMM)
    print "alpha ee, 1j, bare: %6.2f / %6.2f = %6.2f / %6.2f (effective events) " % (
             report1jMC[-1][1][EE][1], Z1jeeMC, report1jMC[-1][1][EE][1]/avgweightEE, Z1jeeMC/avgweightEE)
print "alpha mm, 0j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha0jmm, alpha0jmmErr*100 )
print "alpha ee, 0j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha0jee, alpha0jeeErr*100 )
if do1:
    print "alpha ee, 1j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha1jee, alpha1jeeErr*100 )
    print "alpha mm, 1j, regularized: %7.4f   (rel. uncertainty %6.2f%%)" % ( alpha1jmm, alpha1jmmErr*100 )

print "\n --- Estimated background in the final selection       --- "
print "N(Z->mm, all cuts, 0j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg0jmm, alpha0jmm*Z0jmmSubDataErr, bg0jmm*alpha0jmmErr)
print "N(Z->ee, all cuts, 0j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg0jee, alpha0jee*Z0jeeSubDataErr, bg0jee*alpha0jeeErr)
if do1:
    print "N(Z->mm, all cuts, 1j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg1jmm, alpha1jmm*Z1jmmSubDataErr, bg1jmm*alpha1jmmErr)
    print "N(Z->ee, all cuts, 1j) = %7.3f +/- %7.3f (stat. data) +/- %7.3f (stat. mc)" % (bg1jee, alpha1jee*Z1jeeSubDataErr, bg1jee*alpha1jeeErr)

