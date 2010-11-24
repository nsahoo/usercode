import FWCore.ParameterSet.Config as cms

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario 

choiceID=None; postfix=""
if len(args) > 1:
    choiceID = args[1]
    postfix  = "."+choiceID


process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TnP_EleID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    InputFileNames = cms.vstring('tnpZ_Data.root'),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTreeElEl"),
    OutputFileName = cms.string("data/TnP_Z_EleID_%s%s.root" % (scenario,postfix)),

    Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "70", "120", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("Probe #eta", "-3", "3", ""),
        phi    = cms.vstring("Probe #phi", "-3.2", "3.2", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "3", ""),
        nVtx = cms.vstring("N(prim. vtx)", "0", "6", ""),
        pair_nJet     = cms.vstring("N(jets)",        "0", "10", ""),
        pair_nJetNoPU = cms.vstring("N(jets w/o PU)", "0", "10", ""),
        run    = cms.vstring("Run number", "0", "9999999", ""),
    ),

    Categories = cms.PSet(
        passID = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        passIso = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
        HLT_Any1E  = cms.vstring("eih",  "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[2,1,3])",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),

    ),

    Efficiencies = cms.PSet(), # will be filled later
)

from WWAnalysis.AnalysisStep.tp_lepton_effs_cff import HLT1Es;
for X in HLT1Es:
    setattr(process.TnP_EleID.Categories,        X, cms.vstring("eih", "dummy[pass=1,fail=0]") )
    setattr(process.TnP_EleID.Categories, "tag_"+X, cms.vstring("eih", "dummy[pass=1,fail=0]") )

if scenario == "data_all":
    process.TnP_EleID.InputFileNames = [ "in/tnpLeptonID_ElectronRun2010B.root", "in/tnpLeptonID_EGSep17.root" ]

if scenario == "mc_all" or scenario == "mc_fast":
    process.TnP_EleID.InputFileNames = [ "in/tnpLeptonID_DYToEEM20CT10Z2powhegBX156.root", ]
if scenario == "mc_some" or scenario == "mc_some_fast":
    process.TnP_EleID.InputFileNames = [ "in/tnpLeptonID_DYToEEM20CT10Z2powhegBX156.root", ]

CONSTRAINTS = cms.PSet(
    pt        = cms.vdouble(20, 100),
    pair_nJet = cms.vdouble(-0.5,0.5)
)
AVERAGE = cms.PSet(
    CONSTRAINTS,
    abseta = cms.vdouble( 0,1.479, 2.5)
)
PT_ETA_BINS = AVERAGE.clone( 
    pt = cms.vdouble(10, 20, 30, 40, 50, 100)
)
VTX_BINS = AVERAGE.clone(
    nVtx = cms.vdouble(*[x-0.5 for x in range(8)])
)
JET_BINS = AVERAGE.clone(
    pair_nJet = cms.vdouble(*[x-0.5 for x in range(8)])
)

from math import pi
ETA_COARSE_BINS = cms.PSet(CONSTRAINTS, eta = cms.vdouble(-2.5,-1.479,1.479,2.5))
ETA_FINE_BINS   = cms.PSet(CONSTRAINTS, eta = cms.vdouble(-2.5,-2.0,-1.479,-0.7,0,0.7,1.479,2.0,2.5))
PHI_FINE_BINS   = cms.PSet(CONSTRAINTS, phi = cms.vdouble(*[pi*x/4. for x in range(-4,4+1)]))
RUN_BINS        = cms.PSet(CONSTRAINTS, run = cms.vdouble(136033,139980,141882,144114,147116,148058,149064,149442))

EFFS = {  #Name:NUM:DEN
   'all':(['passID','passIso'],[]), # Full WW efficiency
   'iso':(['passIso'],['passID']),  # Monitoring
   'id':(['passID'],['passIso']),  # Monitoring
   'isoNoID':(['passIso'],[]),  # Monitoring
   'idNoIso' :(['passID'], []),        # Monitoring
   'hlt':(['HLT_Any1E'],['passID','passIso']), # HLT, catch all
}
ALLBINS = [("pt_eta",PT_ETA_BINS),("average",AVERAGE),("vtx",VTX_BINS)]
ALLBINS += [("njet", JET_BINS)]
ALLBINS += [("eta", ETA_COARSE_BINS)]
# ALLBINS += [("eta", ETA_FINE_BINS), ("phi", PHI_FINE_BINS) ]
if "data" in scenario: ALLBINS += [ ("run", RUN_BINS) ]
for EN,E in EFFS.items(): 
    if choiceID != None and EN != choiceID: continue
    for  BN,BV in ALLBINS:
        # Num
        NUM = cms.vstring()
        for X in E[0]: NUM += [ X, "pass" ]
        # Den
        DEN = BV.clone()
        for X in E[1]: setattr(DEN, X, cms.vstring("pass"))
        # Shape
        S = "voigtPlusExpo" if "pt_eta" in BN else "vpvPlusExpo"
        # Eff
        setattr(process.TnP_EleID.Efficiencies, EN+"_"+BN, cms.PSet(
            EfficiencyCategoryAndState = NUM,
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = DEN,
            BinToPDFmap = cms.vstring(S)
        ))
        if "fast" in scenario: 
            del getattr(process.TnP_EleID.Efficiencies, EN+"_"+BN).BinToPDFmap

process.p = cms.Path(process.TnP_EleID)

