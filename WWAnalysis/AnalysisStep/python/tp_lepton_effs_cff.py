import FWCore.ParameterSet.Config as cms

from WWAnalysis.AnalysisStep.wwMuons_cfi     import MUON_ISO_CUT, MUON_ID_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_ISO_CUT,  ELE_ID_CUT

HLT1Es = [
    "HLT_Ele10_LW_L1R",                     # MC && [136033,139980]
    "HLT_Ele15_SW_L1R",                     # MC && [140058,141882]
    "HLT_Ele15_SW_CaloEleId_L1R", 	    # [141956,144114] (~3/pb)
    "HLT_Ele17_SW_CaloEleId_L1R", 	    # [146428,147116] (~7/pb)
    "HLT_Ele17_SW_TightEleId_L1R", 	    # [147196,148058]
    "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2", # [148819,149064]
    "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3", # [149181,149442] 
]
#HLT2Es = [
#    "HLT_DoubleEle10_SW_L1R",                 # 2E31
#    "HLT_DoubleEle15_SW_L1R",                 # 6E31
#    "HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R",  # 6E31
#    "HLT_DoubleEle17_SW_L1R",                 # 2E32
#    "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R", # 2E32
#]
HLT1Ms = [
    "HLT_Mu9",     # unpresc 2E31
    #"HLT_IsoMu9",  # unpresc 6E31
    #"HLT_Mu11",    # unpresc 6E31
    #"HLT_IsoMu13", # unpresc 2E32
    "HLT_Mu15_v1",    # unpresc 2E32
]
#HLT2Ms = [
#    "HLT_DoubleMu3", # unpresc any
#    "HLT_DoubleMu5", # unpresc any
#]
HLT_Any1E = "||".join(["!triggerObjectMatchesByPath('%s').empty" % (X,) for X in HLT1Es])
HLT_Any1M = "||".join(["!triggerObjectMatchesByPath('%s').empty" % (X,) for X in HLT1Ms])

tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.5 " +
                      "&& ("+ HLT_Any1E   +")" +
                      "&& ("+ ELE_ID_CUT  +")" +
                      "&& ("+ ELE_ISO_CUT +")")),
)

tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.1 "+
                      "&& (" + HLT_Any1M + ")" +
                      "&& ("+MUON_ID_CUT + ")" +
                      "&& ("+MUON_ISO_CUT+ ")")),
)

tp_lepton_effs_reco = cms.Sequence( 
    tagMuons     +
    tagElectrons
)

nVerticesElectrons = cms.EDProducer("VertexMultiplicityCounter",
    probes  = cms.InputTag("boostedElectrons"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"), 
)
nVerticesMuons = nVerticesElectrons.clone(probes = "boostedMuons")
nJetsElectrons = cms.EDProducer("CandCleanedMultiplicityCounter", # this is a property of the *pair*
    pairs   = cms.InputTag("tpElEl"),
    objects = cms.InputTag("cleanPatJets"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 25"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)
nJetsElectronsNoPU = nJetsElectrons.clone(objects = "cleanPatJetsNoPU")
nJetsMuons     = nJetsElectrons.clone(pairs = "tpMuMu")
nJetsMuonsNoPU = nJetsMuons.clone(objects = "cleanPatJetsNoPU")
varsElectrons = cms.Sequence(nVerticesElectrons + nJetsElectrons + nJetsElectronsNoPU)
varsMuons = cms.Sequence(nVerticesMuons + nJetsMuons + nJetsMuonsNoPU)

tpElEl = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagElectrons@+ boostedElectrons@-'),
    cut = cms.string('70 < mass < 120 && daughter(1).pt > 10'),
)
tpMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('tagMuons@+ boostedMuons@-'),
    cut = cms.string('70 < mass < 120 && daughter(1).pt > 10'),
)

TPCommonVariables = cms.PSet(
    pt    = cms.string("pt"),
    eta   = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    phi    = cms.string("phi"),
)

EleTriggers = cms.PSet(HLT_Any1E = cms.string(HLT_Any1E))
MuTriggers  = cms.PSet(HLT_Any1M = cms.string(HLT_Any1M))
for X in HLT1Es: 
    setattr(EleTriggers, X, cms.string("!triggerObjectMatchesByPath('%s').empty" % (X,)))
for X in HLT1Ms: 
    setattr(MuTriggers,  X, cms.string("!triggerObjectMatchesByPath('%s').empty" % (X,)))
MuIDFlags = cms.PSet(
    passID  = cms.string(MUON_ID_CUT),
    passIso = cms.string(MUON_ISO_CUT),
    passGlb = cms.string("isGlobalMuon"),
    passTM  = cms.string("isTrackerMuon"),
)
EleIDFlags = cms.PSet(
    passID      = cms.string(ELE_ID_CUT),
    passIso     = cms.string(ELE_ISO_CUT),
    passTkBar   = cms.string("dr03TkSumPt/pt < 0.09"),
    passEcalBar = cms.string("dr03EcalRecHitSumEt/pt < 0.07"),
    passHcalBar = cms.string("dr03HcalTowerSumEt/pt < 0.10"),
    passTkEnd   = cms.string("dr03TkSumPt/pt < 0.04"),
    passEcalEnd = cms.string("dr03EcalRecHitSumEt/pt < 0.05"),
    passHcalEnd = cms.string("dr03HcalTowerSumEt/pt < 0.025"),
)

tpTreeElEl = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpElEl"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        TPCommonVariables,
        nVtx     = cms.InputTag("nVerticesElectrons"),
    ),
    flags = cms.PSet(
        EleIDFlags,
        EleTriggers
    ),
    tagVariables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
        EleTriggers,
    ),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsElectrons"),
        nJetNoPU = cms.InputTag("nJetsElectronsNoPU"),
    ),
    pairFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)
tpTreeMuMu = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpMuMu"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        TPCommonVariables,
        nVtx     = cms.InputTag("nVerticesMuons"),
    ),
    flags = cms.PSet(
        MuIDFlags,
        MuTriggers,
    ),
    tagVariables = cms.PSet(
        pt     = cms.string("pt"),
        eta    = cms.string("eta"),
        abseta = cms.string("abs(eta)"),
    ),
    tagFlags = cms.PSet(
        MuIDFlags,
        MuTriggers,
    ),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsMuons"),
        nJetNoPU = cms.InputTag("nJetsMuonsNoPU"),
    ),
    pairFlags = cms.PSet(),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

tnpSimpleSequenceElEl = cms.Sequence(
    tp_lepton_effs_reco *
    tpElEl              *
    varsElectrons       *
    tpTreeElEl
)
tnpSimpleSequenceMuMu = cms.Sequence(
    tp_lepton_effs_reco *
    tpMuMu              *
    varsMuons           *
    tpTreeMuMu
)


