import FWCore.ParameterSet.Config as cms

#from WWAnalysis.AnalysisStep.wwMuons_cfi     import MUON_ISO_CUT, MUON_ID_CUT, MUON_IP_CUT
#from WWAnalysis.AnalysisStep.electronIDs_cff import ELE_ISO_LH_90_2011 as ELE_ISO_CUT
#from WWAnalysis.AnalysisStep.electronIDs_cff import ELE_ID_LH_90_2011  as ELE_ID_CUT
#from WWAnalysis.AnalysisStep.electronIDs_cff import ELE_NOCONV, ELE_IP

from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT_OLD as MUON_ID_2010_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_ISO as MUON_ISO_CUT
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_MERGE_IP as MUON_IP_CUT
#
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID as ELE_ID_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ID2 as ELE_ID2_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_ISO as ELE_ISO_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_IP as ELE_IP_CUT
from WWAnalysis.AnalysisStep.wwElectrons_cfi import ELE_MERGE_CONV as ELE_NOCONV


HLT1Es = [ 'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT' ]
HLT1Ms = [ 'HLT_IsoMu17', 'HLT_Mu24' ]
HLT_Any1E = "||".join(["!triggerObjectMatchesByPath('%s_v*').empty" % (X,) for X in HLT1Es])
HLT_Any1M = "||".join(["!triggerObjectMatchesByPath('%s_v*').empty" % (X,) for X in HLT1Ms])

tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("boostedElectrons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.5 " +
                      "&& ("+ HLT_Any1E   +")" +
                      "&& ("+ ELE_ID_CUT  +")" +
                      "&& ("+ ELE_IP_CUT      +")" +
                      "&& ("+ ELE_NOCONV  +")" +
                      "&& ("+ ELE_ISO_CUT +")")),

)

tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("boostedMuons"),
    cut = cms.string(("pt > 20 && abs(eta) < 2.1 "+
                      "&& (" + HLT_Any1M + ")" +
                      "&& ("+MUON_ID_2010_CUT + ")" +
                      "&& ("+MUON_IP_CUT + ")" +
                      "&& ("+MUON_ISO_CUT+ ")")),

)

tp_lepton_effs_reco = cms.Sequence( 
    tagMuons     +
    tagElectrons
)

nVerticesElectrons = cms.EDProducer("VertexMultiplicityCounter",
    probes  = cms.InputTag("boostedElectrons"),
    objects = cms.InputTag("goodPrimaryVertices"),
    objectSelection = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"), 
)
nVerticesMuons = nVerticesElectrons.clone(probes = "boostedMuons")
nJetsElectrons = cms.EDProducer("CandCleanedMultiplicityCounter", # this is a property of the *pair*
    pairs   = cms.InputTag("tpElEl"),
    objects = cms.InputTag("slimPatJetsTriggerMatch"),
    objectSelection = cms.string("abs(eta) < 5 && pt > 30"), 
    minTagObjDR   = cms.double(0.3),
    minProbeObjDR = cms.double(0.3),
)
nJetsMuons     = nJetsElectrons.clone(pairs = "tpMuMu")
varsElectrons = cms.Sequence(nVerticesElectrons + nJetsElectrons)
varsMuons = cms.Sequence(nVerticesMuons + nJetsMuons)

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

EleTriggers = cms.PSet(HLT_Any1E = cms.string(HLT_Any1E)); EleTriggersMC = cms.PSet()
MuTriggers  = cms.PSet(HLT_Any1M = cms.string(HLT_Any1M)); MuTriggersMC = cms.PSet()

EleTriggers.HLT_Ele17_CaloIdL_CaloIsoVL = cms.string("!triggerObjectMatchesByFilter('hltEle17CaloIdLCaloIsoVLPixelMatchFilter').empty")
EleTriggers.HLT_Ele8_CaloIdL_CaloIsoVL  = cms.string("!triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*').empty")
EleTriggers.HLT_Ele17_CaloIdL_CaloIsoVL_Unseeded = cms.string(
    "!triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*').empty          &&" +
    " triggerObjectMatchesByPath('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*').at(0).pt >= 17 "
)
EleTriggersMC.HLT_Ele17_CaloIdL_CaloIsoVL = cms.string("!triggerObjectMatchesByFilter('hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEPixelMatchFilter').empty")
EleTriggersMC.HLT_Ele8_CaloIdL_CaloIsoVL  = cms.string("!triggerObjectMatchesByPath('HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*').empty")
EleTriggersMC.HLT_Ele17_CaloIdL_CaloIsoVL_Unseeded = cms.string(
    "!triggerObjectMatchesByPath('HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*').empty          &&" +
    " triggerObjectMatchesByPath('HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*').at(0).pt >= 17 "
)


#MuTriggers.HLT_DoubleMu7 = cms.string("!triggerObjectMatchesByPath('HLT_DoubleMu7_v*').empty")
#MuTriggers.HLT_Mu24      = cms.string("!triggerObjectMatchesByPath('HLT_Mu24_v*').empty")
#MuTriggers.HLT_Mu17_EMu = cms.string(
#    "!triggerObjectMatchesByFilter('hltDiMuon3L2PreFiltered0').empty           && " +
#    " triggerObjectMatchesByFilter('hltDiMuon3L2PreFiltered0').at(0).pt >= 12  && "+
#    "!triggerObjectMatchesByPath('HLT_DoubleMu7_v*').empty                     && "+
#    " triggerObjectMatchesByPath('HLT_DoubleMu7_v*').at(0).pt >= 17"
#)
#MuTriggers.HLT_Mu8_EMu = cms.string(
#    "!triggerObjectMatchesByFilter('hltDiMuon3L2PreFiltered0').empty           && " +
#    " triggerObjectMatchesByFilter('hltDiMuon3L2PreFiltered0').at(0).pt >= 5   && "+
#    "!triggerObjectMatchesByPath('HLT_DoubleMu7_v*').empty                     && "+
#    " triggerObjectMatchesByPath('HLT_DoubleMu7_v*').at(0).pt >= 8"
#)
#MuTriggersMC.HLT_DoubleMu7 = cms.string(
#    "!triggerObjectMatchesByPath('HLT_DoubleMu5_v*').empty  && " +
#    " triggerObjectMatchesByPath('HLT_DoubleMu5_v*').at(0).pt >= 7"
#)
#MuTriggersMC.HLT_Mu24      = cms.string(
#    "!triggerObjectMatchesByPath('HLT_Mu21_v*').empty       && "+
#    " triggerObjectMatchesByPath('HLT_Mu21_v*').at(0).pt >= 24"
#)
#MuTriggersMC.HLT_Mu17_EMu = cms.string(
#    "!triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty           && " +
#    " triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').at(0).pt >= 12  && "+
#    "!triggerObjectMatchesByPath('HLT_DoubleMu5_v*').empty                     && "+
#    " triggerObjectMatchesByPath('HLT_DoubleMu5_v*').at(0).pt >= 17"
#)
#MuTriggersMC.HLT_Mu8_EMu = cms.string(
#    "!triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').empty          && " +
#    " triggerObjectMatchesByFilter('hltDiMuonL2PreFiltered0').at(0).pt >= 5  && "+
#    "!triggerObjectMatchesByPath('HLT_DoubleMu5_v*').empty                   && "+
#    " triggerObjectMatchesByPath('HLT_DoubleMu5_v*').at(0).pt >= 8"
#)
 

MuIDFlags = cms.PSet(
    passID  = cms.string(MUON_ID_CUT),
    passID_2010  = cms.string(MUON_ID_2010_CUT),
    passIP  = cms.string(MUON_IP_CUT),
    passIso = cms.string(MUON_ISO_CUT),
    passGlb = cms.string("isGlobalMuon"),
    passTM  = cms.string("isTrackerMuon"),
)
EleIDFlags = cms.PSet(
    passID      = cms.string(ELE_ID_CUT),
    passID2      = cms.string(ELE_ID2_CUT),
    passIso     = cms.string(ELE_ISO_CUT),
    passIP      = cms.string(ELE_IP_CUT),
    passConvR   = cms.string(ELE_NOCONV),
)

EleTriggers = cms.PSet(); 
MuTriggers  = cms.PSet(); 


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
    tagVariables = cms.PSet(),
    tagFlags = cms.PSet(
        EleTriggers
    ),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsElectrons"),
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
    tagVariables = cms.PSet(),
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(
        nJet     = cms.InputTag("nJetsMuons"),
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
