import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/8ED4CCDF-2CC9-DF11-BD8E-001D09F24F65.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/8E9C2BC5-52C9-DF11-B549-0030487C8E00.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/885FCB74-32C9-DF11-83A6-003048F024C2.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/86FF8497-4EC9-DF11-B371-001D09F24493.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/80B39692-34C9-DF11-B640-001D09F29114.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/7EC7B364-29C9-DF11-8F5E-0030487CAF0E.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/7E4BE0FA-64C9-DF11-9748-001D09F23A3E.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/7C969105-38C9-DF11-9C1F-000423D33970.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/7A2FAADD-91C9-DF11-B4AA-003048F118AA.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/76B5143A-85C9-DF11-9BF2-0030487CD17C.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/74D3E3B4-38C9-DF11-80DD-003048F1BF66.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/66FAFE2F-4DC9-DF11-B663-001D09F24D4E.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/643257BB-38C9-DF11-A7FA-003048F118E0.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/64157A10-60C9-DF11-A6AF-003048CFB40C.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/62AE9531-69C9-DF11-87CE-0030487CD700.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/528B390E-3FC9-DF11-B0F6-001617C3B6C6.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/4C60ED02-36C9-DF11-8307-001617C3B69C.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/4C292275-61C9-DF11-9F72-0019B9F72D71.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/4C0D2F46-22C9-DF11-B0C2-0030487CD718.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/480F85BD-10C9-DF11-880C-001D09F2AD84.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/46EFFBFA-64C9-DF11-8AE3-001D09F2932B.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/46681065-1DC9-DF11-B177-0030487A18A4.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/441EA234-3AC9-DF11-B8B1-003048F1C58C.root',
	'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/644/42DB4042-64C9-DF11-B164-0030487C635A.root',
    ),
)
process.source.fileNames = [ 'file:WZEG_CentralSkim.root' ]
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')


from WWAnalysis.Filters.electronSelections_cff import *

#process.idElectrons = cms.EDFilter("GsfElectronSelector",
#    src = cms.InputTag("gsfElectrons"),
#    cut = cms.string(ELE_ID_CUT + " && " +ELE_ISO_CUT)
#)
from WWAnalysis.Filters.convRejectionSelection_cfi import convRejectionSelection
#process.goodElectrons = convRejectionSelection.clone(src = "idElectrons", filter = False)
process.convRejectionElectrons = convRejectionSelection.clone(src = "gsfElectrons", filter = False)

process.expectedHitsComputer = cms.EDProducer("ExpectedHitsComputer",
    #inputColl       = cms.InputTag("goodElectrons"),
    inputColl       = cms.InputTag("gsfElectrons"),
    useGsfTrack     = cms.bool(True),
    objectSelection = cms.string(""),
    propagator         = cms.string('PropagatorWithMaterialOpposite'),
    navigationSchool   = cms.string('SimpleNavigationSchool'),
    measurementTracker = cms.string(''),
)

process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
process.patElectrons.addElectronID = False
process.patElectrons.addGenMatch   = False
process.patElectrons.userData.userInts.src = cms.VInputTag(cms.InputTag("expectedHitsComputer"))
#process.patElectrons.electronSource = "goodElectrons"

process.patElectronsWithConvR = cms.EDProducer("PATElectronCleaner",
    src = cms.InputTag("patElectrons"), 
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
        convR = cms.PSet(
           src       = cms.InputTag("convRejectionElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""), # no preselection 
           deltaR              = cms.double(0.0001), # they're the same!
           checkRecoComponents = cms.bool(False), # doesn't work
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False), # overlaps don't cause the electron to be discared
        )
    ),
    finalCut = cms.string(''),
)

## Trigger matching (here comes the evil)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patElectronsWithConvR" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL1IsoRecoEcalCandidate', 'hltL1NonIsoRecoEcalCandidate' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)
process.eleIdTriggerMatchHLT = process.eleTriggerMatchHLT.clone(collectionTags = cms.vstring('hltPixelMatchElectronsL1Iso', 'hltPixelMatchElectronsL1NonIso') )
process.patElectronsWithTrigger = cms.EDProducer( "PATTriggerMatchElectronEmbedder",
    src     = cms.InputTag(  "patElectronsWithConvR" ),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'), cms.InputTag('eleIdTriggerMatchHLT'))
)
IDTriggerPSet = cms.PSet(
    passID  = cms.string(ELE_ID_CUT),
    passIso = cms.string(ELE_ISO_CUT),
    passConvR = cms.string("hasOverlaps('convR')"),
    hasHLTSC = cms.string("!triggerObjectMatchesByCollection('hltL1IsoRecoEcalCandidate').empty || !triggerObjectMatchesByCollection('hltL1NonIsoRecoEcalCandidate').empty"),
    hasHLTID = cms.string("!triggerObjectMatchesByCollection('hltPixelMatchElectronsL1Iso').empty || !triggerObjectMatchesByCollection('hltPixelMatchElectronsL1NonIso').empty"),
)
HLTs = [
    "HLT_Ele10_SW_EleId_L1R",             # MC
    "HLT_Ele15_SW_L1R",                   # MC
    "HLT_Ele12_SW_EleIdIsol_L1R",         # 2E31
    "HLT_Ele17_SW_LooseEleId_L1R",        # 2E31
    "HLT_Ele17_SW_CaloEleId_L1R",         # 2E31
    "HLT_DoubleEle10_SW_L1R",             # 2E31, MC
    "HLT_Ele12_SW_TighterEleIdIsol_L1R",  # 6E31
    "HLT_Ele17_SW_TightEleId_L1R",        # 6E31
    "HLT_DoubleEle15_SW_L1R",             # 6E31
]
for X in HLTs:
    setattr(IDTriggerPSet, X, cms.string("!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" %(X,X)))
IDTriggerPSet.HLT_Any = cms.string("||".join(["!triggerObjectMatchesByPath('%s').empty || !triggerObjectMatchesByPath('%s_v1').empty" % (X,X) for X in HLTs]))

process.tagElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string("pt > 20 && %s && %s" % (ELE_ID_CUT, ELE_ISO_CUT)),
)

process.probeElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string("!(1.4442 < abs(eta) < 1.5666)"), 
)

process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 140'),
    decay = cms.string('tagElectrons@+ probeElectrons@-')
)

process.tagElectronsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    src = cms.InputTag("tagElectrons"),
    matched = cms.InputTag("genParticles"),
    pdgId = cms.vint32(11),
    distMin = cms.double(0.3),
)
process.probeElectronsMCMatch = process.tagElectronsMCMatch.clone(src = "probeElectrons")

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
        fbrem = cms.string("fbrem"),
        nEihOld = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfLostHits"),
        nEihNew = cms.string("userInt('expectedHitsComputer')"),
    ),
    flags = cms.PSet(
        IDTriggerPSet,
        passEihOld = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfLostHits == 0"),
        passEihNew = cms.string("userInt('expectedHitsComputer') == 0"),
    ),
    tagVariables =  cms.PSet(
        pt    = cms.string("pt"),
        eta   = cms.string("eta"),
    ),
    tagFlags = cms.PSet(IDTriggerPSet),
    isMC = cms.bool(False),
    tagMatches   = cms.InputTag("tagElectronsMCMatch"),
    probeMatches = cms.InputTag("probeElectronsMCMatch"),
    motherPdgId = cms.vint32(22, 23),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff  = cms.bool(True),
    allProbes = cms.InputTag("probeElectrons"),
)
process.tnpSimpleSequence = cms.Sequence(
    process.tagElectrons   * #process.tagElectronsMCMatch   +
    process.probeElectrons * #process.probeElectronsMCMatch +
    process.tpPairs    +
    process.tpTree
)

process.tagAndProbe = cms.Path( 
    process.convRejectionElectrons *
    process.expectedHitsComputer *
    process.patElectrons *
    process.patElectronsWithConvR *
    (process.patTrigger * (process.eleTriggerMatchHLT + process.eleIdTriggerMatchHLT)) *
    process.patElectronsWithTrigger *
    process.tnpSimpleSequence
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpZee_Data.root"))
