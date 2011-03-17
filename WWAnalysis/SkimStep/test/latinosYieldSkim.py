import FWCore.ParameterSet.Config as cms

process = cms.Process("Yield")

isMC = False
# isMC = RMMEMC

#  _____               _____                               _                
# |  __ \             |  __ \                             | |               
# | |__) |   _ _ __   | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
# |  _  / | | | '_ \  |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | \ \ |_| | | | | | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# |_|  \_\__,_|_| |_| |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
# 

from PhysicsTools.PatAlgos.tools.coreTools import *
 
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# process.GlobalTag.globaltag = 'RMMEGlobalTag'
process.GlobalTag.globaltag = 'START39_V8::All'


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'RMMEFN'
        'file:/nfs/bluearc/group/edm/hww/Winter10.Flat/hww.flat.root',
#         'file:/nfs/bluearc/group/edm/hww/Winter10.Flat/ttbar.flat.root',
    )
)

#from glob import glob
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/oct29Skim/WWFull/*/*.root')]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('latinosYieldSkim.root'),
#     fileName = cms.untracked.string('RMMEFN'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
        # Leptons
        'keep *_patMuonsWithTrigger_*_*',
        'keep *_patElectronsWithTrigger_*_*',
        'keep *_selectedPatTaus_*_*',
        # Jets
        'keep patJets_cleanPatJets_*_*',
        'keep patJets_cleanPatJetsPF_*_*',
        'keep patJets_cleanPatJetsNoPU_*_*',
        'keep patJets_cleanPatJetsJPT_*_*',
        # Tracking
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlinePrimaryVerticesWithBS_*_*',
        'keep *_offlineBeamSpot_*_*',
        # MET
        'keep *_tcMet_*_*',
        'keep *_met_*_*',
        'keep *_pfMet_*_*',
        # MC
        'keep *_prunedGen_*_*',
        'keep *_genMetTrue_*_*',
        # Trigger
        'keep *_TriggerResults_*_*',
        'keep *_vertexMapProd_*_*',
        # Misc
        'keep *_addPileupInfo_*_*',
        'keep *_chargedMetProducer_*_*',
        'keep *_mergedSuperClusters_*_'+process.name_(),
        'keep *_kt6PFJets_rho_'+process.name_(),
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'wwAllPath',)),
)
process.outpath = cms.EndPath(process.out)



#  _______   _                         ______ _ _ _            
# |__   __| (_)                       |  ____(_) | |           
#    | |_ __ _  __ _  __ _  ___ _ __  | |__   _| | |_ ___ _ __ 
#    | | '__| |/ _` |/ _` |/ _ \ '__| |  __| | | | __/ _ \ '__|
#    | | |  | | (_| | (_| |  __/ |    | |    | | | ||  __/ |   
#    |_|_|  |_|\__, |\__, |\___|_|    |_|    |_|_|\__\___|_|   
#               __/ | __/ |                                    
#              |___/ |___/                                  
# 

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

#Trigger Matching (hopefully?!?)
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.onlyStandAlone = True
process.patTrigger.processName  = '*' 
process.patTriggerEvent.processName = '*'

#
#  _____            _                _              
# |  __ \          | |              | |             
# | |__) | __ ___  | |     ___ _ __ | |_ ___  _ __  
# |  ___/ '__/ _ \ | |    / _ \ '_ \| __/ _ \| '_ \ 
# | |   | | |  __/ | |___|  __/ |_) | || (_) | | | |
# |_|   |_|  \___| |______\___| .__/ \__\___/|_| |_|
#                             | |                   
#                             |_|                   
#   _____                                      
#  / ____|                                     
# | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
#  \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
#  ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                 | |                          
#                 |_|                          
# 
  
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)

process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("WWAnalysis.SkimStep.cutsInCategoriesHWWElectronIdentificationV04_cfi")

process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.egammaIDLikelihood = process.eidLikelihoodExt.clone()

process.eIdSequence = cms.Sequence(
    process.eidVeryLoose +
    process.eidLoose +
    process.eidMedium +
    process.eidTight +
    process.eidSuperTight +
    process.eidHyperTight1 +
    process.eidHWWVeryLoose +
    process.eidHWWLoose +
    process.eidHWWMedium +
    process.eidHWWTight +
    process.eidHWWSuperTight +
    process.eidHWWHyperTight1 +
    process.eidHWWHyperTight2 +
    process.eidHWWHyperTight3 +
    process.egammaIDLikelihood 
)

process.load("WWAnalysis.Tools.betaValueMapProducer_cfi")
process.load("WWAnalysis.Tools.rhoValueMapProducer_cfi")

process.betaMu = process.betaValueMapProducer.clone()
process.betaEl = process.betaValueMapProducer.clone()
process.betaEl.leptonTag = "gsfElectrons"
process.betaEl.dEtaVeto = 0.015
process.betaEl.dRVeto = 0.0

process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJets","rho",process.name_()))
process.rhoEl = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJets","rho",process.name_()))

process.valueMaps = cms.Sequence(
    process.betaMu +
    process.betaEl +
    process.rhoMu +
    process.rhoEl 
)
    
process.preLeptonSequence = cms.Sequence(
    process.kt6PFJets * process.valueMaps +
    process.offlinePrimaryVertices +
    process.eIdSequence 
)

#  _____               _____ _    _           
# / ____|             / ____| |  (_)          
#| |  __  ___ _ __   | (___ | | ___ _ __ ___  
#| | |_ |/ _ \ '_ \   \___ \| |/ / | '_ ` _ \ 
#| |__| |  __/ | | |  ____) |   <| | | | | | |
# \_____|\___|_| |_| |_____/|_|\_\_|_| |_| |_|
#                                             

process.prunedGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop  *  ",
        "keep++ pdgId =   {Z0}",
        "++keep pdgId =   {Z0}",
        "keep++ pdgId =   {W+}",
        "++keep pdgId =   {W+}",
        "keep++ pdgId =   {W-}",
        "++keep pdgId =   {W-}",
        "keep++ pdgId =   {h0}",
        "++keep pdgId =   {h0}",
        "keep++ pdgId =   {e+}",
        "++keep pdgId =   {e+}",
        "keep++ pdgId =   {e-}",
        "++keep pdgId =   {e-}",
        "keep++ pdgId =  {mu+}",
        "++keep pdgId =  {mu+}",
        "keep++ pdgId =  {mu-}",
        "++keep pdgId =  {mu-}",
        "++keep pdgId =      5",
        "++keep pdgId =     -5",
        "++keep pdgId =      4",
        "++keep pdgId =     -4",
        "++keep pdgId =     12",
        "++keep pdgId =     14",
        "++keep pdgId =     16",
        "++keep pdgId =    -12",
        "++keep pdgId =    -14",
        "++keep pdgId =    -16",
        "++keep pdgId = {tau+}",
        "++keep pdgId = {tau-}",
    )
)

#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 

process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")

#Code for the liklihood eID
# process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
# process.liklihoodID = process.eidLikelihoodExt.clone() 
# process.liklihoodID.reducedEndcapRecHitCollection = "reducedEcalRecHitsEB"
# process.liklihoodID.reducedBarrelRecHitCollection = "reducedEcalRecHitsEE"

process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patElectrons" ),
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
    src     = cms.InputTag(  "patElectrons" ),
    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'), cms.InputTag('eleIdTriggerMatchHLT'))
) 

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = False
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = True
process.electronMatch.matched = "prunedGen"
process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dcot"),
    cms.InputTag("betaEl"),
    cms.InputTag("rhoEl"),
#     cms.InputTag("liklihoodID")
#     cms.InputTag("kt6PFJets","rho",process.name_())
)

#Set the Pat Electrons to use the eID
from PhysicsTools.PatAlgos.tools.helpers import listModules
for module in listModules(process.eIdSequence):
    setattr(process.patElectrons.electronIDSources,module.label(),cms.InputTag(module.label()))

process.load("WWAnalysis.Tools.convValueMapProd_cfi")

if isMC: 
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd + 
#         process.liklihoodID +
        process.patTrigger + 
        process.prunedGen * 
        process.electronMatch 
    )
else:
    removeMCMatching(process, ['Electrons'])
    process.preElectronSequence = cms.Sequence(
        process.convValueMapProd + 
#         process.liklihoodID +
        process.patTrigger 
    )


process.elSeq = cms.Sequence( 
    process.preElectronSequence * 
    process.patElectrons *
    (process.eleTriggerMatchHLT + process.eleIdTriggerMatchHLT) *
    process.patElectronsWithTrigger 
)

#  __  __                     _____      _   _
# |  \/  |                   |  __ \    | | | |
# | \  / |_   _  ___  _ __   | |__) |_ _| |_| |__
# | |\/| | | | |/ _ \| '_ \  |  ___/ _` | __| '_ \
# | |  | | |_| | (_) | | | | | |  | (_| | |_| | | |
# |_|  |_|\__,_|\___/|_| |_| |_|   \__,_|\__|_| |_|
#

process.muonTriggerMatchHLT = cms.EDProducer( 'PATTriggerMatcherDRDPtLessByR',
    src     = cms.InputTag( 'patMuons' ),
    matched = cms.InputTag( 'patTrigger' ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltL3MuonCandidates' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.1 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)
process.patMuonsWithTrigger = cms.EDProducer( 'PATTriggerMatchMuonEmbedder',
    src     = cms.InputTag(  'patMuons' ),
    matches = cms.VInputTag('muonTriggerMatchHLT')
)

process.patMuons.embedPFCandidate = False
process.patMuons.embedTrack = True
process.patMuons.userData.userFloats.src = cms.VInputTag(
#     cms.InputTag("kt6PFJets","rho",process.name_())
    cms.InputTag("betaMu"),
    cms.InputTag("rhoMu"),
)

process.muonMatch.matched = "prunedGen"

if isMC: 
    if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
        process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
        from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
        addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
        process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger + process.muonClassificationByHits)
        process.MessageLogger.suppressWarning += ['classByHitsGlb'] # kill stupid RPC hit associator warning
    else:
        process.preMuonSequence = cms.Sequence(process.prunedGen * process.muonMatch + process.patTrigger )
else:
    removeMCMatching(process, ['Muons'])
    process.preMuonSequence = cms.Sequence(process.patTrigger)


process.muSeq = cms.Sequence( 
    process.preMuonSequence *
    process.patMuons *
    process.muonTriggerMatchHLT *
    process.patMuonsWithTrigger
)


#   _____ _    _             _____      _   _     
#  / ____| |  (_)           |  __ \    | | | |    
# | (___ | | ___ _ __ ___   | |__) |_ _| |_| |__  
#  \___ \| |/ / | '_ ` _ \  |  ___/ _` | __| '_ \ 
#  ____) |   <| | | | | | | | |  | (_| | |_| | | |
# |_____/|_|\_\_|_| |_| |_| |_|   \__,_|\__|_| |_|
#                                                 

process.allLeps = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("patMuonsWithTrigger"), cms.InputTag("patElectronsWithTrigger"))
)

process.allDiLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('allLeps allLeps'),
    cut = cms.string(
        'deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.05 && ' + 
        'min(daughter(0).pt,daughter(1).pt) >  8 && ' +
        'max(daughter(0).pt,daughter(1).pt) > 18'
    ),
    checkCharge = cms.bool(False)
)

process.countWWAll  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("allDiLep"), minNumber = cms.uint32(1))
process.wwAllSeq = cms.Sequence( process.allLeps * process.allDiLep * process.countWWAll )

#               _               _____      _ _           _   _                 
#    /\        | |             / ____|    | | |         | | (_)                
#   /  \  _   _| |_ _ __ ___  | |     ___ | | | ___  ___| |_ _  ___  _ __  ___ 
#  / /\ \| | | | __| '__/ _ \ | |    / _ \| | |/ _ \/ __| __| |/ _ \| '_ \/ __|
# / ____ \ |_| | |_| | |  __/ | |___| (_) | | |  __/ (__| |_| | (_) | | | \__ \
#/_/    \_\__,_|\__|_|  \___|  \_____\___/|_|_|\___|\___|\__|_|\___/|_| |_|___/
#                                                                              

from PhysicsTools.PatAlgos.tools.jetTools import *

process.load("PhysicsTools.PFCandProducer.pfNoPileUp_cff")  
process.pfPileUp.Enable = True    
process.ak5PFJetsNoPU = process.ak5PFJets.clone(    
    src =   cms.InputTag("pfNoPileUp")  
)
process.jetSequence = cms.Sequence(   
    process.pfNoPileUpSequence  +   
    process.ak5PFJetsNoPU   
)  

if not isMC:
    removeMCMatching(process, ['Jets'])

addJetCollection(
    process,
    cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
    algoLabel    = "JPT",
    typeLabel    = "",
    doJTA        = False,
    doBTagging   = True,
    #jetCorrLabel = ('AK5', 'PF'),
    jetCorrLabel = ('AK5JPT',cms.vstring(['L1JPTOffset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

addJetCollection(
    process,
    cms.InputTag("ak5PFJetsNoPU"),
    algoLabel    = "NoPU",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5', 'PF'),
    jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

addJetCollection(
    process,
    cms.InputTag("ak5PFJets"),
    algoLabel    = "PF",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5', 'PF'),
    jetCorrLabel = ('AK5PF',cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

if not isMC:
    process.patJets.addGenJetMatch = False
    process.patJetsPF.addGenJetMatch = False
    process.patJetsPF.addGenPartonMatch = False
    process.patJetsPF.getJetMCFlavour = False
    process.patJetsNoPU.addGenJetMatch = False
    process.patJetsNoPU.addGenPartonMatch = False
    process.patJetsNoPU.getJetMCFlavour = False
    process.patJetsJPT.addGenJetMatch = False
    process.patJetsJPT.addGenPartonMatch = False
    process.patJetsJPT.getJetMCFlavour = False

process.patJets.addTagInfos = False
process.patJetsPF.addTagInfos = False
process.patJetsNoPU.addTagInfos = False
process.patJetsJPT.addTagInfos = False
process.patJets.embedPFCandidates = False
process.patJetsPF.embedPFCandidates = False
process.patJetsNoPU.embedPFCandidates = False
process.patJetsJPT.embedPFCandidates = False
process.patJets.addAssociatedTracks = False
process.patJetsPF.addAssociatedTracks = False
process.patJetsNoPU.addAssociatedTracks = False
process.patJetsJPT.addAssociatedTracks = False

process.patJetCorrFactors.levels = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])

process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("patJets"),
    preselection = cms.string('pt>10'),
    checkOverlaps = cms.PSet(
        ele = cms.PSet(
           src       = cms.InputTag("patElectronsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False),
        ),
        mu = cms.PSet(
           src       = cms.InputTag("patMuonsWithTrigger"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False),
           pairCut             = cms.string(""),
           requireNoOverlaps = cms.bool(False),
        ),
    ),
    finalCut = cms.string(''),
)
process.cleanPatJetsPF = process.cleanPatJets.clone( src = cms.InputTag("patJetsPF") )
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )
process.cleanPatJetsJPT = process.cleanPatJets.clone( src = cms.InputTag("patJetsJPT") )

process.load("WWAnalysis.Tools.vertexSumPtMapProd_cfi")

process.mergedSuperClusters = cms.EDProducer("SuperClusterCombiner",
    labels = cms.VInputTag(
        cms.InputTag("correctedHybridSuperClusters"),
        cms.InputTag("correctedMulti5x5SuperClustersWithPreshower")
    )
)


process.autreSeq = cms.Sequence( 
    process.jetSequence * 
    process.makePatJets * ( 
        process.cleanPatJets + 
        process.cleanPatJetsPF + 
        process.cleanPatJetsNoPU + 
        process.cleanPatJetsJPT 
    ) + 
    process.vertexMapProd + 
    process.mergedSuperClusters
)

#  _____ _                              _   __  __ ______ _______ 
# / ____| |                            | | |  \/  |  ____|__   __|
#| |    | |__   __ _ _ __ __ _  ___  __| | | \  / | |__     | |   
#| |    | '_ \ / _` | '__/ _` |/ _ \/ _` | | |\/| |  __|    | |   
#| |____| | | | (_| | | | (_| |  __/ (_| | | |  | | |____   | |   
# \_____|_| |_|\__,_|_|  \__, |\___|\__,_| |_|  |_|______|  |_|   
#                         __/ |                                   
#                        |___/                                    

process.load("WWAnalysis.Tools.interestingVertexRefProducer_cfi")
process.load("WWAnalysis.Tools.chargedMetProducer_cfi")

process.patMuonsWithTriggerNoSA = cms.EDFilter("PATMuonRefSelector",
    cut = cms.string("type!=8"),
    src = cms.InputTag("patMuonsWithTrigger"),
    filter = cms.bool(False)
)

process.lepsForMET = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("patMuonsWithTriggerNoSA"), cms.InputTag("patElectronsWithTrigger"))
)


process.lowPtLeps = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("lepsForMET"),
    cut = cms.string("pt>8"),
    filter = cms.bool(False),
)

process.interestingVertexRefProducer.leptonTags = [cms.InputTag("lowPtLeps")]

process.chargedParticleFlow = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("charge != 0"),
    filter = cms.bool(False),
)

process.chargedMetProducer.collectionTag = "chargedParticleFlow"
process.chargedMetProducer.vertexTag = "interestingVertexRefProducer"

process.chargedMetSeq = cms.Sequence( ( (
            process.patMuonsWithTriggerNoSA *
            process.lepsForMET * 
            process.lowPtLeps *
            process.interestingVertexRefProducer 
        ) + 
        process.chargedParticleFlow 
    ) * 
    process.chargedMetProducer
)


# _______              
#|__   __|             
#   | | __ _ _   _ ___ 
#   | |/ _` | | | / __|
#   | | (_| | |_| \__ \
#   |_|\__,_|\__,_|___/
#                      

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(
   process,
   pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
   pfTauLabelNew = cms.InputTag('hpsPFTauProducer'),
   postfix=""
)

process.selectedPatTaus.cut = (
   "pt > 15 && " +
   "tauID('leadingTrackFinding') > 0.5 && tauID('byLooseIsolation') > 0.5"
)
if not isMC:
    removeMCMatching(process, ['Taus'])

process.tauSeq = cms.Sequence(process.makePatTaus * process.selectedPatTaus)


#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            


process.wwAllPath = cms.Path( 
    process.preLeptonSequence * ( 
        process.elSeq + 
        process.muSeq + 
        process.tauSeq ) *  
    process.wwAllSeq + 
    process.autreSeq +
    process.chargedMetSeq 
)
process.scrap = cms.Path( 
    process.noscraping 
)

