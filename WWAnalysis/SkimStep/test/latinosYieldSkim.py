import FWCore.ParameterSet.Config as cms

# CMSSW Regular Stuff
process = cms.Process("Yield")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

#Options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

#   _____  __  __ __  __ ______     
#  |  __ \|  \/  |  \/  |  ____|    
#  | |__) | \  / | \  / | |__   ___ 
#  |  _  /| |\/| | |\/| |  __| / __|
#  | | \ \| |  | | |  | | |____\__ \
#  |_|  \_\_|  |_|_|  |_|______|___/
#                                   

isMC = RMMEMC
# isMC = True

# doPF2PATAlso = RMMEPF2PAT
doPF2PATAlso = False

is41XRelease = RMME41X
# is41XRelease = True

process.GlobalTag.globaltag = 'RMMEGlobalTag'
# process.GlobalTag.globaltag = 'START41_V0::All' #'START311_V2::All' #'GR_R_311_V2::All'

# doFakeRates = RMMEFAKE # 'only', 'also' or None
doFakeRates = None  
doBorisGenFilter = False

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
# process.source.fileNames = ['file:/home/mwlebour/data/hww/Hww2l2nu.Spring11.root']

process.out = cms.OutputModule("PoolOutputModule", outputCommands =  cms.untracked.vstring(), fileName = cms.untracked.string('RMMEFN') )
# process.out = cms.OutputModule("PoolOutputModule", outputCommands =  cms.untracked.vstring(), fileName = cms.untracked.string('latinosYieldSkim.root') )


# Gives us preFakeFilter and preYieldFilter
from WWAnalysis.SkimStep.preFilterPaths_cff import addPreFilters
addPreFilters(process,doBorisGenFilter and isMC)

# pat sequence
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *

# Trigger filter
from WWAnalysis.SkimStep.triggerInformation_cff import addTriggerPaths
jetTrigMatches = addTriggerPaths(process)

process.preLeptonSequence = cms.Sequence()

if is41XRelease:
    process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
    process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
    process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
    process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
    process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)
    process.preLeptonSequence += process.offlinePrimaryVertices


# Rho calculations
from WWAnalysis.SkimStep.rhoCalculations_cff import addRhoVariables
addRhoVariables(process,process.preLeptonSequence)
# added them above, so can remove them here
process.patDefaultSequence.remove( process.pfPileUp )
process.patDefaultSequence.remove( process.pfNoPileUp )

# electron ID stuff
from WWAnalysis.SkimStep.electronIDs_cff import addElectronIDs
eidModules = addElectronIDs(process,process.preLeptonSequence)

# generator stuff
from WWAnalysis.SkimStep.generatorInformation_cff import addGeneratorInfo
if isMC: addGeneratorInfo(process,process.preLeptonSequence)



#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 

process.patElectrons.embedPFCandidate = False
process.patElectrons.embedSuperCluster = True
process.patElectrons.embedTrack = True
process.patElectrons.addElectronID = True
process.electronMatch.matched = "prunedGen"
process.patElectrons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","dist"),
    cms.InputTag("convValueMapProd","passVtxConvert"),
    cms.InputTag("betaEl"),
    cms.InputTag("rhoEl"),
    cms.InputTag("rhoElNoPU"),
)
process.patElectrons.isolationValues = cms.PSet(
#     pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutralIso"),
#     pfChargedHadrons = cms.InputTag("isoValElectronWithChargedIso"),
#     pfPhotons = cms.InputTag("isoValElectronWithPhotonIso")
)


#Set the Pat Electrons to use the eID
for module in eidModules:
    setattr(process.patElectrons.electronIDSources,module.label(),cms.InputTag(module.label()))

process.load("WWAnalysis.Tools.convValueMapProd_cfi")
if not is41XRelease:  process.convValueMapProd.conversionLabel = "allConversions"
process.preElectronSequence = cms.Sequence(process.convValueMapProd)




#  __  __                     _____      _   _
# |  \/  |                   |  __ \    | | | |
# | \  / |_   _  ___  _ __   | |__) |_ _| |_| |__
# | |\/| | | | |/ _ \| '_ \  |  ___/ _` | __| '_ \
# | |  | | |_| | (_) | | | | | |  | (_| | |_| | | |
# |_|  |_|\__,_|\___/|_| |_| |_|   \__,_|\__|_| |_|
#


process.patMuons.embedPFCandidate = False
process.patMuons.embedTrack = True
process.patMuons.userData.userFloats.src = cms.VInputTag(
    cms.InputTag("betaMu"),
    cms.InputTag("rhoMu"),
    cms.InputTag("rhoMuNoPU"),
)
process.patMuons.isolationValues = cms.PSet()
process.muonMatch.matched = "prunedGen"
process.preMuonSequence = cms.Sequence()

# Not implemented yet in 41X:
# if isMC: 
#     if False: ## Turn this on to get extra info on muon MC origin, on GEN-SIM-RECO
#         process.load("MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi")
#         from MuonAnalysis.MuonAssociators.muonClassificationByHits_cfi import addUserData as addClassByHits
#         addClassByHits(process.patMuons, labels=['classByHitsGlb'], extraInfo = True)
#         process.muonClassificationByHits = cms.Sequence(process.mix * process.trackingParticlesNoSimHits * process.classByHitsGlb )
#         process.preMuonSequence += process.muonClassificationByHits
#         process.MessageLogger.suppressWarning += ['classByHitsGlb'] # kill stupid RPC hit associator warning

#  _   _                             _     _   _____  ______ ___  _____     _______ 
# | \ | |                   /\      | |   | | |  __ \|  ____|__ \|  __ \ /\|__   __|
# |  \| | _____      __    /  \   __| | __| | | |__) | |__     ) | |__) /  \  | |   
# | . ` |/ _ \ \ /\ / /   / /\ \ / _` |/ _` | |  ___/|  __|   / /|  ___/ /\ \ | |   
# | |\  | (_) \ V  V /   / ____ \ (_| | (_| | | |    | |     / /_| |  / ____ \| |   
# |_| \_|\___/ \_/\_/   /_/    \_\__,_|\__,_| |_|    |_|    |____|_| /_/    \_\_|   
                                                                                     
def addFastJetCorrection(process,label,seq="patDefaultSequence",thisRho="kt6PFJets"):
    corrFact = getattr(process,"patJetCorrFactors"+label)
    setattr(process,"patJetCorrFactorsFastJet"+label,corrFact.clone())
    getattr(process,"patJetCorrFactorsFastJet"+label).levels[0] = 'L1FastJet'
    getattr(process,"patJetCorrFactorsFastJet"+label).rho = cms.InputTag(thisRho,"rho")
    getattr(process,seq).replace(
        getattr(process,"patJetCorrFactors"+label),
        getattr(process,"patJetCorrFactors"+label) +
        getattr(process,"patJetCorrFactorsFastJet"+label) 
    )
    getattr(process,"patJets"+label).jetCorrFactorsSource = cms.VInputTag(
        cms.InputTag("patJetCorrFactorsFastJet"+label) ,
        cms.InputTag("patJetCorrFactors"+label) 
    )

from PhysicsTools.PatAlgos.tools.pfTools import *

if doPF2PATAlso:
    usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=isMC, postfix="PFlow") 
    process.pfNoTauPFlow.enable = False

    if not isMC:
        removeMCMatchingPF2PAT( process, '' )

    # For some reason, with the other functions that I have called, this still needs to be setup:
    process.patPF2PATSequencePFlow.replace(
        process.selectedPatCandidateSummaryPFlow,
        process.selectedPatCandidateSummaryPFlow +
        process.cleanPatMuonsPFlow + 
        process.cleanPatElectronsPFlow + 
        process.cleanPatTausPFlow +
        process.cleanPatJetsPFlow 
    )
    delattr(process.cleanPatJetsPFlow.checkOverlaps,"photons")
    process.patPF2PATSequencePFlow.remove(process.cleanPatPhotonsTriggerMatchPFlow)
    process.patPF2PATSequencePFlow.remove(process.cleanPhotonTriggerMatchHLTPhoton26IsoVLPhoton18PFlow)
    process.patJetsPFlow.embedCaloTowers = False
    process.patJetsPFlow.addTagInfos = False
    process.patJetsPFlow.embedPFCandidates = False
    process.patJetsPFlow.addAssociatedTracks = False
    #Tell PF2PAT to recluster w/ Area calculation on:
    process.pfJetsPFlow.doAreaFastjet = True
    process.pfJetsPFlow.Rho_EtaMax = cms.double(4.5)
    # Turn on secondary JEC w/ FastJet
    addFastJetCorrection(process,"PFlow","patPF2PATSequencePFlow","kt6PFJetsNoPU")

else:
    if not isMC:
        removeMCMatching(process)




#      _      _      _____                                      
#     | |    | |    / ____|                                     
#     | | ___| |_  | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
# _   | |/ _ \ __|  \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
#| |__| |  __/ |_   ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# \____/ \___|\__| |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                                  | |                          
#                                  |_|                          

#Add L2L3Residual if on data:
if isMC:
    myCorrLabels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute')
else:
    if is41XRelease:
        myCorrLabels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual')
    else:
        myCorrLabels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute')

#all the other jets:
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
)

addJetCollection(
    process,
    cms.InputTag("ak5PFJetsNoPU"),
    algoLabel    = "NoPU",
    typeLabel    = "",
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF',myCorrLabels),
    doL1Cleaning = False,
    doL1Counters = True,                 
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True,
    jetIdLabel   = 'ak5',
)

# Some stuff to save space
process.patJets.embedCaloTowers = False
process.patJetsNoPU.embedCaloTowers = False
process.patJets.addTagInfos = False
process.patJetsNoPU.addTagInfos = False
process.patJets.embedPFCandidates = False
process.patJetsNoPU.embedPFCandidates = False
process.patJets.addAssociatedTracks = False
process.patJetsNoPU.addAssociatedTracks = False

# Not set up correctly by PAT:
process.cleanPatJetsNoPU = process.cleanPatJets.clone( src = cms.InputTag("patJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanPatJets,
    process.cleanPatJets +
    process.cleanPatJetsNoPU 
)

for X in jetTrigMatches:
    oldmatch = getattr(process,X)
    newmatch = oldmatch.clone( src = cms.InputTag("cleanPatJetsNoPU") )
    setattr(process, X+'NoPU', newmatch)
    process.patDefaultSequence.replace(oldmatch, oldmatch+newmatch)

process.cleanPatJetsTriggerMatchNoPU = process.cleanPatJetsTriggerMatch.clone( matches = [X+"NoPU" for X in jetTrigMatches], src = cms.InputTag("cleanPatJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanPatJetsTriggerMatch,
    process.cleanPatJetsTriggerMatch +
    process.cleanPatJetsTriggerMatchNoPU
)

process.slimPatJetsTriggerMatch = cms.EDProducer("PATJetSlimmer",
    src = cms.InputTag("cleanPatJetsTriggerMatch"),
    clearJetVars = cms.bool(True),
    clearDaughters = cms.bool(True),
    dropSpecific = cms.bool(False),
)
process.slimPatJetsTriggerMatchNoPU = process.slimPatJetsTriggerMatch.clone( src = "cleanPatJetsTriggerMatchNoPU" ) 
process.patDefaultSequence += (
    process.slimPatJetsTriggerMatch     +
    process.slimPatJetsTriggerMatchNoPU
)

# Other stuff to do for fun:
if doPF2PATAlso:
    process.slimPatJetsTriggerMatchPFlow = process.slimPatJetsTriggerMatch.clone( src = "cleanPatJetsTriggerMatchPFlow" )
    process.patPF2PATSequencePFlow += process.slimPatJetsTriggerMatchPFlow

# Add the fast jet correction:
addFastJetCorrection(process,"")
addFastJetCorrection(process,"NoPU","patDefaultSequence","kt6PFJetsNoPU")

#               _               _____      _ _           _   _                 
#    /\        | |             / ____|    | | |         | | (_)                
#   /  \  _   _| |_ _ __ ___  | |     ___ | | | ___  ___| |_ _  ___  _ __  ___ 
#  / /\ \| | | | __| '__/ _ \ | |    / _ \| | |/ _ \/ __| __| |/ _ \| '_ \/ __|
# / ____ \ |_| | |_| | |  __/ | |___| (_) | | |  __/ (__| |_| | (_) | | | \__ \
#/_/    \_\__,_|\__|_|  \___|  \_____\___/|_|_|\___|\___|\__|_|\___/|_| |_|___/
#                                                                              

process.load("WWAnalysis.Tools.vertexSumPtMapProd_cfi")

# process.mergedSuperClusters = cms.EDProducer("SuperClusterCombiner",
#     labels = cms.VInputTag(
#         cms.InputTag("correctedHybridSuperClusters"),
#         cms.InputTag("correctedMulti5x5SuperClustersWithPreshower")
#     )
# )

process.autreSeq = cms.Sequence( 
    process.vertexMapProd 
#     process.mergedSuperClusters
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
    src = cms.InputTag("boostedMuons"),
    filter = cms.bool(False)
)

process.lepsForMET = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("patMuonsWithTriggerNoSA"), cms.InputTag("boostedElectrons"))
)


process.lowPtLeps = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("lepsForMET"),
    cut = cms.string("pt>8"),
    filter = cms.bool(False),
)

process.interestingVertexRefProducer.leptonTags = [cms.InputTag("lowPtLeps")]

process.chargedMetProducer.collectionTag = "particleFlow"
process.chargedMetProducer.vertexTag = "interestingVertexRefProducer"
process.trackMetProducer = process.chargedMetProducer.clone(minNeutralPt = 99999., maxNeutralEta = 0)

process.chargedMetSeq = cms.Sequence( ( 
        process.patMuonsWithTriggerNoSA *
        process.lepsForMET * 
        process.lowPtLeps *
        process.interestingVertexRefProducer ) * 
    process.chargedMetProducer +
    process.trackMetProducer 
)


# _______              
#|__   __|             
#   | | __ _ _   _ ___ 
#   | |/ _` | | | / __|
#   | | (_| | |_| \__ \
#   |_|\__,_|\__,_|___/
#                      

switchToPFTauHPS(
   process,
   pfTauLabelOld = cms.InputTag('shrinkingConePFTauProducer'),
   pfTauLabelNew = cms.InputTag('hpsPFTauProducer'),
   postfix=""
)

# process.selectedPatTaus.cut = (
#     "pt > 8"
#    "pt > 8 && " +
#    "tauID('leadingTrackFinding') > 0.2 && tauID('byLooseIsolation') > 0.2"
# )

#  _____  ______ _____              _____                                      
# |  __ \|  ____|_   _|            / ____|                                     
# | |__) | |__    | |  ___  ___   | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
# |  ___/|  __|   | | / __|/ _ \   \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
# | |    | |     _| |_\__ \ (_) |  ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# |_|    |_|    |_____|___/\___/  |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                                                 | |                          
#                                                 |_|                          

# Select good muons to remove from cone
from WWAnalysis.AnalysisStep.wwMuons_cfi import MUON_ID_CUT, MUON_IP_CUT
process.goodMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("preBoostedMuons"),
    cut = cms.string( 'pt > 10 && ' + MUON_ID_CUT + " && " + MUON_IP_CUT )
)

# Select good electrons to remove from cone
from WWAnalysis.AnalysisStep.electronIDs_cff import ELE_NOCONV, ELE_IP, ELE_ID_LH_90_2011
process.goodElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("preBoostedElectrons"),
    cut = cms.string( 'pt > 10 && ' + ELE_NOCONV + " && " + ELE_IP + " && " + ELE_ID_LH_90_2011 )
)


# create isolation 'deposits'
process.pfIsoChargedHadronsDZ = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 2212, -2212),
    src = cms.InputTag("particleFlow")
)
process.pfIsoNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUp")
)
process.pfIsoChargedHadrons = process.pfIsoNeutralHadrons.clone ( pdgId = [211, -211, 321, -321, 999211, 2212, -2212] )
process.pfIsoPhotons        = process.pfIsoNeutralHadrons.clone ( pdgId = [22] )

# make the actual IsoDeposits
process.isoDepMuonWithChargedIsoDZ = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("preBoostedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(0.1),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfIsoChargedHadronsDZ"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)

process.isoDepMuonWithChargedIso = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("preBoostedMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfIsoChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)
process.isoDepMuonWithNeutralIso = process.isoDepMuonWithChargedIso.clone()
process.isoDepMuonWithPhotonIso = process.isoDepMuonWithChargedIso.clone()
process.isoDepMuonWithNeutralIso.ExtractorPSet.inputCandView = "pfIsoNeutralHadrons"
process.isoDepMuonWithPhotonIso.ExtractorPSet.inputCandView = "pfIsoPhotons"

process.isoDepElectronWithChargedIsoDZ = process.isoDepMuonWithChargedIsoDZ.clone( src = "preBoostedElectrons" )
process.isoDepElectronWithChargedIso = process.isoDepMuonWithChargedIso.clone( src = "preBoostedElectrons" )
process.isoDepElectronWithNeutralIso = process.isoDepMuonWithNeutralIso.clone( src = "preBoostedElectrons" )
process.isoDepElectronWithPhotonIso  = process.isoDepMuonWithPhotonIso.clone( src = "preBoostedElectrons" )

# Mu vetos
muVetos       = cms.vstring('0.01')
elChargedVeto = cms.vstring('0.01')
elNeutralVeto = cms.vstring('0.07')
elPhotonVeto  = cms.vstring('RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)')

# insert them into the pat leptons
# ha, made you look, they are actually down below in the electron and muon sections

# make the crazy sequence
process.pfIsoSequence = cms.Sequence( 
    process.goodMuons +
    process.goodElectrons +
    process.pfIsoNeutralHadrons +
    process.pfIsoChargedHadrons +
    process.pfIsoChargedHadronsDZ +
    process.pfIsoPhotons * (
        process.isoDepMuonWithChargedIsoDZ +
        process.isoDepMuonWithChargedIso +
        process.isoDepMuonWithNeutralIso +
        process.isoDepMuonWithPhotonIso +
        process.isoDepElectronWithChargedIsoDZ +
        process.isoDepElectronWithChargedIso +
        process.isoDepElectronWithNeutralIso +
        process.isoDepElectronWithPhotonIso 
    )
)
    

#  _____                        _                _              
# / ____|                      | |              | |             
#| (___   ___  _ __ ___   ___  | |     ___ _ __ | |_ ___  _ __  
# \___ \ / _ \| '_ ` _ \ / _ \ | |    / _ \ '_ \| __/ _ \| '_ \ 
# ____) | (_) | | | | | |  __/ | |___|  __/ |_) | || (_) | | | |
#|_____/ \___/|_| |_| |_|\___| |______\___| .__/ \__\___/|_| |_|
#                                         | |                   
#                                         |_|                   
# ____                  _   _             
#|  _ \                | | (_)            
#| |_) | ___   ___  ___| |_ _ _ __   __ _ 
#|  _ < / _ \ / _ \/ __| __| | '_ \ / _` |
#| |_) | (_) | (_) \__ \ |_| | | | | (_| |
#|____/ \___/ \___/|___/\__|_|_| |_|\__, |
#                                    __/ |
#                                   |___/ 

# First boost to get the IP values
# Then boost to add the PF isolation and the 


# add track IP information?
process.load("WWAnalysis.AnalysisStep.leptonBoosting_cff")
process.preBoostedElectrons = process.boostedElectrons.clone( electronTag = cms.InputTag("cleanPatElectronsTriggerMatch") )
process.preBoostedMuons = process.boostedMuons.clone( muonTag = cms.InputTag("cleanPatMuonsTriggerMatch") )
process.patDefaultSequence += process.preBoostedElectrons
process.patDefaultSequence += process.preBoostedMuons

if doPF2PATAlso:
    print "========================================================="
    print "__          __     _____  _   _ _____ _   _  _____ _ _ _ "
    print "\ \        / /\   |  __ \| \ | |_   _| \ | |/ ____| | | |"
    print " \ \  /\  / /  \  | |__) |  \| | | | |  \| | |  __| | | |"
    print "  \ \/  \/ / /\ \ |  _  /| . ` | | | | . ` | | |_ | | | |"
    print "   \  /\  / ____ \| | \ \| |\  |_| |_| |\  | |__| |_|_|_|"
    print "    \/  \/_/    \_\_|  \_\_| \_|_____|_| \_|\_____(_|_|_)"
    print "========================================================="
    print "                                                         "
    print "The new pf based isolation hasn't been adapted for PF2PAT"
    print "                                                         "
    print "========================================================="
    print "__          __     _____  _   _ _____ _   _  _____ _ _ _ "
    print "\ \        / /\   |  __ \| \ | |_   _| \ | |/ ____| | | |"
    print " \ \  /\  / /  \  | |__) |  \| | | | |  \| | |  __| | | |"
    print "  \ \/  \/ / /\ \ |  _  /| . ` | | | | . ` | | |_ | | | |"
    print "   \  /\  / ____ \| | \ \| |\  |_| |_| |\  | |__| |_|_|_|"
    print "    \/  \/_/    \_\_|  \_\_| \_|_____|_| \_|\_____(_|_|_)"
    print "========================================================="
    process.preBoostedElectronsPFlow = process.boostedElectrons.clone( muonTag = cms.InputTag("cleanPatElectronsTriggerMatchPFlow") )
    process.preBoostedMuonsPFlow     = process.boostedMuons.clone( muonTag = cms.InputTag("cleanPatMuonsTriggerMatchPFlow") )
    process.patPF2PATSequencePFlow += process.preBoostedElectronsPFlow
    process.patPF2PATSequencePFlow += process.preBoostedMuonsPFlow


# run the iso deposit producer for hcal
process.eleIsoDepositHcalFromTowers.src = "preBoostedElectrons"
process.patDefaultSequence += process.eleIsoDepositHcalFromTowers

# this is a mess, i am doing this so as to keep the same final branch names
# basically i am swapping the boosted leptons and the iso added leptons but
# this stuff below, has to come after that stuff above
process.load("WWAnalysis.AnalysisStep.isoAdding_cff")
process.boostedElectrons = process.isoAddedElectrons.clone( electronTag = "preBoostedElectrons" )
process.boostedMuons = process.isoAddedMuons.clone( muonTag = "preBoostedMuons" )

# add hcal information in full cone
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string("hcalFull")
process.boostedElectrons.deposits[-1].deltaR = 0.3
process.boostedElectrons.deposits[-1].vetos = []

# add the pf isolation values
# muons
# ucsd guess
process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithChargedIso"
process.boostedMuons.deposits[-1].label = cms.string("pfCharged")
process.boostedMuons.deposits[-1].deltaR = 0.4
process.boostedMuons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedMuons.deposits[-1].vetos += [ veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elChargedVeto ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithNeutralIso"
process.boostedMuons.deposits[-1].label = cms.string("pfNeutral")
process.boostedMuons.deposits[-1].deltaR = 0.4
process.boostedMuons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedMuons.deposits[-1].vetos += [ veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elNeutralVeto ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithPhotonIso"
process.boostedMuons.deposits[-1].label = cms.string("pfPhoton")
process.boostedMuons.deposits[-1].deltaR = 0.4
process.boostedMuons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedMuons.deposits[-1].vetos += [ veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedMuons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elPhotonVeto ]

# electrons
# ucsd guess
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithChargedIso"
process.boostedElectrons.deposits[-1].label = cms.string("pfCharged")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedElectrons.deposits[-1].vetos += [ veto for veto in elChargedVeto ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elChargedVeto ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithNeutralIso"
process.boostedElectrons.deposits[-1].label = cms.string("pfNeutral")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedElectrons.deposits[-1].vetos += [ veto for veto in elNeutralVeto ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elNeutralVeto ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithPhotonIso"
process.boostedElectrons.deposits[-1].label = cms.string("pfPhoton")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'Threshold(0.5)' ]
process.boostedElectrons.deposits[-1].vetos += [ veto for veto in elPhotonVeto ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodMuons:'+veto for veto in muVetos ]
process.boostedElectrons.deposits[-1].vetos += [ 'goodElectrons:'+veto for veto in elPhotonVeto ]

process.patDefaultSequence += process.pfIsoSequence
process.patDefaultSequence += process.boostedElectrons
process.patDefaultSequence += process.boostedMuons

# muons
# smurf default
process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithChargedIsoDZ"
process.boostedMuons.deposits[-1].label = cms.string("smurfCharged")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos = [ '0.01' ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithNeutralIso"
process.boostedMuons.deposits[-1].label = cms.string("smurfNeutral")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos  = [ 'Threshold(1.0)' ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithPhotonIso"
process.boostedMuons.deposits[-1].label = cms.string("smurfPhoton")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos = [ ]

# electrons
# smurf default
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithChargedIsoDZ"
process.boostedElectrons.deposits[-1].label = cms.string("smurfCharged")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ '0.01' ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithNeutralIso"
process.boostedElectrons.deposits[-1].label = cms.string("smurfNeutral")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'Threshold(1.0)','0.07' ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithPhotonIso"
process.boostedElectrons.deposits[-1].label = cms.string("smurfPhoton")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)' ]

# muons
# smurf default no overlap
process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithChargedIsoDZ"
process.boostedMuons.deposits[-1].label = cms.string("smurfNoOverCharged")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos = [ '0.01' ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithNeutralIso"
process.boostedMuons.deposits[-1].label = cms.string("smurfNoOverNeutral")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos  = [ 'Threshold(1.0)' ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

process.boostedMuons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedMuons.deposits[-1].src = "isoDepMuonWithPhotonIso"
process.boostedMuons.deposits[-1].label = cms.string("smurfNoOverPhoton")
process.boostedMuons.deposits[-1].deltaR = 0.3
process.boostedMuons.deposits[-1].vetos = [ ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

# electrons
# smurf default no overlap
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithChargedIsoDZ"
process.boostedElectrons.deposits[-1].label = cms.string("smurfNoOverCharged")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ '0.01' ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithNeutralIso"
process.boostedElectrons.deposits[-1].label = cms.string("smurfNoOverNeutral")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'Threshold(1.0)','0.07' ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]

process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].src = "isoDepElectronWithPhotonIso"
process.boostedElectrons.deposits[-1].label = cms.string("smurfNoOverPhoton")
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits[-1].vetos  = [ 'RectangularEtaPhiVeto(-0.025,0.025,-0.5,0.5)' ]
process.boostedMuons.deposits[-1].vetos += [ 'goodMuons:0.01','goodElectrons:0.01' ]


process.patDefaultSequence += process.pfIsoSequence
process.patDefaultSequence += process.boostedElectrons
process.patDefaultSequence += process.boostedMuons

if doPF2PATAlso:
    process.eleIsoDepositHcalFromTowersPFlow = process.eleIsoDepositHcalFromTowers( src = "cleanPatElectronsTriggerMatchPFlow" )
    process.patPF2PATSequencePFlow += process.eleIsoDepositHcalFromTowersPFlow

    process.boostedElectronsPFlow = process.isoAddedElectrons.clone( electronTag = "preBoostedElectronsPFlow" )
    process.boostedMuonsPFlow = process.isoAddedElectrons.clone( electronTag = "preBoostedMuonsPFlow" )

    process.boostedElectronsPFlow.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
    process.boostedElectronsPFlow.deposits[-1].src = "eleIsoDepositHcalFromTowersPFlow"
    process.boostedElectronsPFlow.deposits[-1].label = cms.string("hcalFull")
    process.boostedElectronsPFlow.deposits[-1].deltaR = 0.3
    process.boostedElectronsPFlow.deposits[-1].vetos = []



#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            

process.out.outputCommands =  cms.untracked.vstring(
    'drop *',
    # Leptons
    'keep *_boostedElectrons*_*_*',
    'keep *_boostedMuons*_*_*',
    'keep *_cleanPatTausTriggerMatch*_*_*',
    # Jets
    'keep patJets_slimPatJetsTriggerMatch_*_*',
    'keep patJets_slimPatJetsTriggerMatchPFlow_*_*',
    'keep patJets_slimPatJetsTriggerMatchNoPU_*_*',
    'keep recoGenJets_patJets_genJets_*',
    'keep recoGenJets_patJetsPFlow_genJets_*',
    'keep recoGenJets_patJetsNoPU_genJets_*',
    'keep recoGenJets_selectedPatJets_genJets_*',
    'keep recoGenJets_selectedPatJetsPFlow_genJets_*',
    'keep recoGenJets_selectedPatJetsNoPU_genJets_*',
#     'keep patJets_slimPatJetsTriggerMatchCalo_*_*',
#     'keep patJets_slimPatJetsTriggerMatchJPT_*_*',
    # Tracking
    #'keep *_offlinePrimaryVertices_*_'+process.name_(),
    'keep *_offlinePrimaryVerticesWithBS_*_*',
    'keep *_offlineBeamSpot_*_*',
    # MET
    'keep *_tcMet_*_*',
    'keep *_met_*_*',
    'keep *_pfMet_*_*',
    # MC
    'keep *_prunedGen_*_*',
    'keep *_genMetTrue_*_*',
    'keep GenEventInfoProduct_generator__HLT',
    # Trigger
    'keep *_TriggerResults_*_*',
    'keep *_vertexMapProd_*_*',
    # Misc
    'keep *_addPileupInfo_*_*',
    'keep *_chargedMetProducer_*_*',
    'keep *_trackMetProducer_*_*',
#     'keep *_mergedSuperClusters_*_'+process.name_(),
    'keep *_kt6PF*_rho_'+process.name_(),
    # Debug info, usually commented out
    #'keep *_patTrigger_*_*',  
    #'keep *_l1extraParticles_*_*',  
)
if is41XRelease:
    process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_'+process.name_())
else:
    process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_*')

process.prePatSequence  = cms.Sequence( process.preLeptonSequence + process.preElectronSequence + process.preMuonSequence)
process.postPatSequence = cms.Sequence( process.autreSeq + process.chargedMetSeq )

if not is41XRelease:
    # In order to use the offline vertices with BS constratint everywhere 
    massSearchReplaceAnyInputTag(process.preYieldFilter,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.prePatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.postPatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    if doPF2PATAlso:
        massSearchReplaceAnyInputTag(process.patPF2PATSequencePFlow,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"))



process.scrap      = cms.Path( process.noscraping ) 
process.outpath    = cms.EndPath(process.out)

if  doPF2PATAlso:
    process.patPath = cms.Path( process.preYieldFilter + process.prePatSequence * process.patDefaultSequence * process.patPF2PATSequencePFlow * process.postPatSequence )
    process.fakPath = cms.Path( process.preFakeFilter + process.prePatSequence * process.patDefaultSequence * process.patPF2PATSequencePFlow * process.postPatSequence )
else:
    process.patPath = cms.Path( process.preYieldFilter + process.prePatSequence * process.patDefaultSequence * process.postPatSequence)
    process.fakPath = cms.Path( process.preFakeFilter + process.prePatSequence * process.patDefaultSequence * process.postPatSequence )

process.out.SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('patPath'))

# from WWAnalysis.SkimStep.skimTools import addIsolationInformation
# addIsolationInformation(process)

if not doFakeRates:
    process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)
elif doFakeRates == 'also':
    process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'patPath', 'fakPath' ))
    process.countOverlaps = cms.Path(process.preYieldFilter + process.preFakeFilter)
    process.schedule = cms.Schedule( process.patPath, process.fakPath, process.countOverlaps, process.scrap, process.outpath)
elif doFakeRates == 'only':
    process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'fakPath' ))
    process.schedule = cms.Schedule( process.fakPath, process.scrap, process.outpath)
    

