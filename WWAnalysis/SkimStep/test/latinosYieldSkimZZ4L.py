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

isMC = True
process.GlobalTag.globaltag = 'START42_V14B::All'
doBorisGenFilter = False
isVV = False

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
process.source.fileNames = ['file:HToZZTo4L_M-120_Fall11S6.00215E21D5C4.root'] #on lxplus423.cern.ch

process.out = cms.OutputModule("PoolOutputModule", outputCommands =  cms.untracked.vstring(), fileName = cms.untracked.string('latinosYieldSkim.root') )

process.maxEvents.input = 1000

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

process.load('WWAnalysis.SkimStep.vertexFiltering_cff')
process.preLeptonSequence += process.firstVertexIsGood
process.preLeptonSequence += process.goodPrimaryVertices




# Rho calculations
from WWAnalysis.SkimStep.rhoCalculations_cff import addRhoVariables
addRhoVariables(process,process.preLeptonSequence)

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
    cms.InputTag("eleSmurfPF"),
    cms.InputTag("rhoEl")
)

#Set the Pat Electrons to use the eID
for module in eidModules:
    setattr(process.patElectrons.electronIDSources,module.label(),cms.InputTag(module.label()))

process.load("WWAnalysis.Tools.electronPFIsoMapProd_cfi")
process.eleSmurfPF = process.electronPFIsoMapProd.clone()
process.preElectronSequence = cms.Sequence(process.eleSmurfPF)




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
    cms.InputTag("muSmurfPF"),
    cms.InputTag("rhoMu"),
)
process.patMuons.userData.userInts.src = cms.VInputTag(
    cms.InputTag("pfMuId"),
)
process.patMuons.isolationValues = cms.PSet()
process.muonMatch.matched = "prunedGen"
process.load("WWAnalysis.Tools.muonPFIsoMapProd_cfi")
process.muSmurfPF = process.muonPFIsoMapProd.clone()
process.load("WWAnalysis.Tools.pfMuId_cfi")
process.preMuonSequence = cms.Sequence(
    process.muSmurfPF +
    process.pfMuId
)

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
    myCorrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute')
else:
    myCorrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

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

# only keep em above 7 GeV as the f'in smurfs
process.selectedPatJets.cut = "correctedJet('Uncorrected').pt > 7"


# Some stuff to save space
process.patJets.embedCaloTowers = False
process.patJets.addTagInfos = False
process.patJets.embedPFCandidates = False
process.patJets.addAssociatedTracks = False

process.boostedPatJets = cms.EDProducer("PatJetBooster",
    jetTag = cms.InputTag("cleanPatJets"),
    vertexTag = cms.InputTag("goodPrimaryVertices"),
)

process.slimPatJets = cms.EDProducer("PATJetSlimmer",
    src = cms.InputTag("boostedPatJets"),
    clearJetVars = cms.bool(True),
    clearDaughters = cms.bool(True),
    dropSpecific = cms.bool(False),
)

process.patDefaultSequence += (
    process.boostedPatJets *
    process.slimPatJets
)



#               _               _____      _ _           _   _                 
#    /\        | |             / ____|    | | |         | | (_)                
#   /  \  _   _| |_ _ __ ___  | |     ___ | | | ___  ___| |_ _  ___  _ __  ___ 
#  / /\ \| | | | __| '__/ _ \ | |    / _ \| | |/ _ \/ __| __| |/ _ \| '_ \/ __|
# / ____ \ |_| | |_| | |  __/ | |___| (_) | | |  __/ (__| |_| | (_) | | | \__ \
#/_/    \_\__,_|\__|_|  \___|  \_____\___/|_|_|\___|\___|\__|_|\___/|_| |_|___/
#                                                                              

process.reducedPFCands = cms.EDProducer("ReducedCandidatesProducer",
    srcCands = cms.InputTag("particleFlow",""),
    srcVertices = cms.InputTag("goodPrimaryVertices"),
    dz = cms.double(0.1),
    ptThresh = cms.double(0.5),
)

process.autreSeq = cms.Sequence(process.reducedPFCands)

process.load("WWAnalysis.SkimStep.hzz4lDetectorIsolation_cff")

# make the crazy sequence
process.hzzIsoSequence = cms.Sequence(process.hzz4lDetectorIsolationSequence)
    

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

process.load("WWAnalysis.AnalysisStep.isoAdding_cff")
process.boostedElectrons = process.isoAddedElectrons.clone( electronTag = "preBoostedElectrons" )
process.boostedMuons = process.isoAddedMuons.clone( muonTag = "preBoostedMuons" )
from WWAnalysis.SkimStep.hzz4lDetectorIsolation_cff import muIsoFromDepsZZ4L, eleIsoFromDepsZZ4L
process.boostedMuons.deposits     += muIsoFromDepsZZ4L
process.boostedElectrons.deposits += eleIsoFromDepsZZ4L

process.patDefaultSequence += process.hzzIsoSequence
process.patDefaultSequence += process.boostedElectrons
process.patDefaultSequence += process.boostedMuons


# Setting up PAT photons
process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = True


#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            

process.out.outputCommands =  cms.untracked.vstring(
    'drop *',
    'keep *_prunedGen_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_addPileupInfo_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_goodPrimaryVertices_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_reducedPFCands_*_*',
    'keep *_boostedElectrons*_*_*',
    'keep *_boostedMuons*_*_*',
    'keep *_cleanPatPhotons_*_*',
    'keep *_correctedHybridSuperClusters_*_*',
    'keep *_correctedMulti5x5SuperClustersWithPreshower_*_*',
    'keep *_slimPatJets_*_*',
    'keep *_pfMet_*_*',
    'keep *_kt6PF*_rho_*',

)

process.prePatSequence  = cms.Sequence( process.preLeptonSequence + process.preElectronSequence + process.preMuonSequence )
process.postPatSequence = cms.Sequence( process.autreSeq )


# In order to use the good primary vertices everywhere (It would be nicer to set the proper inputTags in the first place)
massSearchReplaceAnyInputTag(process.prePatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
massSearchReplaceAnyInputTag(process.postPatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)
process.firstVertexIsGood.vertices = cms.InputTag("offlinePrimaryVertices")
process.goodPrimaryVertices.src = cms.InputTag("offlinePrimaryVertices")

process.scrap      = cms.Path( process.noscraping ) 
process.outpath    = cms.EndPath(process.out)

process.patPath = cms.Path( process.prePatSequence * process.patDefaultSequence * process.postPatSequence)

process.out.SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('patPath'))

process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)

