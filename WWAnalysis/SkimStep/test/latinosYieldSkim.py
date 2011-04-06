import FWCore.ParameterSet.Config as cms

process = cms.Process("Yield")

#  _____               _____                               _                
# |  __ \             |  __ \                             | |               
# | |__) |   _ _ __   | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
# |  _  / | | | '_ \  |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | \ \ |_| | | | | | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# |_|  \_\__,_|_| |_| |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
# 

#Change me depending on your needs
isMC = True
# isMC = False
# isMC = RMMEMC
doPF2PATAlso = True
# doPF2PATAlso = False
# doPF2PATAlso = RMMEPF2PAT
doGenFilter = False

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

#Options
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#Global Tag Stuff
process.GlobalTag.globaltag = 'START311_V2::All'
# process.GlobalTag.globaltag = 'RMMEGlobalTag'

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

#Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/nfs/bluearc/group/edm/hww/Winter10.Flat/hww.flat.root',
#        'file:/nfs/bluearc/group/edm/hww/Winter10.Flat/hww.flat.root',
#        'file:/data/mangano/MC/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV_Spring11_AOD.root'
#        'file:/nfs/bluearc/group/edm/hww/Winter10.Flat/hww.flat.root',
#        'file:/home/mangano/skim/CMSSW_4_1_3/src/workingDirPU/lowPU.root'
#        'file:/home/mangano/skim/CMSSW_4_1_3/src/workingDirPU/highPU.root'
#        'file:/home/mangano/skim/CMSSW_4_1_3/src/workingDirPU/veryHighPU.root'
#        'file:/data/mangano/MC/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV_Spring11_AOD.root'
#        'file:/data/mangano/MC/Spring11/TTJets_madgraph_Spring11_AOD.root'
#        'file:/data/mangano/MC/Spring11/WJets_madgraph_Spring11_AOD.root'
#         'file:/home/mwlebour/data/WW.38XMC.Samples/DYToEEM20CT10Z2powheg.root'
        'file:/home/mwlebour/data/Winter10/Hww160.root'
#        'RMMEFN'
    )
)

# Good if you have a bunch of files you want to run on
#from glob import glob
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
#process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/oct29Skim/WWFull/*/*.root')]

#You need an output module before calling any of the pat functions, changed down below
process.out = cms.OutputModule("PoolOutputModule",fileName=cms.untracked.string(''),outputCommands=cms.untracked.vstring( 'drop *',))


#  _____                 ______ _ _ _            
# |  __ \               |  ____(_) | |           
# | |__) | __ ___ ______| |__   _| | |_ ___ _ __ 
# |  ___/ '__/ _ \______|  __| | | | __/ _ \ '__|
# | |   | | |  __/      | |    | | | ||  __/ |   
# |_|   |_|  \___|      |_|    |_|_|\__\___|_|   
#                                                

process.nonSTAMuons = cms.EDFilter("MuonRefSelector",
    cut = cms.string("type!=8"),
    src = cms.InputTag("muons"),
    filter = cms.bool(False)
)

process.cleanRecoTaus = cms.EDFilter("PFTauSelector",
    src = cms.InputTag("hpsPFTauProducer"),
    discriminators = cms.VPSet(
        cms.PSet( 
          discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
          selectionCut=cms.double(0.5)
        )
    )
)

process.allLeps = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(
        cms.InputTag("gsfElectrons"), 
        cms.InputTag("nonSTAMuons"), 
        cms.InputTag("cleanRecoTaus")
    )
)

process.noTauLeps = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(
        cms.InputTag("gsfElectrons"), 
        cms.InputTag("nonSTAMuons"), 
    )
)

process.allDiLep = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('allLeps noTauLeps'),
    cut = cms.string(
        'deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.05 && ' + 
        'min(daughter(0).pt,daughter(1).pt) >  8 && ' +
        'max(daughter(0).pt,daughter(1).pt) > 18'
    ),
    checkCharge = cms.bool(False)
)

process.countDiLeps  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("allDiLep"), minNumber = cms.uint32(1))
process.preFilter = cms.Sequence( process.nonSTAMuons * process.cleanRecoTaus * process.allLeps * process.noTauLeps * process.allDiLep * process.countDiLeps )



#  _____            _____           _   _      _        ______ _ _ _            
# / ____|          |  __ \         | | (_)    | |      |  ____(_) | |           
#| |  __  ___ _ __ | |__) |_ _ _ __| |_ _  ___| | ___  | |__   _| | |_ ___ _ __ 
#| | |_ |/ _ \ '_ \|  ___/ _` | '__| __| |/ __| |/ _ \ |  __| | | | __/ _ \ '__|
#| |__| |  __/ | | | |  | (_| | |  | |_| | (__| |  __/ | |    | | | ||  __/ |   
# \_____|\___|_| |_|_|   \__,_|_|   \__|_|\___|_|\___| |_|    |_|_|\__\___|_|   
        
process.genLepFromW10 = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    filter = cms.bool(True),
    cut = cms.string("(abs(pdgId)==11 || abs(pdgId)==13) && abs(mother.mother.pdgId)==24 &&" +
                     "abs(eta)<2.5 && pt>10.")
)

process.genLep10CountFilter=cms.EDFilter("CandViewCountFilter", 
                  src =cms.InputTag("genLepFromW10"), 
                  minNumber = cms.uint32(2)
)

process.genLepFromW20 = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genLepFromW10"),
    filter = cms.bool(True),
    cut = cms.string("pt>20.")
)

process.genFilter = cms.Sequence(process.genLepFromW10*process.genLep10CountFilter*process.genLepFromW20)

if isMC:
    if doGenFilter:                                
        process.preFilter.replace(
            process.nonSTAMuons,
            process.nonSTAMuons*
            process.genFilter
        )


#  _____     _______    _____                                      
# |  __ \ /\|__   __|  / ____|                                     
# | |__) /  \  | |    | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
# |  ___/ /\ \ | |     \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
# | |  / ____ \| |     ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# |_| /_/    \_\_|    |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                                     | |                          
#                                     |_|    

process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.trigTools import *

#  _______   _                         ______ _ _ _            
# |__   __| (_)                       |  ____(_) | |           
#    | |_ __ _  __ _  __ _  ___ _ __  | |__   _| | |_ ___ _ __ 
#    | | '__| |/ _` |/ _` |/ _ \ '__| |  __| | | | __/ _ \ '__|
#    | | |  | | (_| | (_| |  __/ |    | |    | | | ||  __/ |   
#    |_|_|  |_|\__, |\__, |\___|_|    |_|    |_|_|\__\___|_|   
#               __/ | __/ |                                    
#              |___/ |___/                                  
# 

process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerMatcher_cfi" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi" )

tempProd = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path( "HLT_Mu20_v*" )'),
    src = cms.InputTag("cleanPatMuons"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)


eleTriggers = [
    "HLT_Ele8_CaloIdL_CaloIsoVL_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
    "HLT_Ele17_CaloIdL_CaloIsoVL_v*",
    "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
    "HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*",
    "HLT_Ele45_CaloIdVT_TrkIdT_v*",
    "HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v*",
    "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v*",
    "HLT_Mu17_Ele8_CaloIdL_v*",
    "HLT_Mu8_Ele17_CaloIdL_v*",
    "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*",
    "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*",
    "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*",
    "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v*",
]
eleTriggerModules = dict(zip([ "cleanElectronTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in eleTriggers ],eleTriggers))
for key in eleTriggerModules:
    setattr(process,key,tempProd.clone(src = "cleanPatElectrons", matchedCuts = 'path("{0}")'.format(eleTriggerModules[key])))

muTriggers = [
    "HLT_Mu5_v*",
    "HLT_Mu8_v*",
    "HLT_Mu12_v*",
    "HLT_Mu15_v*",
    "HLT_Mu20_v*",
    "HLT_Mu24_v*",
    "HLT_Mu30_v*",
    "HLT_IsoMu12_v*",
    "HLT_IsoMu15_v*",
    "HLT_IsoMu17_v*",
    "HLT_IsoMu24_v*",
    "HLT_IsoMu30_v*",
    "HLT_Mu17_Ele8_CaloIdL_v*",
    "HLT_Mu8_Ele17_CaloIdL_v*",
    "HLT_Mu8_Jet40_v*",
    "HLT_Mu8_Photon20_CaloIdVT_IsoT_v*",
    "HLT_IsoMu12_LooseIsoPFTau10_v*",
    "HLT_Mu15_LooseIsoPFTau20_v*",
    "HLT_DoubleMu7_v*",
]
muTriggerModules = dict(zip([ "cleanMuonTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in muTriggers ],muTriggers))
for key in muTriggerModules:
    setattr(process,key,tempProd.clone(src = "cleanPatMuons", matchedCuts = 'path("{0}")'.format(muTriggerModules[key])))

tauTriggers = [
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*",
    "HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v*",
    "HLT_Mu15_LooseIsoPFTau20_v*",
    "HLT_IsoMu12_LooseIsoPFTau10_v*",
    "HLT_IsoMu15_LooseIsoPFTau20_v*",
]
tauTriggerModules = dict(zip([ "cleanTauTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in tauTriggers ],tauTriggers))
for key in tauTriggerModules:
    setattr(process,key,tempProd.clone(src = "cleanPatTaus", matchedCuts = 'path("{0}")'.format(tauTriggerModules[key])))

myDefaultTriggerMatchers = eleTriggerModules.keys()[:] + muTriggerModules.keys()[:] + tauTriggerModules.keys()[:] + [
    'cleanPhotonTriggerMatchHLTPhoton26IsoVLPhoton18',
    'cleanJetTriggerMatchHLTJet240',
    'metTriggerMatchHLTMET100',
]

switchOnTriggerMatchEmbedding(process,triggerMatchers=myDefaultTriggerMatchers)

process.patTrigger.onlyStandAlone = True
process.patTrigger.processName  = '*' 

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

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
process.kt6PFJets.Rho_EtaMax = cms.double(5.0)

process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)
process.kt6PFJetsForIso.Ghost_EtaMax = cms.double(2.5)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)

process.load("WWAnalysis.SkimStep.cutsInCategoriesHWWElectronIdentificationV04_cfi")
process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
process.egammaIDLikelihood = process.eidLikelihoodExt.clone()

process.eIdSequence = cms.Sequence(
    process.eidVeryLoose +
    process.eidLoose +
    process.eidMedium +
    process.eidTight +
    process.eidSuperTight +
    process.eidHyperTight1 +
    process.eidHyperTight2 +
    process.eidHyperTight3 +
    process.eidVeryLooseMC +
    process.eidLooseMC +
    process.eidMediumMC +
    process.eidTightMC +
    process.eidSuperTightMC +
    process.eidHyperTight1MC +
    process.eidHyperTight2MC +
    process.eidHyperTight3MC +
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

process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho",process.name_()))
process.rhoEl = process.rhoMu.clone(leptonTag = "gsfElectrons")

process.valueMaps = cms.Sequence(
    process.betaMu +
    process.betaEl +
    process.rhoMu +
    process.rhoEl 
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

if isMC:
    process.preLeptonSequence = cms.Sequence(process.prunedGen)
else:
    process.preLeptonSequence = cms.Sequence()

process.preLeptonSequence += (process.kt6PFJets * process.kt6PFJetsForIso * process.valueMaps +
                              process.offlinePrimaryVertices +
                              process.eIdSequence )





#  ______ _           _                     _____      _   _      
# |  ____| |         | |                   |  __ \    | | | |     
# | |__  | | ___  ___| |_ _ __ ___  _ __   | |__) |_ _| |_| |__   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \  |  ___/ _` | __| '_ \  
# | |____| |  __/ (__| |_| | | (_) | | | | | |  | (_| | |_| | | | 
# |______|_|\___|\___|\__|_|  \___/|_| |_| |_|   \__,_|\__|_| |_| 
#                                                                 

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
)

#Set the Pat Electrons to use the eID
for module in listModules(process.eIdSequence):
    setattr(process.patElectrons.electronIDSources,module.label(),cms.InputTag(module.label()))

process.load("WWAnalysis.Tools.convValueMapProd_cfi")
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
)

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



#  _   _                             _     _ 
# | \ | |                   /\      | |   | |
# |  \| | _____      __    /  \   __| | __| |
# | . ` |/ _ \ \ /\ / /   / /\ \ / _` |/ _` |
# | |\  | (_) \ V  V /   / ____ \ (_| | (_| |
# |_| \_|\___/ \_/\_/   /_/    \_\__,_|\__,_|
#                                               
#  _____  ______ ___  _____     _______ 
# |  __ \|  ____|__ \|  __ \ /\|__   __|
# | |__) | |__     ) | |__) /  \  | |   
# |  ___/|  __|   / /|  ___/ /\ \ | |   
# | |    | |     / /_| |  / ____ \| |   
# |_|    |_|    |____|_| /_/    \_\_|   
#                                       
def addFastJetCorrection(process,label,seq="patDefaultSequence"):
    corrFact = getattr(process,"patJetCorrFactors"+label)
    setattr(process,"patJetCorrFactorsFastJet"+label,corrFact.clone())
    getattr(process,"patJetCorrFactorsFastJet"+label).levels[0] = 'L1FastJet'
    getattr(process,seq).replace(
        getattr(process,"patJetCorrFactors"+label),
        getattr(process,"patJetCorrFactors"+label) +
        getattr(process,"patJetCorrFactorsFastJet"+label) 
    )
    getattr(process,"patJets"+label).jetCorrFactorsSource.append( cms.InputTag("patJetCorrFactorsFastJet"+label) )
  


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
    process.pfJetsPFlow.Rho_EtaMax = cms.double(5.0)
    # Turn on secondary JEC w/ FastJet
    addFastJetCorrection(process,"PFlow","patPF2PATSequencePFlow")

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

# Re-cluster ak5PFJets w/ Area calculation on
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(5.0)

# Re-cluster jets w/ pfNoPileUp
process.ak5PFJetsNoPU = process.ak5PFJets.clone(    
    src =   cms.InputTag("pfNoPileUp")  
)
process.patDefaultSequence.replace(process.pfNoPileUp,process.ak5PFJets + process.pfNoPileUp*process.ak5PFJetsNoPU)


#all the other jets:
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
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
    jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
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

process.cleanJetTriggerMatchHLTJet240NoPU = process.cleanJetTriggerMatchHLTJet240.clone( src = cms.InputTag("cleanPatJetsNoPU") )
process.patDefaultSequence.replace(
    process.cleanJetTriggerMatchHLTJet240,
    process.cleanJetTriggerMatchHLTJet240 + 
    process.cleanJetTriggerMatchHLTJet240NoPU 
)

process.cleanPatJetsTriggerMatchNoPU = process.cleanPatJetsTriggerMatch.clone( matches = ["cleanJetTriggerMatchHLTJet240NoPU"], src = cms.InputTag("cleanPatJetsNoPU") )
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
addFastJetCorrection(process,"NoPU")

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
    src = cms.InputTag("cleanPatMuonsTriggerMatch"),
    filter = cms.bool(False)
)

process.lepsForMET = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag(cms.InputTag("patMuonsWithTriggerNoSA"), cms.InputTag("cleanPatElectronsTriggerMatch"))
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


# _____                        _                _              
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

process.load("WWAnalysis.AnalysisStep.leptonBoosting_cff")
process.boostedElectrons = process.boostedElectrons.clone( muonTag = cms.untracked.InputTag("cleanPatElectronsTriggerMatch") )
process.boostedMuons = process.boostedMuons.clone( muonTag = cms.untracked.InputTag("cleanPatMuonsTriggerMatch") )
process.patDefaultSequence += process.boostedElectrons
process.patDefaultSequence += process.boostedMuons

process.boostedElectronsPFlow = process.boostedElectrons.clone( muonTag = cms.untracked.InputTag("cleanPatElectronsTriggerMatchPFlow") )
process.boostedMuonsPFlow     = process.boostedMuons.clone( muonTag = cms.untracked.InputTag("cleanPatMuonsTriggerMatchPFlow") )
process.patPF2PATSequencePFlow += process.boostedElectronsPFlow
process.patPF2PATSequencePFlow += process.boostedMuonsPFlow




#   _____      _              _       _      
#  / ____|    | |            | |     | |     
# | (___   ___| |__   ___  __| |_   _| | ___ 
#  \___ \ / __| '_ \ / _ \/ _` | | | | |/ _ \
#  ____) | (__| | | |  __/ (_| | |_| | |  __/
# |_____/ \___|_| |_|\___|\__,_|\__,_|_|\___|
#                                            

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('latinosYieldSkim.root'),
#     fileName = cms.untracked.string('RMMEFN'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
        # Leptons
        'keep *_boostedElectrons*_*_*',
        'keep *_boostedMuons*_*_*',
        'keep *_cleanPatTausTriggerMatch*_*_*',
        # Jets
        'keep patJets_slimPatJetsTriggerMatch_*_*',
        'keep patJets_slimPatJetsTriggerMatchPFlow_*_*',
#         'keep patJets_slimPatJetsTriggerMatchCalo_*_*',
        'keep patJets_slimPatJetsTriggerMatchNoPU_*_*',
#         'keep patJets_slimPatJetsTriggerMatchJPT_*_*',
        # Tracking
        'keep *_offlinePrimaryVertices_*_'+process.name_(),
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
#         'keep *_mergedSuperClusters_*_'+process.name_(),
        'keep *_kt6PFJetsForIso_rho_'+process.name_(),
        'keep *_kt6PFJets_rho_'+process.name_(),
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('patPath' )),
)


process.prePatSequence  = cms.Sequence( process.preLeptonSequence + process.preElectronSequence + process.preMuonSequence )
process.postPatSequence = cms.Sequence( process.autreSeq + process.chargedMetSeq )

process.scrap      = cms.Path( process.noscraping ) 
process.outpath    = cms.EndPath(process.out)

if  doPF2PATAlso:
    process.patPath = cms.Path( process.preFilter * process.prePatSequence * process.patDefaultSequence * process.patPF2PATSequencePFlow * process.postPatSequence)
else:
    process.patPath = cms.Path( process.preFilter * process.prePatSequence * process.patDefaultSequence * process.postPatSequence)

from WWAnalysis.SkimStep.skimTools import addIsolationInformation
addIsolationInformation(process)

process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)

