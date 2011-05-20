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
isMC = RMMEMC
# isMC = True
# isMC = False
# doPF2PATAlso = RMMEPF2PAT
# doPF2PATAlso = True
doPF2PATAlso = False
doGenFilter = False
is41XRelease = False

doFakeRates = RMMEFAKE # 'only', 'also' or None
# doFakeRates = None # 'only', 'also' or None
# doFakeRates = 'only' # 'only', 'also' or None


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

#Global Tag Stuff
process.GlobalTag.globaltag = 'RMMEGlobalTag'
# process.GlobalTag.globaltag = 'START311_V2::All'
# process.GlobalTag.globaltag = 'GR_R_311_V2::All'

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

#Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'RMMEFN'
#        'file:/nfs/bluearc/group/edm/hww/Spring11.Flat/GluGluToHToWWTo2L2Nu_M-160_7TeV_Spring11_AOD.root'
#       ##  some real data from good run 163339
#	'/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v2/000/163/339/D0E8C932-7F70-E011-8137-001D09F252E9.root',  # events:   20957
#	'/store/data/Run2011A/DoubleMu/AOD/PromptReco-v2/000/163/339/DA578B6F-2B70-E011-8A36-003048F0258C.root',        # events:   21840
#	'/store/data/Run2011A/SingleMu/AOD/PromptReco-v2/000/163/339/F4E034EB-2170-E011-B8A2-001617E30D52.root',        # events:   25917
#	'/store/data/Run2011A/SingleElectron/AOD/PromptReco-v2/000/163/339/B4D937BB-5E70-E011-A3FB-001617C3B654.root',  # events:   23922
#	'/store/data/Run2011A/MuEG/AOD/PromptReco-v2/000/163/339/5E40BD3B-7F70-E011-805B-001D09F26C5C.root',            # events:   19775
        #'file:/data/gpetrucc/7TeV/hwww/DYToMuMu_Spring11_AODSIM.root'
        #'file:/data/gpetrucc/7TeV/hwww/SingleMu_2011Av2_r163339.root'
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
process.preFilter = cms.Sequence( process.nonSTAMuons * (process.cleanRecoTaus * process.allLeps + process.noTauLeps) * process.allDiLep * process.countDiLeps )

### Pre-filters for fake rate studies
from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter
process.hltFilter4Fakes = triggerResultsFilter.clone(
    l1tResults = '',
    hltResults = cms.InputTag( "TriggerResults", "", "HLT"),
    throw = True,
    triggerConditions = [ 'HLT_Mu8_v*',
                          'HLT_Mu15_v*',
                          'HLT_Mu24_v*',
                          'HLT_Ele8_v*',
                          'HLT_Ele8_CaloIdL_CaloIsoVL_v*',
                          'HLT_Ele17_CaloIdL_CaloIsoVL_v*' ]
)

process.jetsPt15 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("ak5PFJets"),
    cut = cms.string("pt > 15 && numberOfDaughters >= 3"),
    filter = cms.bool(True),
)
process.leptonPlusJet = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('noTauLeps jetsPt15'),
    cut = cms.string(
        'deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 0.7 && ' +
        'daughter(0).pt >= 10' 
    ),
    checkCharge = cms.bool(False),
)
process.countLeptonPlusJet = process.countDiLeps.clone(src = "leptonPlusJet")

process.diLeptons4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('noTauLeps@+ noTauLeps@-'),
    cut = cms.string('abs(daughter(0).pdgId) == abs(daughter(1).pdgId) && (mass < 12 || abs(mass-91.1876) > 15)')
)
process.diLeptons4VetoFilter = process.countDiLeps.clone(src = "diLeptons4Veto")

process.metVeto20 = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("pfMet"),
    cut = cms.string("pt < 20"),
    filter = cms.bool(True),
)
process.recoW4Veto = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('noTauLeps pfMet'),
    cut = cms.string("sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi))) > 20"),
    checkCharge = cms.bool(False)
)
process.recoW4VetoFilter = process.countDiLeps.clone(src = "recoW4Veto")

process.pre4Fakes = cms.Sequence( 
    process.hltFilter4Fakes +
    process.metVeto20 +
    process.nonSTAMuons * process.noTauLeps *
    ( process.jetsPt15 * process.leptonPlusJet * process.countLeptonPlusJet +
      process.diLeptons4Veto * ~process.diLeptons4VetoFilter +
      process.recoW4Veto     * ~process.recoW4VetoFilter )
)

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

process.patL1Trigger = cms.EDProducer("L1Extra2PAT", src = cms.InputTag("l1extraParticles"))
process.patDefaultSequence.replace(process.patTrigger, process.patTrigger + process.patL1Trigger)

tempProd = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
    matchedCuts = cms.string('path( "HLT_Mu20_v*" )'),
    src = cms.InputTag("cleanPatMuons"),
    maxDPtRel = cms.double(0.5),
    resolveByMatchQuality = cms.bool(True),
    maxDeltaR = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("patTrigger")
)
tempProdNoPt = tempProd.clone()
setattr(tempProdNoPt, '_TypedParameterizable__type', 'PATTriggerMatcherDRLessByR')
delattr(tempProdNoPt, 'maxDPtRel')

## List of lists of exclusive e/gamma trigger objetcs
eleTriggerColls = [ 
    [ 'hltL1IsoRecoEcalCandidate',       'hltL1NonIsoRecoEcalCandidate' ],
    [ 'hltPixelMatchElectronsL1Iso',     'hltPixelMatchElectronsL1NonIso' ],
    [ 'hltPixelMatch3HitElectronsL1Iso', 'hltPixelMatch3HitElectronsL1NonIso' ],
    [ 'hltPixelMatchElectronsActivity' ] ,
    [ 'hltPixelMatch3HitElectronsActivity' ] ,
    [ 'hltRecoEcalSuperClusterActivityCandidate' ] ,
]
eleTriggers = [
    "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*",
    "HLT_Mu11_Ele8_v*",
    "HLT_Mu5_Ele17_v*",
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
    "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*",
    "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v*",
    "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v*",
    "HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v*",
    "HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v*",
]
eleTriggerModules = dict(zip([ "cleanElectronTriggerMatch{0}".format(k.replace('v*','').replace('HLT_','').replace('_','')) for k in eleTriggers ],eleTriggers))
for key in eleTriggerModules:
    setattr(process,key,tempProd.clone(src = "cleanPatElectrons", matchedCuts = 'path("{0}")'.format(eleTriggerModules[key])))

eleTriggerCollModules = dict(zip([ "cleanElectronTriggerMatch{0}".format(k[0]) for k in eleTriggerColls ],eleTriggerColls))
for key in eleTriggerCollModules:
    setattr(process,key,tempProd.clone(src = "cleanPatElectrons", matchedCuts = " || ".join(['coll("%s")' % c for c in eleTriggerCollModules[key]])))

process.cleanElectronTriggerMatchL1 = tempProdNoPt.clone(
    src = "cleanPatElectrons",
    matched = "patL1Trigger",
    matchedCuts = "coll('l1extraParticles:Isolated') || coll('l1extraParticles:NonIsolated')",
)

jetTrigMatches = [ "cleanJetTriggerMatchHLTJet240", "cleanJetTriggerMatchL3Mu", "cleanJetTriggerMatchL1EG" ]   # HLTJet240 is from PAT default, I suppose
process.cleanJetTriggerMatchL1EG = tempProdNoPt.clone(
    src = "cleanPatJets",
    matched = "patL1Trigger",
    matchedCuts = "coll('l1extraParticles:Isolated') || coll('l1extraParticles:NonIsolated')",
)
process.cleanJetTriggerMatchL3Mu = tempProdNoPt.clone(
    src = "cleanPatJets",
    matched = "patTrigger",
    matchedCuts = "coll('hltL3MuonCandidates')",
    maxDeltaR = 0.5,
)

muTriggers = [
    "HLT_Mu5_v*",
    "HLT_Mu8_v*",
    "HLT_Mu12_v*",
    "HLT_Mu15_v*",
    "HLT_Mu20_v*",
    "HLT_Mu21_v*",
    "HLT_Mu24_v*",
    "HLT_Mu30_v*",
    "HLT_IsoMu12_v*",
    "HLT_IsoMu15_v*",
    "HLT_IsoMu17_v*",
    "HLT_IsoMu24_v*",
    "HLT_IsoMu30_v*",
    "HLT_Mu11_Ele8_v*",
    "HLT_Mu5_Ele17_v*",
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

process.cleanMuonTriggerMatchByObject   = tempProd.clone(src = "cleanPatMuons", matchedCuts = 'coll("hltL3MuonCandidates")')
process.cleanMuonTriggerMatchByL2Object = tempProdNoPt.clone(src = "cleanPatMuons", matchedCuts = 'coll("hltL2MuonCandidates")', maxDeltaR=0.7)

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

myDefaultTriggerMatchers = eleTriggerModules.keys()[:] + eleTriggerCollModules.keys()[:] + muTriggerModules.keys()[:] + tauTriggerModules.keys()[:] + jetTrigMatches + [
    'cleanMuonTriggerMatchByObject',
    'cleanMuonTriggerMatchByL2Object',
    'cleanPhotonTriggerMatchHLTPhoton26IsoVLPhoton18',
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

if is41XRelease:
    process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
    process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
    process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
    process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
    process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)

#  _____             _____  _           _             
# |  __ \           |  __ \| |         (_)            
# | |__) |___ ______| |__) | |__   ___  _ _ __   __ _ 
# |  _  // _ \______|  _  /| '_ \ / _ \| | '_ \ / _` |
# | | \ \  __/      | | \ \| | | | (_) | | | | | (_| |
# |_|  \_\___|      |_|  \_\_| |_|\___/|_|_| |_|\__, |
#                                                __/ |
#                                               |___/ 
# 
  
process.load('RecoJets.JetProducers.kt4PFJets_cfi')

process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(4.5)

process.kt6PFJetsForIso = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)
process.kt6PFJetsForIso.Ghost_EtaMax = cms.double(2.5)

process.kt6PFJetsNoPU = process.kt6PFJets.clone( src = "pfNoPileUp" )
process.kt6PFJetsForIsoNoPU = process.kt6PFJetsForIso.clone( src = "pfNoPileUp" )

# Re-cluster ak5PFJets w/ Area calculation on
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.5)

# Re-cluster jets w/ pfNoPileUp
process.ak5PFJetsNoPU = process.ak5PFJets.clone( src = "pfNoPileUp" )

process.load("WWAnalysis.Tools.betaValueMapProducer_cfi")
process.load("WWAnalysis.Tools.rhoValueMapProducer_cfi")

process.betaMu = process.betaValueMapProducer.clone()
process.betaEl = process.betaValueMapProducer.clone()
process.betaEl.leptonTag = "gsfElectrons"
process.betaEl.dEtaVeto = 0.015
process.betaEl.dRVeto = 0.0

process.rhoMu = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIso","rho",process.name_()), leptonTag = "muons")
process.rhoEl = process.rhoMu.clone(leptonTag = "gsfElectrons")
process.rhoMuNoPU = process.rhoValueMapProducer.clone(rhoTag = cms.untracked.InputTag("kt6PFJetsForIsoNoPU","rho",process.name_()), leptonTag = "muons")
process.rhoElNoPU = process.rhoMuNoPU.clone(leptonTag = "gsfElectrons")

process.valueMaps = cms.Sequence(
    process.betaMu +
    process.betaEl +
    process.rhoMu +
    process.rhoEl +
    process.rhoMuNoPU +
    process.rhoElNoPU 
)
    

#  ______ _____ _____     _____                                      
# |  ____|_   _|  __ \   / ____|                                     
# | |__    | | | |  | | | (___   ___  __ _ _   _  ___ _ __   ___ ___ 
# |  __|   | | | |  | |  \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \
# | |____ _| |_| |__| |  ____) |  __/ (_| | |_| |  __/ | | | (_|  __/
# |______|_____|_____/  |_____/ \___|\__, |\__,_|\___|_| |_|\___\___|
#                                       | |                          
#                                       |_|                         

process.eIdSequence = cms.Sequence()

from WWAnalysis.SkimStep.simpleCutBasedElectronIDSpring11_cfi import simpleCutBasedElectronID
process.vbtf11WP60 = simpleCutBasedElectronID.clone( electronQuality = '60' )
process.vbtf11WP70 = simpleCutBasedElectronID.clone( electronQuality = '70' )
process.vbtf11WP80 = simpleCutBasedElectronID.clone( electronQuality = '80' )
process.vbtf11WP85 = simpleCutBasedElectronID.clone( electronQuality = '85' )
process.vbtf11WP90 = simpleCutBasedElectronID.clone( electronQuality = '90' )
process.vbtf11WP95 = simpleCutBasedElectronID.clone( electronQuality = '95' )
process.eIdSequence += process.vbtf11WP60
process.eIdSequence += process.vbtf11WP70
process.eIdSequence += process.vbtf11WP80
process.eIdSequence += process.vbtf11WP85
process.eIdSequence += process.vbtf11WP90
process.eIdSequence += process.vbtf11WP95

from WWAnalysis.SkimStep.cutsInCategoriesHWWElectronIdentificationV04_cfi import *
process.cicVeryLooseHWW   = eidHWWVeryLoose.clone()
process.cicLooseHWW       = eidHWWLoose.clone()
process.cicMediumHWW      = eidHWWMedium.clone()
process.cicTightHWW       = eidHWWTight.clone()
process.cicSuperTightHWW  = eidHWWSuperTight.clone()
process.cicHyperTight1HWW = eidHWWHyperTight1.clone()
process.cicHyperTight2HWW = eidHWWHyperTight2.clone()
process.cicHyperTight3HWW = eidHWWHyperTight3.clone()
process.eIdSequence += process.cicVeryLooseHWW
process.eIdSequence += process.cicLooseHWW
process.eIdSequence += process.cicMediumHWW
process.eIdSequence += process.cicTightHWW
process.eIdSequence += process.cicSuperTightHWW 
process.eIdSequence += process.cicHyperTight1HWW 
process.eIdSequence += process.cicHyperTight2HWW 
process.eIdSequence += process.cicHyperTight3HWW 

process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
process.egammaIDLikelihood = process.eidLikelihoodExt.clone()
process.eIdSequence += process.egammaIDLikelihood 

from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
process.cicVeryLoose   = eidVeryLoose.clone()
process.cicLoose       = eidLoose.clone()
process.cicMedium      = eidMedium.clone()
process.cicTight       = eidTight.clone()
process.cicSuperTight  = eidSuperTight.clone()
process.cicHyperTight1 = eidHyperTight1.clone()
process.cicHyperTight2 = eidHyperTight2.clone()
process.cicHyperTight3 = eidHyperTight3.clone()
process.eIdSequence += process.cicVeryLoose 
process.eIdSequence += process.cicLoose 
process.eIdSequence += process.cicMedium 
process.eIdSequence += process.cicTight 
process.eIdSequence += process.cicSuperTight 
process.eIdSequence += process.cicHyperTight1 
process.eIdSequence += process.cicHyperTight2 
process.eIdSequence += process.cicHyperTight3 

from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi import *
process.cicVeryLooseMC   = eidVeryLooseMC.clone()
process.cicLooseMC       = eidLooseMC.clone()
process.cicMediumMC      = eidMediumMC.clone()
process.cicTightMC       = eidTightMC.clone()
process.cicSuperTightMC  = eidSuperTightMC.clone()
process.cicHyperTight1MC = eidHyperTight1MC.clone()
process.cicHyperTight2MC = eidHyperTight2MC.clone()
process.cicHyperTight3MC = eidHyperTight3MC.clone()
process.eIdSequence += process.cicVeryLooseMC 
process.eIdSequence += process.cicLooseMC 
process.eIdSequence += process.cicMediumMC 
process.eIdSequence += process.cicTightMC 
process.eIdSequence += process.cicSuperTightMC 
process.eIdSequence += process.cicHyperTight1MC 
process.eIdSequence += process.cicHyperTight2MC 
process.eIdSequence += process.cicHyperTight3MC 



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

if is41XRelease:
    process.preLeptonSequence += (process.offlinePrimaryVertices)


process.preLeptonSequence += ( 
    process.eIdSequence + 
    process.ak5PFJets + 
    process.kt6PFJets +
    process.kt6PFJetsForIso + 
    process.pfPileUp *
    process.pfNoPileUp * (
        process.ak5PFJetsNoPU +
        process.kt6PFJetsNoPU + 
        process.kt6PFJetsForIsoNoPU ) *
    process.valueMaps 
)


process.patDefaultSequence.remove( process.pfPileUp )
process.patDefaultSequence.remove( process.pfNoPileUp )



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
    cms.InputTag("convValueMapProd","dcot"),
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
    cms.InputTag("rhoMuNoPU"),
)
process.patMuons.isolationValues = cms.PSet(
#     pfNeutralHadrons = cms.InputTag("isoValMuonWithNeutralIso"),
#     pfChargedHadrons = cms.InputTag("isoValMuonWithChargedIso"),
#     pfPhotons = cms.InputTag("isoValMuonWithPhotonIso")
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
def addFastJetCorrection(process,label,seq="patDefaultSequence",thisRho="kt6PFJets"):
    if is41XRelease:
        print "==========================================================================================="
        print " _   _                   _____       _                      _       _     _   _      _ _ _ "
        print "| \ | |                 |  __ \     | |                    (_)     (_)   | | (_)    | | | |"
        print "|  \| | ___  _ __ ______| |  | | ___| |_ ___ _ __ _ __ ___  _ _ __  _ ___| |_ _  ___| | | |"
        print "| . ` |/ _ \| '_ \______| |  | |/ _ \ __/ _ \ '__| '_ ` _ \| | '_ \| / __| __| |/ __| | | |"
        print "| |\  | (_) | | | |     | |__| |  __/ ||  __/ |  | | | | | | | | | | \__ \ |_| | (__|_|_|_|"
        print "|_| \_|\___/|_| |_|     |_____/ \___|\__\___|_|  |_| |_| |_|_|_| |_|_|___/\__|_|\___(_|_|_)"
        print "==========================================================================================="

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
process.pfIsoNeutralHadrons = cms.EDFilter("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUp")
)
process.pfIsoChargedHadrons = process.pfIsoNeutralHadrons.clone ( pdgId = [211, -211, 321, -321, 999211, 2212, -2212] )
process.pfIsoPhotons        = process.pfIsoNeutralHadrons.clone ( pdgId = [22] )

# make the actual IsoDeposits
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
    process.pfIsoPhotons * (
        process.isoDepMuonWithChargedIso +
        process.isoDepMuonWithNeutralIso +
        process.isoDepMuonWithPhotonIso +
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
#muons
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

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#      fileName = cms.untracked.string('latinosYieldSkim.root'), 
    outputCommands =  cms.untracked.vstring(
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
#         'keep patJets_slimPatJetsTriggerMatchCalo_*_*',
#         'keep patJets_slimPatJetsTriggerMatchJPT_*_*',
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
#         'keep *_mergedSuperClusters_*_'+process.name_(),
        'keep *_kt6PF*_rho_'+process.name_(),
        # Debug info, usually commented out
        #'keep *_patTrigger_*_*',  
        #'keep *_l1extraParticles_*_*',  
    ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('patPath' )),
)
if is41XRelease:
    process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_'+process.name_())
else:
    process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_*')

process.prePatSequence  = cms.Sequence( process.preLeptonSequence + process.preElectronSequence + process.preMuonSequence)
process.postPatSequence = cms.Sequence( process.autreSeq + process.chargedMetSeq )



if not is41XRelease:
    # In order to use the offline vertices with BS constratint everywhere 
    massSearchReplaceAnyInputTag(process.preFilter,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.prePatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    massSearchReplaceAnyInputTag(process.postPatSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"),True)
    if doPF2PATAlso:
        massSearchReplaceAnyInputTag(process.patPF2PATSequencePFlow,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("offlinePrimaryVerticesWithBS"))



process.scrap      = cms.Path( process.noscraping ) 
process.outpath    = cms.EndPath(process.out)

if  doPF2PATAlso:
    process.patPath = cms.Path( process.preFilter + process.prePatSequence * process.patDefaultSequence * process.patPF2PATSequencePFlow * process.postPatSequence )
    process.fakPath = cms.Path( process.pre4Fakes + process.prePatSequence * process.patDefaultSequence * process.patPF2PATSequencePFlow * process.postPatSequence )
else:
    process.patPath = cms.Path( process.preFilter + process.prePatSequence * process.patDefaultSequence * process.postPatSequence)
    process.fakPath = cms.Path( process.pre4Fakes + process.prePatSequence * process.patDefaultSequence * process.postPatSequence )

# from WWAnalysis.SkimStep.skimTools import addIsolationInformation
# addIsolationInformation(process)

if not doFakeRates:
    process.schedule = cms.Schedule( process.patPath, process.scrap, process.outpath)
elif doFakeRates == 'also':
    process.out.SelectEvents.SelectEvents = [ 'patPath', 'fakPath' ]
    process.countOverlaps = cms.Path(process.preFilter + process.pre4Fakes)
    process.schedule = cms.Schedule( process.patPath, process.fakPath, process.countOverlaps, process.scrap, process.outpath)
elif doFakeRates == 'only':
    process.out.SelectEvents.SelectEvents = [ 'fakPath' ]
    process.schedule = cms.Schedule( process.fakPath, process.scrap, process.outpath)
    

