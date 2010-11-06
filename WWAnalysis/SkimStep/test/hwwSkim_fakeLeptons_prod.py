from hwwSkim import *

process.jetTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patJets" ),
    matched = cms.InputTag( "patTrigger" ),
    andOr          = cms.bool( False ),
    filterIdsEnum  = cms.vstring( '*' ),
    filterIds      = cms.vint32( 0 ),
    filterLabels   = cms.vstring( '*' ),
    pathNames      = cms.vstring( '*' ),
    collectionTags = cms.vstring( 'hltMCJetCorJetIcone5HF07' ),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( False )
)
process.jetTriggerMatchHLTNoPU = process.jetTriggerMatchHLT.clone(src = "patJetsNoPU")
process.patJetsWithTrigger = cms.EDProducer( "PATTriggerMatchJetEmbedder",
    src     = cms.InputTag(  "patJets" ),
    matches = cms.VInputTag("jetTriggerMatchHLT")
)
process.patJetsWithTriggerNoPU = cms.EDProducer( "PATTriggerMatchJetEmbedder",
    src     = cms.InputTag(  "patJetsNoPU" ),
    matches = cms.VInputTag("jetTriggerMatchHLTNoPU")
)
process.cleanPatJets.src = "patJetsWithTrigger"
process.cleanPatJetsNoPU.src = "patJetsWithTriggerNoPU"
process.autreSeq.replace(process.cleanPatJets,     process.jetTriggerMatchHLT *     process.patJetsWithTrigger *     process.cleanPatJets)
process.autreSeq.replace(process.cleanPatJetsNoPU, process.jetTriggerMatchHLTNoPU * process.patJetsWithTriggerNoPU * process.cleanPatJetsNoPU)

HLT_JetU = " || ".join([ " || ".join(["!triggerObjectMatchesByPath('HLT_Jet%dU%s').empty()" % (I,P) for P in "", "_v2", "_v3"]) for I in 15, 30, 50, 70])
process.patJetsPassingTriggerNoPU = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("cleanPatJetsNoPU"),
    cut = cms.string(HLT_JetU),
)

process.muonsJet = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("patJetsPassingTriggerNoPU patMuonsWithTrigger"),
    cut = cms.string("deltaR(daughter(0).eta,daughter(0).phi,daughter(1).eta,daughter(1).phi) > 1.0 && daughter(0).pt > 10 && daughter(1).pt > 10"),
    checkCharge = cms.bool(False),
    checkOverlap = cms.bool(False),
)
process.electronsJet    = process.muonsJet.clone(decay = cms.string("patJetsPassingTriggerNoPU patElectronsWithTrigger"))

process.muonsJetFilter        = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("muonsJet"),        minNumber = cms.uint32(1))
process.electronsJetFilter    = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("electronsJet"),    minNumber = cms.uint32(1))

process.electronsUnbiasTrigger = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("patElectronsWithTrigger"),
    cut = cms.string(("pt > 10 && "+
                     "(!triggerObjectMatchesByPath('HLT_Photon10_Cleaned_L1R').empty() ||"+
                     " !triggerObjectMatchesByPath('HLT_L1SingleEG8').empty() ||"+
                     " !triggerObjectMatchesByPath('HLT_Ele10_SW_L1R').empty() ||"+
                     " !triggerObjectMatchesByPath('HLT_Ele17_SW_L1R').empty() ||"+
                     " !triggerObjectMatchesByPath('HLT_Ele22_SW_L1R').empty())")),
    filter = cms.bool(True),
)

process.muonsUnbiasTrigger = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 10 && (!triggerObjectMatchesByPath('HLT_Mu9').empty() || !triggerObjectMatchesByPath('HLT_Mu15_v1').empty()) "),
    filter = cms.bool(True),
)

## Quick filters to avoid running most of the sequence if there are no leptons of interest
process.ffMu  = cms.EDFilter("MuonSelector",        src = cms.InputTag("muons"),        cut = cms.string("isGlobalMuon && pt > 10"), filter = cms.bool(True))
process.ffEl  = cms.EDFilter("GsfElectronSelector", src = cms.InputTag("gsfElectrons"), cut = cms.string("pt > 10"),                 filter = cms.bool(True))

del process.sched
del process.wwAllPath

process.wwAllSeq.remove(process.countWWAll)
process.recoSeq = cms.Sequence(
    ( process.elSeq + process.muSeq ) *  
    process.wwAllSeq + 
    process.autreSeq +
    process.patJetsPassingTriggerNoPU +
    process.muonsJet +
    process.electronsJet 
)

if True:
    process.GlobalTag.globaltag = 'GR_R_38X_V14::All'
    removeMCMatching(process, ['Electrons','Muons','Jets'])
    process.preElectronSequence.remove(process.prunedGen) 
    process.preElectronSequence.remove(process.electronMatch) 
    process.preMuonSequence.remove(process.prunedGen)
    process.preMuonSequence.remove(process.muonMatch)
    process.makePatJets.remove(process.patJetPartonMatchNoPU)
    process.makePatJets.remove(process.patJetPartonMatchJPT)
    process.makePatJets.remove(process.patJetGenJetMatchNoPU)
    process.makePatJets.remove(process.patJetGenJetMatchJPT)
    process.patJetsNoPU.addGenJetMatch = False
    process.patJetsNoPU.addGenJetMatch = False
    process.patJetsNoPU.addGenJetMatch = False
    process.patJetsNoPU.addGenPartonMatch = False
    process.patJetsNoPU.getJetMCFlavour = False
    process.patJetsJPT.addGenJetMatch = False
    process.patJetsJPT.addGenPartonMatch = False
    process.patJetsJPT.getJetMCFlavour = False

process.out.SelectEvents.SelectEvents = []
for X,P in [ ("muonsJetFilter",         1), 
             ("electronsJetFilter",     1), 
             ("muonsUnbiasTrigger",     2),  ## <- NOTE the prescale
             ("electronsUnbiasTrigger", 1)] :
    setattr(process, X+"Prescale", cms.EDFilter("Prescaler",  prescaleFactor = cms.int32(P), prescaleOffset = cms.int32(0)))
    preFilter = process.ffMu if "muon" in X else process.ffEl
    setattr(process, "Skim_"+X,    cms.Path(getattr(process,X+"Prescale") + preFilter + process.recoSeq + getattr(process,X)))
    process.out.SelectEvents.SelectEvents += [ "Skim_"+X ] 
