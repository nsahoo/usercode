import FWCore.ParameterSet.Config as cms

skimEventProducer = cms.EDProducer('SkimEventProducer',
   muTag       = cms.InputTag("wwMuons"),
   elTag       = cms.InputTag("wwElectrons"),
   extraMuTag  = cms.InputTag("wwMuons4Veto"),
   extraElTag  = cms.InputTag("wwElectrons"),
#   jetTag      = cms.InputTag("cleanPatJets"), ## for new skim
   jetTag      = cms.InputTag("patJets"),
   pfMetTag    = cms.InputTag("pfMet"),
   tcMetTag    = cms.InputTag("tcMet"),
   vtxTag      = cms.InputTag("offlinePrimaryVertices"),
   branchAlias = cms.string("ww"),
   hypoType = cms.string("WWELMU"),
)
