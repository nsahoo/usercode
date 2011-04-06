import FWCore.ParameterSet.Config as cms

skimEventProducer = cms.EDProducer('SkimEventProducer',
   muTag       = cms.InputTag("wwMuons"),
   elTag       = cms.InputTag("wwElectrons"),
   softMuTag  = cms.InputTag("wwMuons4Veto"),
#   extraElTag  = cms.InputTag("wwElectrons"),
   jetTag      = cms.InputTag("slimPatJetsTriggerMatch"), 
   tagJetTag   = cms.InputTag("slimPatJetsTriggerMatch"), 
   pfMetTag    = cms.InputTag("pfMet"),
   tcMetTag    = cms.InputTag("tcMet"),
   vtxTag      = cms.InputTag("offlinePrimaryVertices"),
   sptTag      = cms.InputTag("vertexMapProd","sumPt"),
   spt2Tag     = cms.InputTag("vertexMapProd","sumPt2"),
   branchAlias = cms.string("ww"),
   hypoType = cms.string("WWELMU"),
)

