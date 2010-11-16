import FWCore.ParameterSet.Config as cms

skimEventProducer = cms.EDProducer('SkimEventProducer',
   muTag       = cms.InputTag("wwMuons"),
   elTag       = cms.InputTag("wwElectrons"),
   extraMuTag  = cms.InputTag("wwMuons4Veto"),
   extraElTag  = cms.InputTag("wwElectrons"),
   jetTag      = cms.InputTag("cleanPatJetsPF"), 
   pfMetTag    = cms.InputTag("pfMet"),
   tagJetTag   = cms.InputTag("cleanPatJets"), 
   tcMetTag    = cms.InputTag("tcMet"),
   vtxTag      = cms.InputTag("offlinePrimaryVertices"),
   sptTag      = cms.InputTag("vertexMapProd","sumPt"),
   spt2Tag     = cms.InputTag("vertexMapProd","sumPt2"),
   branchAlias = cms.string("ww"),
   hypoType = cms.string("WWELMU"),
)
