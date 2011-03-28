import FWCore.ParameterSet.Config as cms

boostedMuons = cms.EDProducer("PatMuonBooster",
  muonTag = cms.untracked.InputTag("cleanPatMuonsTriggerMatch"), 
  trackTag = cms.untracked.InputTag("generalTracks"),       
  vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
  deposits = cms.VPSet()
)

boostedElectrons = cms.EDProducer("PatElectronBooster",
    electronTag = cms.untracked.InputTag("cleanPatElectronsTriggerMatch"),
    trackTag = cms.untracked.InputTag("generalTracks"),       
    vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    electronIDSources = cms.PSet(),
    deposits = cms.VPSet(),
)
