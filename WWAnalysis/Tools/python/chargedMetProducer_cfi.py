import FWCore.ParameterSet.Config as cms

chargedMetProducer = cms.EDProducer('ChargedMetProducer',
    collectionTag = cms.untracked.InputTag("chargedParticleFlow"),
    vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    dzCut = cms.untracked.double(0.2),
)

