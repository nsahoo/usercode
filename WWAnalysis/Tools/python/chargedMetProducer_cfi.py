import FWCore.ParameterSet.Config as cms

chargedMetProducer = cms.EDProducer('ChargedMetProducer',
    collectionTag = cms.untracked.InputTag("chargedParticleFlow"),
    vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
    dzCut = cms.untracked.double(0.1),
    minNeutralPt = cms.untracked.double(4.0),
    maxNeutralEta = cms.untracked.double(3.0),
)

