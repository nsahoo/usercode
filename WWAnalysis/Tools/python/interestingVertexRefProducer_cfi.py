import FWCore.ParameterSet.Config as cms

interestingVertexRefProducer = cms.EDProducer('InterestingVertexRefProducer',
    leptonTags = cms.untracked.VInputTag("muons","gsfElectrons"),
    vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
)

