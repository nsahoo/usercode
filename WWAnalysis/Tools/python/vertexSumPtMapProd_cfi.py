import FWCore.ParameterSet.Config as cms

vertexMapProd = cms.EDProducer('VertexSumPtMapProd',
    vertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
)
