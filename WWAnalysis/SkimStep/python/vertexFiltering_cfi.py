import FWCore.ParameterSet.Config as cms

goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVerticesWithBS"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(False),
)
