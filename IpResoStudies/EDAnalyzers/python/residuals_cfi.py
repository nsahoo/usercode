import FWCore.ParameterSet.Config as cms

residuals = cms.EDAnalyzer("Residuals",
  TrackLabel = cms.InputTag("generalTracks"),     
  VertexLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),     

  TkMinPt = cms.double(0.5),
  TkMinNHits = cms.int32(11),

  VtxTracksSizeMin = cms.int32(10),
  VtxTracksSizeMax = cms.int32(16),
                          
)
