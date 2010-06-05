import FWCore.ParameterSet.Config as cms

residuals = cms.EDAnalyzer("Residuals",

  # Selection of Tracks                             
  TrackLabel = cms.InputTag("generalTracks"),     
  TkMinPt = cms.double(0.3),
  TkMinXLayers = cms.int32(7),
  TkMaxMissedOuterLayers = cms.int32(4),
  TkMaxMissedInnerLayers = cms.int32(0),
                           

  # Selection of Vertices                         
  #VertexLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),       
  VertexLabel = cms.InputTag("offlinePrimaryVertices"),       
  VtxTracksSizeMin = cms.int32(2),
  VtxTracksSizeMax = cms.int32(300),
  VtxErrorXMin = cms.double(0.0018),
  VtxErrorXMax = cms.double(0.0024),
  VtxErrorYMin = cms.double(0.0018),
  VtxErrorYMax = cms.double(0.0024),
  VtxErrorZMin = cms.double(0.0022),
  VtxErrorZMax = cms.double(0.0029),                           
)
