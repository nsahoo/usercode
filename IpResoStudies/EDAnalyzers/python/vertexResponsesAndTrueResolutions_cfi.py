import FWCore.ParameterSet.Config as cms

vertexResponsesAndTrueResolutions = cms.EDAnalyzer("VertexResponsesAndTrueResolutions",

  # Selection of Tracks                             
  TrackLabel = cms.InputTag("generalTracks"),     
  TkMinPt = cms.double(0.3),
  TkMinXLayers = cms.int32(9),
  TkMaxMissedOuterLayers = cms.int32(4),
  TkMaxMissedInnerLayers = cms.int32(0),
                           

  # Selection of Vertices                         
  #VertexLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),       
  VertexLabel = cms.InputTag("offlinePrimaryVertices"),       
  VtxTracksSizeMin = cms.int32(2),
  VtxTracksSizeMax = cms.int32(300),
  VtxErrorXMin = cms.double(0.0027),
  VtxErrorXMax = cms.double(0.0033),
  VtxErrorYMin = cms.double(0.0027),
  VtxErrorYMax = cms.double(0.0033),
  VtxErrorZMin = cms.double(0.0030),
  VtxErrorZMax = cms.double(0.0045),

                           
)
