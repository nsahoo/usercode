import FWCore.ParameterSet.Config as cms

muPlusSelection1 = cms.EDFilter("MuonSelection1",
  src = cms.InputTag("muons"),
  filter = cms.bool(False),
  ptMin     = cms.double(20.0),                              
  absEtaMax = cms.double(2.4),                              
  charge    = cms.int32(1),                              
)
