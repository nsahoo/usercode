import FWCore.ParameterSet.Config as cms


muonSelection1 = cms.EDFilter("MuonSelection1",
  src = cms.InputTag("muons"),
#  filter = cms.bool(False),
  filter = cms.bool(True),
  ptMin     = cms.double(20.0),                              
  absEtaMax = cms.double(2.4),                              
  charge    = cms.int32(1),                              
)
