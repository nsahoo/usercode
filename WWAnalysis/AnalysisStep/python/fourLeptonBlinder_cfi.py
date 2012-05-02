import FWCore.ParameterSet.Config as cms

fourLeptonBlinder = cms.EDFilter("FourLeptonBlinder",
   electrons = cms.InputTag("boostedElectrons"), 
   muons = cms.InputTag("boostedMuons"), 
   electronSelection = cms.InputTag("pt > 5"),
   muonSelection = cms.InputTag("pt > 3"),
   blindMC   = cms.bool(False),
   blindData = cms.bool(False),
   #firstBlindedDataRun = cms.uint32(190000),
   firstBlindedDataRun = cms.uint32(0),
)
