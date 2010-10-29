import FWCore.ParameterSet.Config as cms

convValueMapProd = cms.EDProducer('ConvValueMapProd',
    gsfLabel = cms.untracked.InputTag("gsfElectrons"),
    tkLabel = cms.untracked.InputTag("generalTracks")
)
