import FWCore.ParameterSet.Config as cms

danieleBenedettiBDTProducer = cms.EDProducer("DanieleBenedettiBDTProducer",
    eleLabel = cms.InputTag("gsfElectrons"),
    vtxLabel = cms.InputTag("offlinePrimaryVertices"),
    bdtName  = cms.string("BDTSimpleCat"),
    bdtXML   = cms.string("WWAnalysis/Tools/data/DanieleBenetti_EleID_TMVA_BDTSimpleCat.weights.xml")
)

