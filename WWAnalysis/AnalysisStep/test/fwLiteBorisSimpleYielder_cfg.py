import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        'file:/data/mangano/MC/Spring11/WW_414_SKIM_V04/ggToH150toWWto2L2Nu.step2.root'
        ), ## mandatory
    maxEvents   = cms.int32(-1),                            ## optional
    outputEvery = cms.uint32(100),                            ## optional
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('histo.fwlite.root'),  ## mandatory
)

process.myAnalyzer = cms.PSet(
    ## input specific for this analyzer
    skimEventTag = cms.InputTag('wwmumuIPCBL')
)
