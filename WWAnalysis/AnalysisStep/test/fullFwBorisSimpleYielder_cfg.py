import FWCore.ParameterSet.Config as cms

#process = cms.PSet()
process = cms.Process("Test")


### INPUT ###
#process.fwliteInput = cms.PSet(
#    fileNames   = cms.vstring(
#        'file:/data/mangano/MC/Spring11/WW_414_SKIM_V04/ggToH150toWWto2L2Nu.step2.root'
#        ), ## mandatory
#    maxEvents   = cms.int32(-1),                            ## optional
#    outputEvery = cms.uint32(100),                            ## optional
#)
 
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        'file:/data/mangano/MC/Spring11/WW_414_SKIM_V04/ggToH150toWWto2L2Nu.step2.root'
  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger = cms.Service("MessageLogger")
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
############## 
  


### OUTPUT ###
#process.fwliteOutput = cms.PSet(
#    fileName  = cms.string('histo.fw.root'),  ## mandatory
#)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('histo.fullfw.root')
)
###############


### MODULE TO RUN ###
#process.myAnalyzer = cms.PSet(
#    ## input specific for this analyzer
#    skimEventTag = cms.InputTag('wwmumuIPCBL')
#)


process.myAnalyzer = cms.EDAnalyzer("EDMWrappedBorisSimpleYielder",
  skimEventTag = cms.InputTag('wwmumuIPCBL'),                                           
)
process.p = cms.Path(process.myAnalyzer)
