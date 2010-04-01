import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGE")

### standard includes
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

### global tag
#process.GlobalTag.globaltag = "GR09_R_35X_V3::All"


### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#put your input files here
   )
)

### avoid complains if the events are not sorted
process.source.noEventSort = cms.untracked.bool(True) 

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


### path
#process.p = cms.Path(
#)


### output
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('merged.root'),
    outputCommands = cms.untracked.vstring('keep *',
#        'keep *_whatever_*_*'
#        'drop *_whatever_*_*'
    ),
    ## Uncomment to activate skimming!
    #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
)
process.e = cms.EndPath(process.out)


