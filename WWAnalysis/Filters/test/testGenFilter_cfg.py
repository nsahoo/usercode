import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [  
'file:input.root',
]);



secFiles.extend( [ ]);


process = cms.Process("GenFilter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START38_V12::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = source

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 


process.load('WWAnalysis.Filters.genFilterWW2L2Nu_cfi')



process.p1 = cms.Path(process.genFilterWW2L2Nu)

# ---- endPath ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("output.root"),
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p1"))
)
process.end = cms.EndPath(process.out)


