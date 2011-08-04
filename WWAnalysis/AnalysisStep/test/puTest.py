import FWCore.ParameterSet.Config as cms

# CMSSW Regular Stuff
process = cms.Process("Yield")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

#Options
process.GlobalTag.globaltag = 'START42_V13::All'

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('RMMEFN'))
process.source.fileNames = ['file:/nfs/bluearc/group/edm/hww/Summer11/RMME']

process.TFileService = cms.Service("TFileService", fileName = cms.string('RMME'))

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.load('WWAnalysis.SkimStep.vertexFiltering_cff')
from WWAnalysis.AnalysisStep.pileupReweighting_cfi import dataWeights, summer11
process.puHists = cms.EDAnalyzer("PUDistroMaker",
    lumi = cms.vdouble(dataWeights[:]),
    mc = cms.vdouble(summer11[:]),

)

process.p = cms.Path(process.goodPrimaryVertices + process.puHists)

