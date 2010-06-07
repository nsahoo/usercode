import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root'
]);



secFiles.extend( [ ]);


process = cms.Process("IpResolutions")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'SET_GLOBALTAG'

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("Validation.RecoTrack.TrackValidation_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = source

process.load('IpResoStudies.EDAnalyzers.vertexResponsesAndTrueResolutions_cfi')


process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("SET_OUTPUT"),
      closeFileFast = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

process.p = cms.Path(process.vertexResponsesAndTrueResolutions)
