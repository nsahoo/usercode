import FWCore.ParameterSet.Config as cms

process = cms.Process("PatrickMET")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.GlobalTag.globaltag = 'START311_V2::All'

#Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/data/mangano/MC/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV_Spring11_AOD.root',
        'file:/data/mangano/MC/Spring11/MinBias.root',
    )
)

process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesDA_cfi")
process.offlinePrimaryVertices = process.offlinePrimaryVerticesDA.clone()
process.offlinePrimaryVertices.useBeamConstraint = cms.bool(True)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.Tmin = cms.double(4.)
process.offlinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = cms.double(0.01)


process.badAssPFCandidates = cms.EDProducer('PatrickMETProducer'
)

process.patrickMET = process.pfMet.clone(
    src = cms.InputTag("badAssPFCandidates")
)


process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('h160.root'),
    fileName = cms.untracked.string('newMet.root'),
    outputCommands =  cms.untracked.vstring(
        'drop *',
        #'keep *_*_*_PatrickMET',
        'keep recoPFMETs_*__*',
        )
)

  
process.p = cms.Path(#process.offlinePrimaryVertices*
                     process.badAssPFCandidates*
                     process.patrickMET)

process.e = cms.EndPath(process.out)
