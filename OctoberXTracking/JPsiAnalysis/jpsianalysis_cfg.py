import FWCore.ParameterSet.Config as cms

process = cms.Process("OctoberXTracking")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MC_31X_V5::All"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/users/mangano/OctoberExercise/CMSSW_3_1_4/src/JPsiMuMu.AOD.root'
    )
)


# should include cfi file
process.myJPsiAnalysis = cms.EDProducer('JPsiAnalysis'
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep edmTriggerResults_TriggerResults_*_*',
    #'keep edmTriggerResults_TriggerResults__HLT8E29',
    #'keep edmTriggerResults_TriggerResults__HLT',
    'keep *_*_*_OctoberXTracking'
)
  
process.p = cms.Path(process.myJPsiAnalysis)
process.e = cms.EndPath(process.out)
