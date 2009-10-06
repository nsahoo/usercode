import FWCore.ParameterSet.Config as cms

process = cms.Process("ProjectCalafuria")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_31X_V5::All"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/users/mangano/Quarkonia/CMSSW_3_1_2/src/JPsi_312.root'
    )
)


# should include cfi file
process.myJPsiAnalysisMC = cms.EDProducer('JPsiAnalysisMC'
)



process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_*_*_ProjectCalafuria'
)
  
process.p = cms.Path(process.myJPsiAnalysisMC)
process.e = cms.EndPath(process.out)
