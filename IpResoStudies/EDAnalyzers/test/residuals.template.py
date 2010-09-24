import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [  #before run introducing prescaling for jetTriggers (135445)
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/FE49901C-9E72-DF11-8384-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/F2BFED72-9E72-DF11-A02B-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/EA91F5D8-9D72-DF11-A90C-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E86BF28A-9E72-DF11-845E-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E673C3FA-9D72-DF11-97F1-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E2EBC24B-9D72-DF11-96B4-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E06FF96A-9D72-DF11-A9C3-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E05C5938-9C72-DF11-83A6-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/DAAFA989-9E72-DF11-B55E-0018F3D0961A.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D8311F50-9C72-DF11-9CF1-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D4D9AA9D-9C72-DF11-B547-001A92811732.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D2079EE2-9C72-DF11-849B-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D055CE02-9C72-DF11-BC06-0018F3D096F8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/C8E21D60-9C72-DF11-8544-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/C2D15BE8-9C72-DF11-A391-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/BED87A53-9C72-DF11-B2E9-003048678FA6.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/B4CE00C4-9B72-DF11-A7B1-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/B481581A-9F72-DF11-82A0-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/AAC54980-9F72-DF11-B3EE-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/9EFBB04F-9E72-DF11-AE61-002618943894.root',
]);



secFiles.extend( [ ]);


process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = source

#======================================
# HLT
#======================================
process.HLTMinBias = cms.EDFilter("HLTHighLevel",
#                               TriggerResultsTag = cms.InputTag("TriggerResults","","REDIGI36"),
                               TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#                               HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'), 
                               HLTPaths = cms.vstring('HLT_L1Jet6U'), 
                               eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                               andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                               throw = cms.bool(True)    # throw exception on unknown path names
                               )


process.load('IpResoStudies.EDAnalyzers.residuals_cfi')


process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("SET_OUTPUT"),
      closeFileFast = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) 

process.p = cms.Path(process.HLTMinBias
                     *process.residuals
                     )
