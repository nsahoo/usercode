import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/g/gpetrucc/scratch0/higgs/CMSSW_4_1_4/src/WWAnalysis/SkimStep/test/latinosYieldSkim.root',
        #'file:/data/gpetrucc/7TeV/hwww/DoubleElectron2011A_WW_414_SKIM_V05-v4_27_1_81X.root'
        #'/store/user/mwlebour/SingleElectron/R414_S1_V06_S2_V02_S3_V00_ID078_SingleElectron2011Av2/1fe9ae3d0313178a1fa1694ec77457b7/SingleElectron2011Av2_17_1_wAU.root',
        #'/store/user/mwlebour/SingleElectron/R414_S1_V06_S2_V02_S3_V00_ID078_SingleElectron2011Av2/1fe9ae3d0313178a1fa1694ec77457b7/SingleElectron2011Av2_38_1_tvj.root',
        #'/store/user/mwlebour/SingleMu/R414_S1_V06_S2_V02_S3_V00_ID079_SingleMuon2011Av2/1fe9ae3d0313178a1fa1694ec77457b7/SingleMuon2011Av2_49_1_XlR.root',
        #'/store/user/mwlebour/SingleMu/R414_S1_V06_S2_V02_S3_V00_ID079_SingleMuon2011Av2/1fe9ae3d0313178a1fa1694ec77457b7/SingleMuon2011Av2_75_1_6Av.root',
        '/store/user/mwlebour/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/R414_S1_V06_S2_V02_S3_V00_ID004_DYtoMuMu/29b3c44db3f5764c74aff4bbe6341c1e/DYtoMuMu_1_2_GZ9.root',
        '/store/user/mwlebour/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/R414_S1_V06_S2_V02_S3_V00_ID003_DYtoElEl/29b3c44db3f5764c74aff4bbe6341c1e/DYtoElEl_2_2_VhM.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')
#process.GlobalTag.globaltag = cms.string('START311_V2::All')

process.load("WWAnalysis.AnalysisStep.tp_lepton_effs_cff")

process.runMuMu = cms.Path(process.tnpSimpleSequenceMuMu)
process.runElEl = cms.Path(process.tnpSimpleSequenceElEl) ## fixme: off for now

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpLeptonID.root"))

