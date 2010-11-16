import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

from glob import glob
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        'file:/nfs/bluearc/group/skims/ww/oct29Skim/MuRun2010B/MuRun2010B_4_1_Ufm.root',
        'file:/nfs/bluearc/group/skims/ww/oct29Skim/MuRun2010B/MuRun2010B_5_1_6ue.root',
        'file:/nfs/bluearc/group/skims/ww/oct29Skim/MuRun2010B/MuRun2010B_6_1_zWG.root',
        'file:/nfs/bluearc/group/skims/ww/oct29Skim/MuRun2010B/MuRun2010B_7_1_8WM.root',
        'file:/nfs/bluearc/group/skims/ww/oct29Skim/MuRun2010B/MuRun2010B_8_1_4B1.root',
    )
)
#process.source.fileNames = [ 'file:WZEG_CentralSkim.root' ]
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')

process.load("WWAnalysis.AnalysisStep.tp_lepton_effs_cff")

process.runMuMu = cms.Path(process.tnpSimpleSequenceMuMu)
##process.runElEl = cms.Path(process.tnpSimpleSequenceElEl) ## fixme: off for now

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpLeptonID.root"))
