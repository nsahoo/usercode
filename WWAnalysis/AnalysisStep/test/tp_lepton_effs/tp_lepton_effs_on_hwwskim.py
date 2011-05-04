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
        'file:/afs/cern.ch/user/g/gpetrucc/scratch0/higgs/CMSSW_4_1_4/src/WWAnalysis/SkimStep/test/latinosYieldSkim.DY_EE.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')
process.GlobalTag.globaltag = cms.string('START311_V2::All')

process.load("WWAnalysis.AnalysisStep.tp_lepton_effs_cff")

##process.runMuMu = cms.Path(process.tnpSimpleSequenceMuMu)
process.runElEl = cms.Path(process.tnpSimpleSequenceElEl) ## fixme: off for now

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpLeptonID_mc.root"))

MC=True
if MC:
    process.tagElectrons.cut = process.tagElectrons.cut.value().replace("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT","HLT_Ele22_SW_TighterCaloIdIsol_L1R");
    process.tpTreeElEl.flags.HLT_Ele8_CaloIdL_CaloIsoVL     = "!triggerObjectMatchesByPath('HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*').empty"
    process.tpTreeElEl.flags.HLT_Ele17_CaloIdL_CaloIsoVL    = "!triggerObjectMatchesByFilter('hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEPixelMatchFilter').empty"
    process.tpTreeElEl.tagFlags.HLT_Ele17_CaloIdL_CaloIsoVL = "!triggerObjectMatchesByFilter('hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEPixelMatchFilter').empty"

