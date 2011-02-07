import FWCore.ParameterSet.Config as cms

process = cms.Process("Iso")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")
process.GlobalTag.globaltag = 'START38_V12::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = ['file:Hww160.skim.root']
# from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('*.root') ]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.root'),
    outputCommands = cms.untracked.vstring(
#         'keep *_goodElectrons_*_*',
#         'keep *_goodMuons_*_*',
        'keep *_muTuple_*_*',
        'keep *_eleTuple_*_*',
#         'keep *_*_*_Iso',
        ),
)

process.load("RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi")
process.load("RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi")
process.muIsoDepositTk.IOPSet.inputMuonCollection = 'goodMuons'

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequencePAT_cff")
process.eleIsoDepositTk.src = 'goodElectrons'

# dZValues   = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 9999.0 ]
# d0Values   = [0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 9999.0 ]
# tkStrips   = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
# tkDeltaRs  = [0.000, 0.005, 0.010, 0.015, 0.020, 0.025]
# outerDRs   = [0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1]
# eneThresh  = [0.02 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
#               0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 
# etThresh   = [0.02 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
#               0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 ] 


process.isoSeq = cms.Sequence(
    process.eleIsoDepositTk +
    process.muIsoDepositTk 
)

process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")
process.boostedElectrons.electronTag = 'goodElectrons'
process.boostedMuons.muonTag = 'goodMuons'


process.load("RecoEgamma.EgammaIsolationAlgos.eleIsoFromDepsModules_cff")
process.boostedElectrons.deposits.append( process.eleIsoFromDepsTk.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('isoTk')
process.boostedElectrons.deposits[-1].deltaR = 0.3
process.boostedElectrons.deposits.append( process.eleIsoFromDepsEcalFromHitsByCrystal.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('isoEcal')
process.boostedElectrons.deposits[-1].deltaR = 0.4
process.boostedElectrons.deposits.append( process.eleIsoFromDepsHcalFromTowers.deposits[0].clone() )
process.boostedElectrons.deposits[-1].label = cms.string('isoHcal')
process.boostedElectrons.deposits[-1].deltaR = 0.3

# NoiseTow_EB = cms.double(0.04),
# NoiseTow_EE = cms.double(0.15),
# Threshold_E = cms.double(0.2),
# DR_Veto_E = cms.double(0.07),
# Noise_HB = cms.double(0.2),
# Noise_HE = cms.double(0.2),
# Threshold_H = cms.double(0.5),
# DR_Veto_H = cms.double(0.1),

# Diff_z = cms.double(0.2),
# Diff_r = cms.double(0.1),
# DR_Veto = cms.double(0.01),
# Pt_Min = cms.double(-1.0),


tempPSet =  cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("muIsoDepositTk"),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring('0.01','Threshold(0.0)'),
        skipDefaultVeto = cms.bool(True),
        label = cms.string('isoTk')
)

process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoTkM')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoEcalM')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
process.boostedMuons.deposits[-1].vetos = [
    'EcalBarrel:ThresholdFromTransverse(0.12)',
    'EcalEndcaps:ThresholdFromTransverse(0.45)',
    'Threshold(0.2)',
    '0.07'
]
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoHcalM')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')
process.boostedMuons.deposits[-1].vetos = [
    'ThresholdFromTransverse(0.6)',
    'Threshold(0.5)',
    '0.1'
]

process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoTk')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositTk')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoEcal')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','ecal')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False
process.boostedMuons.deposits.append( tempPSet.clone() )
process.boostedMuons.deposits[-1].label = cms.string('isoHcal')
process.boostedMuons.deposits[-1].src = cms.InputTag('muIsoDepositCalByAssociatorTowers','hcal')
process.boostedMuons.deposits[-1].vetos = []
process.boostedMuons.deposits[-1].skipDefaultVeto = False


process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.isoDepMuonWithNeutralHadrons = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("goodMuons"),
    MultipleDepositsFlag = cms.bool(False),
    trackType = cms.string('candidate'),
    ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
    )
)

process.pfAllPhotonsNoPU = process.pfAllPhotons.clone()
process.pfAllNeutralHadronsNoPU = process.pfAllNeutralHadrons.clone()
process.pfAllChargedHadronsNoPU = process.pfAllChargedHadrons.clone()
process.pfAllPhotons.src = "particleFlow"
process.pfAllNeutralHadrons.src = "particleFlow"
process.pfAllChargedHadrons.src = "particleFlow"

process.isoDepMuonWithPhotons = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithPhotons.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotons")
process.isoDepMuonWithChargedPFlow = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithChargedPFlow.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadrons")

process.isoDepMuonWithNeutralHadronsNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithNeutralHadronsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllNeutralHadronsNoPU")
process.isoDepMuonWithPhotonsNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithPhotonsNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllPhotonsNoPU")
process.isoDepMuonWithChargedPFlowNoPU = process.isoDepMuonWithNeutralHadrons.clone()
process.isoDepMuonWithChargedPFlowNoPU.ExtractorPSet.inputCandView = cms.InputTag("pfAllChargedHadronsNoPU")

process.pfSeq = cms.Sequence( (
        process.pfAllPhotons + 
        process.pfAllNeutralHadrons + 
        process.pfAllChargedHadrons  
    ) * (
        process.isoDepMuonWithNeutralHadrons +
        process.isoDepMuonWithPhotons + 
        process.isoDepMuonWithChargedPFlow
    ) + 
    process.pfPileUp * 
    process.pfNoPileUp * ( 
        process.pfAllPhotonsNoPU + 
        process.pfAllNeutralHadronsNoPU + 
        process.pfAllChargedHadronsNoPU 
    ) * (
        process.isoDepMuonWithNeutralHadronsNoPU +
        process.isoDepMuonWithPhotonsNoPU + 
        process.isoDepMuonWithChargedPFlowNoPU
    )
)

process.eleTuple = cms.EDAnalyzer("BetterElectronTupleProducer", 
    src = cms.untracked.InputTag("boostedElectrons"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("eid80"),              quantity = cms.untracked.string("electronID('eidVBTFRel80')")),
        cms.PSet( tag = cms.untracked.string("eid95"),              quantity = cms.untracked.string("electronID('eidVBTFRel95')")),
        cms.PSet( tag = cms.untracked.string("liklihoodID"),        quantity = cms.untracked.string("electronID('liklihoodID')")),
        cms.PSet( tag = cms.untracked.string("isoTk"),              quantity = cms.untracked.string("userFloat('isoTk')")),
        cms.PSet( tag = cms.untracked.string("pt"),                 quantity = cms.untracked.string("pt")),
        cms.PSet( tag = cms.untracked.string("eta"),                quantity = cms.untracked.string("eta")),
        cms.PSet( tag = cms.untracked.string("dxyPV"),              quantity = cms.untracked.string("userFloat('dxyPV')")),
        cms.PSet( tag = cms.untracked.string("dzPV"),               quantity = cms.untracked.string("userFloat('dzPV')")),
        cms.PSet( tag = cms.untracked.string("dcot"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dcot')")),
        cms.PSet( tag = cms.untracked.string("dist"),               quantity = cms.untracked.string("userFloat('convValueMapProd:dist')")),
    ), 
    weight = cms.untracked.double(1)
)  

process.muTuple = cms.EDAnalyzer("BetterMuonTupleProducer", 
    src = cms.untracked.InputTag("boostedMuons"),
    variables = cms.untracked.VPSet(
        cms.PSet( tag = cms.untracked.string("pt"),                 quantity = cms.untracked.string("pt")),
        cms.PSet( tag = cms.untracked.string("eta"),                quantity = cms.untracked.string("eta")),
        cms.PSet( tag = cms.untracked.string("isoTk"),              quantity = cms.untracked.string("userFloat('isoTk')")),
        cms.PSet( tag = cms.untracked.string("dxyPV"),              quantity = cms.untracked.string("userFloat('dxyPV')")),
        cms.PSet( tag = cms.untracked.string("dzPV"),               quantity = cms.untracked.string("userFloat('dzPV')")),
    ), 
    weight = cms.untracked.double(1)
)  


process.TFileService = cms.Service("TFileService",fileName = cms.string("hwwTuple.root"))

process.p = cms.Path( (
        process.isoSeq + 
        process.pfSeq
    ) * (
        ( process.boostedElectrons * process.eleTuple ) + 
        ( process.boostedMuons * process.muTuple ) 
    )
)

# process.e = cms.EndPath(process.out)


