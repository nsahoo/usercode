import FWCore.ParameterSet.Config as cms

process = cms.Process("WWEvents")

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('TrackingTools.Configuration.TrackingTools_cff')

# process.GlobalTag.globaltag = 'START311_V2::All'
process.GlobalTag.globaltag = 'RMMEGlobalTag'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
    inputCommands = cms.untracked.vstring( "keep *" )
)
##############################################

#from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
#process.source.fileNames += [ 'file:%s'%x for x in glob('/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.root') ]

# process.source.fileNames = ['file:hwwSkim.root']

# import commands
# myDir = "/castor/cern.ch/user/m/mwlebour/hww/trees/WW_414_SKIM_V00/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V00_ID026_WJetsToLNuMad_Xfer/2bd6c1b7661d2d45dc2e9f49eb14f680/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")


# 0
process.wwelmu0 = process.skimEventProducer.clone(hypoType = "WWELMU", muTag="wwMuMatch", elTag = "wwEleMatch" )
process.wwmuel0 = process.skimEventProducer.clone(hypoType = "WWMUEL", muTag="wwMuMatch", elTag = "wwEleMatch" )
process.wwelel0 = process.skimEventProducer.clone(hypoType = "WWELEL", muTag="wwMuMatch", elTag = "wwEleMatch" )
process.wwmumu0 = process.skimEventProducer.clone(hypoType = "WWMUMU", muTag="wwMuMatch", elTag = "wwEleMatch" )

# LHL
process.wwelmuIDLHL   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHL" )
process.wwmuelIDLHL   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHL" )
process.wwelelIDLHL   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHL" )
process.wwmumuIDLHL   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHL" )

process.wwelmuISOLHL  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHL" )
process.wwmuelISOLHL  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHL" )
process.wwelelISOLHL  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHL" )
process.wwmumuISOLHL  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHL" )

process.wwelmuCONVLHL = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHL" )
process.wwmuelCONVLHL = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHL" )
process.wwelelCONVLHL = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHL" )
process.wwmumuCONVLHL = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHL" )

process.wwelmuIPLHL   = process.wwelmu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHL" )
process.wwmuelIPLHL   = process.wwmuel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHL" )
process.wwelelIPLHL   = process.wwelel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHL" )
process.wwmumuIPLHL   = process.wwmumu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHL" )

# LHT
process.wwelmuIDLHT   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmuelIDLHT   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwelelIDLHT   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )
process.wwmumuIDLHT   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDLHT" )

process.wwelmuISOLHT  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmuelISOLHT  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwelelISOLHT  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )
process.wwmumuISOLHT  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOLHT" )

process.wwelmuCONVLHT = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmuelCONVLHT = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwelelCONVLHT = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )
process.wwmumuCONVLHT = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVLHT" )

process.wwelmuIPLHT   = process.wwelmu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHT" )
process.wwmuelIPLHT   = process.wwmuel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHT" )
process.wwelelIPLHT   = process.wwelel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHT" )
process.wwmumuIPLHT   = process.wwmumu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPLHT" )

# CBL
process.wwelmuIDCBL   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBL" )
process.wwmuelIDCBL   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBL" )
process.wwelelIDCBL   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBL" )
process.wwmumuIDCBL   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBL" )

process.wwelmuISOCBL  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBL" )
process.wwmuelISOCBL  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBL" )
process.wwelelISOCBL  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBL" )
process.wwmumuISOCBL  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBL" )

process.wwelmuCONVCBL = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBL" )
process.wwmuelCONVCBL = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBL" )
process.wwelelCONVCBL = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBL" )
process.wwmumuCONVCBL = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBL" )

process.wwelmuIPCBL   = process.wwelmu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBL" )
process.wwmuelIPCBL   = process.wwmuel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBL" )
process.wwelelIPCBL   = process.wwelel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBL" )
process.wwmumuIPCBL   = process.wwmumu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBL" )

# CBT
process.wwelmuIDCBT   = process.wwelmu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBT" )
process.wwmuelIDCBT   = process.wwmuel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBT" )
process.wwelelIDCBT   = process.wwelel0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBT" )
process.wwmumuIDCBT   = process.wwmumu0.clone(muTag = "wwMuonsID",  elTag = "wwEleIDCBT" )

process.wwelmuISOCBT  = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBT" )
process.wwmuelISOCBT  = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBT" )
process.wwelelISOCBT  = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBT" )
process.wwmumuISOCBT  = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleISOCBT" )

process.wwelmuCONVCBT = process.wwelmu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBT" )
process.wwmuelCONVCBT = process.wwmuel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBT" )
process.wwelelCONVCBT = process.wwelel0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBT" )
process.wwmumuCONVCBT = process.wwmumu0.clone(muTag = "wwMuonsISO", elTag = "wwEleCONVCBT" )

process.wwelmuIPCBT   = process.wwelmu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBT" )
process.wwmuelIPCBT   = process.wwmuel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBT" )
process.wwelelIPCBT   = process.wwelel0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBT" )
process.wwmumuIPCBT   = process.wwmumu0.clone(muTag = "wwMuonsIP",  elTag = "wwEleIPCBT" )

##############################################

#0
process.skimElMu0 = cms.EDFilter("SkimEventSelector",
   src = cms.InputTag("wwelmu0"),
   filter = cms.bool(True),
   cut = cms.string("nLep >=2 "),                                   
)
process.skimMuEl0 = process.skimElMu0.clone( src = "wwmuel0" )
process.skimElEl0 = process.skimElMu0.clone( src = "wwelel0" )
process.skimMuMu0 = process.skimElMu0.clone( src = "wwmumu0" )

# LHL
process.skimElMuIDLHL = process.skimElMu0.clone( src = "wwelmuIDLHL" )
process.skimMuElIDLHL = process.skimElMu0.clone( src = "wwmuelIDLHL" )
process.skimElElIDLHL = process.skimElMu0.clone( src = "wwelelIDLHL" )
process.skimMuMuIDLHL = process.skimElMu0.clone( src = "wwmumuIDLHL" )

process.skimElMuISOLHL = process.skimElMu0.clone( src = "wwelmuISOLHL" )
process.skimMuElISOLHL = process.skimElMu0.clone( src = "wwmuelISOLHL" )
process.skimElElISOLHL = process.skimElMu0.clone( src = "wwelelISOLHL" )
process.skimMuMuISOLHL = process.skimElMu0.clone( src = "wwmumuISOLHL" )

process.skimElMuCONVLHL = process.skimElMu0.clone( src = "wwelmuCONVLHL" )
process.skimMuElCONVLHL = process.skimElMu0.clone( src = "wwmuelCONVLHL" )
process.skimElElCONVLHL = process.skimElMu0.clone( src = "wwelelCONVLHL" )
process.skimMuMuCONVLHL = process.skimElMu0.clone( src = "wwmumuCONVLHL" )

process.skimElMuIPLHL = process.skimElMu0.clone( src = "wwelmuIPLHL" )
process.skimMuElIPLHL = process.skimElMu0.clone( src = "wwmuelIPLHL" )
process.skimElElIPLHL = process.skimElMu0.clone( src = "wwelelIPLHL" )
process.skimMuMuIPLHL = process.skimElMu0.clone( src = "wwmumuIPLHL" )

# LHT
process.skimElMuIDLHT = process.skimElMu0.clone( src = "wwelmuIDLHL" )
process.skimMuElIDLHT = process.skimElMu0.clone( src = "wwmuelIDLHL" )
process.skimElElIDLHT = process.skimElMu0.clone( src = "wwelelIDLHL" )
process.skimMuMuIDLHT = process.skimElMu0.clone( src = "wwmumuIDLHL" )

process.skimElMuISOLHT = process.skimElMu0.clone( src = "wwelmuISOLHL" )
process.skimMuElISOLHT = process.skimElMu0.clone( src = "wwmuelISOLHL" )
process.skimElElISOLHT = process.skimElMu0.clone( src = "wwelelISOLHL" )
process.skimMuMuISOLHT = process.skimElMu0.clone( src = "wwmumuISOLHL" )

process.skimElMuCONVLHT = process.skimElMu0.clone( src = "wwelmuCONVLHL" )
process.skimMuElCONVLHT = process.skimElMu0.clone( src = "wwmuelCONVLHL" )
process.skimElElCONVLHT = process.skimElMu0.clone( src = "wwelelCONVLHL" )
process.skimMuMuCONVLHT = process.skimElMu0.clone( src = "wwmumuCONVLHL" )

process.skimElMuIPLHT = process.skimElMu0.clone( src = "wwelmuIPLHL" )
process.skimMuElIPLHT = process.skimElMu0.clone( src = "wwmuelIPLHL" )
process.skimElElIPLHT = process.skimElMu0.clone( src = "wwelelIPLHL" )
process.skimMuMuIPLHT = process.skimElMu0.clone( src = "wwmumuIPLHL" )

# CBL
process.skimElMuIDCBL = process.skimElMu0.clone( src = "wwelmuIDCBL" )
process.skimMuElIDCBL = process.skimElMu0.clone( src = "wwmuelIDCBL" )
process.skimElElIDCBL = process.skimElMu0.clone( src = "wwelelIDCBL" )
process.skimMuMuIDCBL = process.skimElMu0.clone( src = "wwmumuIDCBL" )

process.skimElMuISOCBL = process.skimElMu0.clone( src = "wwelmuISOCBL" )
process.skimMuElISOCBL = process.skimElMu0.clone( src = "wwmuelISOCBL" )
process.skimElElISOCBL = process.skimElMu0.clone( src = "wwelelISOCBL" )
process.skimMuMuISOCBL = process.skimElMu0.clone( src = "wwmumuISOCBL" )

process.skimElMuCONVCBL = process.skimElMu0.clone( src = "wwelmuCONVCBL" )
process.skimMuElCONVCBL = process.skimElMu0.clone( src = "wwmuelCONVCBL" )
process.skimElElCONVCBL = process.skimElMu0.clone( src = "wwelelCONVCBL" )
process.skimMuMuCONVCBL = process.skimElMu0.clone( src = "wwmumuCONVCBL" )

process.skimElMuIPCBL = process.skimElMu0.clone( src = "wwelmuIPCBL" )
process.skimMuElIPCBL = process.skimElMu0.clone( src = "wwmuelIPCBL" )
process.skimElElIPCBL = process.skimElMu0.clone( src = "wwelelIPCBL" )
process.skimMuMuIPCBL = process.skimElMu0.clone( src = "wwmumuIPCBL" )

# CBT
process.skimElMuIDCBT = process.skimElMu0.clone( src = "wwelmuIDCBT" )
process.skimMuElIDCBT = process.skimElMu0.clone( src = "wwmuelIDCBT" )
process.skimElElIDCBT = process.skimElMu0.clone( src = "wwelelIDCBT" )
process.skimMuMuIDCBT = process.skimElMu0.clone( src = "wwmumuIDCBT" )

process.skimElMuISOCBT = process.skimElMu0.clone( src = "wwelmuISOCBT" )
process.skimMuElISOCBT = process.skimElMu0.clone( src = "wwmuelISOCBT" )
process.skimElElISOCBT = process.skimElMu0.clone( src = "wwelelISOCBT" )
process.skimMuMuISOCBT = process.skimElMu0.clone( src = "wwmumuISOCBT" )

process.skimElMuCONVCBT = process.skimElMu0.clone( src = "wwelmuCONVCBT" )
process.skimMuElCONVCBT = process.skimElMu0.clone( src = "wwmuelCONVCBT" )
process.skimElElCONVCBT = process.skimElMu0.clone( src = "wwelelCONVCBT" )
process.skimMuMuCONVCBT = process.skimElMu0.clone( src = "wwmumuCONVCBT" )

process.skimElMuIPCBT = process.skimElMu0.clone( src = "wwelmuIPCBT" )
process.skimMuElIPCBT = process.skimElMu0.clone( src = "wwmuelIPCBT" )
process.skimElElIPCBT = process.skimElMu0.clone( src = "wwelelIPCBT" )
process.skimMuMuIPCBT = process.skimElMu0.clone( src = "wwmumuIPCBT" )


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep edmTriggerResults_TriggerResults_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        # Misc
        'keep *_addPileupInfo_*_*',
        # skimEvent
        'keep *_wwelmu*_*_*',
        'keep *_wwmuel*_*_*',
        'keep *_wwelel*_*_*',
        'keep *_wwmumu*_*_*',
        ),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMu*','selMuEl*','selElEl*','selMuMu*' ))
)

# 0
process.selElMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 * process.skimElMu0)
process.selMuEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 * process.skimMuEl0)
process.selElEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 * process.skimElEl0)
process.selMuMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 * process.skimMuMu0)

# LHL
process.selElMuIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDLHL * process.skimElMuIDLHL)
process.selMuElIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDLHL * process.skimMuElIDLHL)
process.selElElIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDLHL * process.skimElElIDLHL)
process.selMuMuIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDLHL * process.skimMuMuIDLHL)

process.selElMuISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOLHL * process.skimElMuISOLHL)
process.selMuElISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOLHL * process.skimMuElISOLHL)
process.selElElISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOLHL * process.skimElElISOLHL)
process.selMuMuISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOLHL * process.skimMuMuISOLHL)

process.selElMuCONVLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuCONVLHL * process.skimElMuCONVLHL)
process.selMuElCONVLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelCONVLHL * process.skimMuElCONVLHL)
process.selElElCONVLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelCONVLHL * process.skimElElCONVLHL)
process.selMuMuCONVLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuCONVLHL * process.skimMuMuCONVLHL)

process.selElMuIPLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIPLHL * process.skimElMuIPLHL)
process.selMuElIPLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIPLHL * process.skimMuElIPLHL)
process.selElElIPLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIPLHL * process.skimElElIPLHL)
process.selMuMuIPLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIPLHL * process.skimMuMuIPLHL)

# LHT
process.selElMuIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDLHT * process.skimElMuIDLHT)
process.selMuElIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDLHT * process.skimMuElIDLHT)
process.selElElIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDLHT * process.skimElElIDLHT)
process.selMuMuIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDLHT * process.skimMuMuIDLHT)

process.selElMuISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOLHT * process.skimElMuISOLHT)
process.selMuElISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOLHT * process.skimMuElISOLHT)
process.selElElISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOLHT * process.skimElElISOLHT)
process.selMuMuISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOLHT * process.skimMuMuISOLHT)

process.selElMuCONVLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuCONVLHT * process.skimElMuCONVLHT)
process.selMuElCONVLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelCONVLHT * process.skimMuElCONVLHT)
process.selElElCONVLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelCONVLHT * process.skimElElCONVLHT)
process.selMuMuCONVLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuCONVLHT * process.skimMuMuCONVLHT)

process.selElMuIPLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIPLHT * process.skimElMuIPLHT)
process.selMuElIPLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIPLHT * process.skimMuElIPLHT)
process.selElElIPLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIPLHT * process.skimElElIPLHT)
process.selMuMuIPLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIPLHT * process.skimMuMuIPLHT)

# CBL
process.selElMuIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDCBL * process.skimElMuIDCBL)
process.selMuElIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDCBL * process.skimMuElIDCBL)
process.selElElIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDCBL * process.skimElElIDCBL)
process.selMuMuIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDCBL * process.skimMuMuIDCBL)

process.selElMuISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOCBL * process.skimElMuISOCBL)
process.selMuElISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOCBL * process.skimMuElISOCBL)
process.selElElISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOCBL * process.skimElElISOCBL)
process.selMuMuISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOCBL * process.skimMuMuISOCBL)

process.selElMuCONVCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuCONVCBL * process.skimElMuCONVCBL)
process.selMuElCONVCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelCONVCBL * process.skimMuElCONVCBL)
process.selElElCONVCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelCONVCBL * process.skimElElCONVCBL)
process.selMuMuCONVCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuCONVCBL * process.skimMuMuCONVCBL)

process.selElMuIPCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIPCBL * process.skimElMuIPCBL)
process.selMuElIPCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIPCBL * process.skimMuElIPCBL)
process.selElElIPCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIPCBL * process.skimElElIPCBL)
process.selMuMuIPCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIPCBL * process.skimMuMuIPCBL)

# CBT
process.selElMuIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDCBT * process.skimElMuIDCBT)
process.selMuElIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDCBT * process.skimMuElIDCBT)
process.selElElIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDCBT * process.skimElElIDCBT)
process.selMuMuIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDCBT * process.skimMuMuIDCBT)

process.selElMuISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOCBT * process.skimElMuISOCBT)
process.selMuElISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOCBT * process.skimMuElISOCBT)
process.selElElISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOCBT * process.skimElElISOCBT)
process.selMuMuISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOCBT * process.skimMuMuISOCBT)

process.selElMuCONVCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuCONVCBT * process.skimElMuCONVCBT)
process.selMuElCONVCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelCONVCBT * process.skimMuElCONVCBT)
process.selElElCONVCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelCONVCBT * process.skimElElCONVCBT)
process.selMuMuCONVCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuCONVCBT * process.skimMuMuCONVCBT)

process.selElMuIPCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIPCBT * process.skimElMuIPCBT)
process.selMuElIPCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIPCBT * process.skimMuElIPCBT)
process.selElElIPCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIPCBT * process.skimElElIPCBT)
process.selMuMuIPCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIPCBT * process.skimMuMuIPCBT)


process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    # 0
    process.selElMu0     ,  process.selMuEl0     ,  process.selElEl0     ,  process.selMuMu0      ,
    # LHL
    process.selElMuIDLHL ,  process.selMuElIDLHL ,  process.selElElIDLHL ,  process.selMuMuIDLHL  , 
    process.selElMuISOLHL,  process.selMuElISOLHL,  process.selElElISOLHL,  process.selMuMuISOLHL , 
    process.selElMuCONVLHL, process.selMuElCONVLHL, process.selElElCONVLHL, process.selMuMuCONVLHL, 
    process.selElMuIPLHL,   process.selMuElIPLHL,   process.selElElIPLHL,   process.selMuMuIPLHL  , 
    # LHT
    process.selElMuIDLHT ,  process.selMuElIDLHT ,  process.selElElIDLHT ,  process.selMuMuIDLHT  , 
    process.selElMuISOLHT,  process.selMuElISOLHT,  process.selElElISOLHT,  process.selMuMuISOLHT , 
    process.selElMuCONVLHT, process.selMuElCONVLHT, process.selElElCONVLHT, process.selMuMuCONVLHT, 
    process.selElMuIPLHT,   process.selMuElIPLHT,   process.selElElIPLHT,   process.selMuMuIPLHT  , 
    # CBL
    process.selElMuIDCBL ,  process.selMuElIDCBL ,  process.selElElIDCBL ,  process.selMuMuIDCBL  , 
    process.selElMuISOCBL,  process.selMuElISOCBL,  process.selElElISOCBL,  process.selMuMuISOCBL , 
    process.selElMuCONVCBL, process.selMuElCONVCBL, process.selElElCONVCBL, process.selMuMuCONVCBL, 
    process.selElMuIPCBL,   process.selMuElIPCBL,   process.selElElIPCBL,   process.selMuMuIPCBL  , 
    # CBT
    process.selElMuIDCBT ,  process.selMuElIDCBT ,  process.selElElIDCBT ,  process.selMuMuIDCBT  , 
    process.selElMuISOCBT,  process.selMuElISOCBT,  process.selElElISOCBT,  process.selMuMuISOCBT , 
    process.selElMuCONVCBT, process.selMuElCONVCBT, process.selElElCONVCBT, process.selMuMuCONVCBT, 
    process.selElMuIPCBT,   process.selMuElIPCBT,   process.selElElIPCBT,   process.selMuMuIPCBT  , 
    # end
    process.e
)




