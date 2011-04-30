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

process.GlobalTag.globaltag = 'RMMEGlobalTag'
# process.GlobalTag.globaltag = 'START311_V2::All'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:RMMEFN'),
#     fileNames = cms.untracked.vstring('file:../SkimStep/latinosYieldSkim.root'),
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

process.wwelmuJetNoPU = process.wwelmuIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmuelJetNoPU = process.wwmuelIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwelelJetNoPU = process.wwelelIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmumuJetNoPU = process.wwmumuIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )


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


# LHT
process.skimElMuIDLHT = process.skimElMu0.clone( src = "wwelmuIDLHT" )
process.skimMuElIDLHT = process.skimElMu0.clone( src = "wwmuelIDLHT" )
process.skimElElIDLHT = process.skimElMu0.clone( src = "wwelelIDLHT" )
process.skimMuMuIDLHT = process.skimElMu0.clone( src = "wwmumuIDLHT" )

process.skimElMuISOLHT = process.skimElMu0.clone( src = "wwelmuISOLHT" )
process.skimMuElISOLHT = process.skimElMu0.clone( src = "wwmuelISOLHT" )
process.skimElElISOLHT = process.skimElMu0.clone( src = "wwelelISOLHT" )
process.skimMuMuISOLHT = process.skimElMu0.clone( src = "wwmumuISOLHT" )

process.skimElMuCONVLHT = process.skimElMu0.clone( src = "wwelmuCONVLHT" )
process.skimMuElCONVLHT = process.skimElMu0.clone( src = "wwmuelCONVLHT" )
process.skimElElCONVLHT = process.skimElMu0.clone( src = "wwelelCONVLHT" )
process.skimMuMuCONVLHT = process.skimElMu0.clone( src = "wwmumuCONVLHT" )

process.skimElMuIPLHT = process.skimElMu0.clone( src = "wwelmuIPLHT" )
process.skimMuElIPLHT = process.skimElMu0.clone( src = "wwmuelIPLHT" )
process.skimElElIPLHT = process.skimElMu0.clone( src = "wwelelIPLHT" )
process.skimMuMuIPLHT = process.skimElMu0.clone( src = "wwmumuIPLHT" )

process.skimElMuJetNoPU = process.skimElMu0.clone( src = "wwelmuJetNoPU" )
process.skimMuElJetNoPU = process.skimElMu0.clone( src = "wwmuelJetNoPU" )
process.skimElElJetNoPU = process.skimElMu0.clone( src = "wwelelJetNoPU" )
process.skimMuMuJetNoPU = process.skimElMu0.clone( src = "wwmumuJetNoPU" )

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('RMMEFN'),
#     fileName = cms.untracked.string('hypoEvents.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep edmTriggerResults_TriggerResults_*_*',
        'keep *_pfMet_*_*',
        'keep *_tcMet_*_*',
        'keep *_slimPatJetsTriggerMatch_*_*',
        'keep *_slimPatJetsTriggerMatchNoPU_*_*',
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

process.selElMuJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuJetNoPU * process.skimElMuJetNoPU)
process.selMuElJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelJetNoPU * process.skimMuElJetNoPU)
process.selElElJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelJetNoPU * process.skimElElJetNoPU)
process.selMuMuJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuJetNoPU * process.skimMuMuJetNoPU)



process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    # 0
    process.selElMu0     ,  process.selMuEl0     ,  process.selElEl0     ,  process.selMuMu0      ,
    # LHT
    process.selElMuIDLHT ,  process.selMuElIDLHT ,  process.selElElIDLHT ,  process.selMuMuIDLHT  , 
    process.selElMuISOLHT,  process.selMuElISOLHT,  process.selElElISOLHT,  process.selMuMuISOLHT , 
    process.selElMuCONVLHT, process.selMuElCONVLHT, process.selElElCONVLHT, process.selMuMuCONVLHT, 
    process.selElMuIPLHT,   process.selMuElIPLHT,   process.selElElIPLHT,   process.selMuMuIPLHT  , 
    process.selElMuJetNoPU,   process.selMuElJetNoPU,   process.selElElJetNoPU,   process.selMuMuJetNoPU  , 
    # end
    process.e
)




