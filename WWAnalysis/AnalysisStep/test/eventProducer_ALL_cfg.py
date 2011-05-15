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

isMC = RMMEMC
# isMC = True
# isMC = False

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
#     fileNames = cms.untracked.vstring('file:/home/mwlebour/higgs/hww/code/newPFIdSelector.414/src/WWAnalysis/SkimStep/latinosYieldSkim.root'),
#     fileNames = cms.untracked.vstring('file:../SkimStep/latinosYieldSkim.root'),
#     fileNames = cms.untracked.vstring('file:/data/mwlebour/MC/Spring11/WW_414_SKIM_V04/ggToH150toWWto2L2Nu.root'),
    inputCommands = cms.untracked.vstring( "keep *" )
)
##############################################

#from glob import glob
# process.source.fileNames += [ 'file:%s'%x for x in glob('/nfs/bluearc/group/skims/ww/nov12Skim/WWTo2L2NuZ26noPU/*.root') ]
# process.source.fileNames += [ 'file:%s'%x for x in glob('/home/mangano/skim/CMSSW_4_1_3/src/workingDir/latinosYieldSkim.root') ]
# process.source.fileNames = ['file:hwwSkim.root']

# import commands
# myDir = "/castor/cern.ch/user/m/mwlebour/hww/trees/WW_414_SKIM_V00/mwlebour/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/WW_414_SKIM_V00_ID026_WJetsToLNuMad_Xfer/2bd6c1b7661d2d45dc2e9f49eb14f680/"
# process.source.fileNames = [ 'rfio:%s'%myDir+x for x in commands.getoutput("rfdir "+myDir+" | awk '{print $9}'").split() ] 


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")

process.onlyHiggsGen = cms.EDProducer( "GenParticlePruner",
    src = cms.InputTag("prunedGen"),
    select = cms.vstring(
        "drop  *  ",
        "keep pdgId =   {h0}",
    )
)

if isMC:
    process.genPath = cms.Path(process.onlyHiggsGen)
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","REDIGI311X")
else:
    process.genPath = cms.Path()
    process.skimEventProducer.triggerTag = cms.InputTag("TriggerResults","","HLT")



# 0
process.wwelmu0 = process.skimEventProducer.clone(hypoType = "WWELMU", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmuel0 = process.skimEventProducer.clone(hypoType = "WWMUEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwelel0 = process.skimEventProducer.clone(hypoType = "WWELEL", muTag="wwMuMatch", elTag = "wwEleMatch")
process.wwmumu0 = process.skimEventProducer.clone(hypoType = "WWMUMU", muTag="wwMuMatch", elTag = "wwEleMatch")

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

# PF LHT
process.wwelmuISOPFLHT  = process.wwelmu0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleISOPFLHT" )
process.wwmuelISOPFLHT  = process.wwmuel0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleISOPFLHT" )
process.wwelelISOPFLHT  = process.wwelel0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleISOPFLHT" )
process.wwmumuISOPFLHT  = process.wwmumu0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleISOPFLHT" )

process.wwelmuCONVPFLHT = process.wwelmu0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleCONVPFLHT" )
process.wwmuelCONVPFLHT = process.wwmuel0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleCONVPFLHT" )
process.wwelelCONVPFLHT = process.wwelel0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleCONVPFLHT" )
process.wwmumuCONVPFLHT = process.wwmumu0.clone(muTag = "wwMuonsISOPF", elTag = "wwEleCONVPFLHT" )

process.wwelmuIPPFLHT   = process.wwelmu0.clone(muTag = "wwMuonsIPPF",  elTag = "wwEleIPPFLHT" )
process.wwmuelIPPFLHT   = process.wwmuel0.clone(muTag = "wwMuonsIPPF",  elTag = "wwEleIPPFLHT" )
process.wwelelIPPFLHT   = process.wwelel0.clone(muTag = "wwMuonsIPPF",  elTag = "wwEleIPPFLHT" )
process.wwmumuIPPFLHT   = process.wwmumu0.clone(muTag = "wwMuonsIPPF",  elTag = "wwEleIPPFLHT" )

process.wwelmuJetNoPU = process.wwelmuIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmuelJetNoPU = process.wwmuelIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwelelJetNoPU = process.wwelelIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmumuJetNoPU = process.wwmumuIPLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                   tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )

process.wwelmuPFJetNoPU = process.wwelmuIPPFLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                       tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmuelPFJetNoPU = process.wwmuelIPPFLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                       tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwelelPFJetNoPU = process.wwelelIPPFLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
                                                       tagJetTag = cms.InputTag("slimPatJetsTriggerMatchNoPU") )
process.wwmumuPFJetNoPU = process.wwmumuIPPFLHT.clone( jetTag    = cms.InputTag("slimPatJetsTriggerMatchNoPU"), 
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

process.skimElMuISOPFLHT = process.skimElMu0.clone( src = "wwelmuISOPFLHT" )
process.skimMuElISOPFLHT = process.skimElMu0.clone( src = "wwmuelISOPFLHT" )
process.skimElElISOPFLHT = process.skimElMu0.clone( src = "wwelelISOPFLHT" )
process.skimMuMuISOPFLHT = process.skimElMu0.clone( src = "wwmumuISOPFLHT" )

process.skimElMuCONVPFLHT = process.skimElMu0.clone( src = "wwelmuCONVPFLHT" )
process.skimMuElCONVPFLHT = process.skimElMu0.clone( src = "wwmuelCONVPFLHT" )
process.skimElElCONVPFLHT = process.skimElMu0.clone( src = "wwelelCONVPFLHT" )
process.skimMuMuCONVPFLHT = process.skimElMu0.clone( src = "wwmumuCONVPFLHT" )

process.skimElMuIPPFLHT = process.skimElMu0.clone( src = "wwelmuIPPFLHT" )
process.skimMuElIPPFLHT = process.skimElMu0.clone( src = "wwmuelIPPFLHT" )
process.skimElElIPPFLHT = process.skimElMu0.clone( src = "wwelelIPPFLHT" )
process.skimMuMuIPPFLHT = process.skimElMu0.clone( src = "wwmumuIPPFLHT" )

process.skimElMuPFJetNoPU = process.skimElMu0.clone( src = "wwelmuPFJetNoPU" )
process.skimMuElPFJetNoPU = process.skimElMu0.clone( src = "wwmuelPFJetNoPU" )
process.skimElElPFJetNoPU = process.skimElMu0.clone( src = "wwelelPFJetNoPU" )
process.skimMuMuPFJetNoPU = process.skimElMu0.clone( src = "wwmumuPFJetNoPU" )

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
        'keep *_onlyHiggsGen_*_*',
        'keep GenEventInfoProduct_generator__HLT',
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
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring( 'selElMu*','selMuEl*','selElEl*','selMuMu*' )),
    # RMME Do we want to add this?
    # dropMetaData = cms.untracked.string('ALL'),
)

# RMME
# do this here or in the step 2 files?
# if isSignal:
#     process.higgsPt = cms.EDProducer("HWWKFactorProducer",
#         inputFilename = cms.untracked.string("WWAnalysis/Misc/Scales/scalefactor.hmRMMEMASS.dat"),
#         inputFilename = cms.untracked.string("WWAnalysis/Misc/Scales/scalefactor.mh160.dat"),
#         ProcessID = cms.untracked.int32(10010),
#         Debug =cms.untracked.bool(False)
#     )
#     process.eventHists.ptWeight = cms.InputTag("higgsPt")
#     process.p = cms.Path(process.higgsPt * process.eventHists)
# else:
#     process.p = cms.Path(process.eventHists)


# 0
process.selElMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 * process.skimElMu0)
process.selMuEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 * process.skimMuEl0)
process.selElEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 * process.skimElEl0)
process.selMuMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 * process.skimMuMu0)


# LHT

process.selElMuISOPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOPFLHT * process.skimElMuISOPFLHT)
process.selMuElISOPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOPFLHT * process.skimMuElISOPFLHT)
process.selElElISOPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOPFLHT * process.skimElElISOPFLHT)
process.selMuMuISOPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOPFLHT * process.skimMuMuISOPFLHT)

process.selElMuCONVPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuCONVPFLHT * process.skimElMuCONVPFLHT)
process.selMuElCONVPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelCONVPFLHT * process.skimMuElCONVPFLHT)
process.selElElCONVPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelCONVPFLHT * process.skimElElCONVPFLHT)
process.selMuMuCONVPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuCONVPFLHT * process.skimMuMuCONVPFLHT)

process.selElMuIPPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIPPFLHT * process.skimElMuIPPFLHT)
process.selMuElIPPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIPPFLHT * process.skimMuElIPPFLHT)
process.selElElIPPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIPPFLHT * process.skimElElIPPFLHT)
process.selMuMuIPPFLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIPPFLHT * process.skimMuMuIPPFLHT)

process.selElMuPFJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuPFJetNoPU * process.skimElMuPFJetNoPU)
process.selMuElPFJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelPFJetNoPU * process.skimMuElPFJetNoPU)
process.selElElPFJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelPFJetNoPU * process.skimElElPFJetNoPU)
process.selMuMuPFJetNoPU = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuPFJetNoPU * process.skimMuMuPFJetNoPU)

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
    process.genPath,
    # 0
    process.selElMu0     ,  process.selMuEl0     ,  process.selElEl0     ,  process.selMuMu0      ,
    # LHT
    process.selElMuIDLHT ,  process.selMuElIDLHT ,  process.selElElIDLHT ,  process.selMuMuIDLHT  , 
    process.selElMuISOLHT,  process.selMuElISOLHT,  process.selElElISOLHT,  process.selMuMuISOLHT , 
    process.selElMuCONVLHT, process.selMuElCONVLHT, process.selElElCONVLHT, process.selMuMuCONVLHT, 
    process.selElMuIPLHT,   process.selMuElIPLHT,   process.selElElIPLHT,   process.selMuMuIPLHT  , 
    #process.selElMuJetNoPU, process.selMuElJetNoPU, process.selElElJetNoPU, process.selMuMuJetNoPU  , 
    # PFLHT
    #process.selElMuISOPFLHT,  process.selMuElISOPFLHT,  process.selElElISOPFLHT,  process.selMuMuISOPFLHT , 
    #process.selElMuCONVPFLHT, process.selMuElCONVPFLHT, process.selElElCONVPFLHT, process.selMuMuCONVPFLHT, 
    #process.selElMuIPPFLHT,   process.selMuElIPPFLHT,   process.selElElIPPFLHT,   process.selMuMuIPPFLHT  , 
    #process.selElMuPFJetNoPU, process.selMuElPFJetNoPU, process.selElElPFJetNoPU, process.selMuMuPFJetNoPU  , 
    # end
    process.e
)




