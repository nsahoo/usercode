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

process.GlobalTag.globaltag = 'START311_V2::All'
# process.GlobalTag.globaltag = 'SET_GLOBALTAG'


process.es_prefer_mag = cms.ESPrefer("AutoMagneticFieldESProducer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### HERE I SET THE SAMPLE I WANT TO RUN ON ###
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('file:id101160.1.V01.root'),
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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("WWAnalysis.AnalysisStep.skimEventProducer_cfi")
process.load("WWAnalysis.AnalysisStep.wwElectrons_cfi")
process.load("WWAnalysis.AnalysisStep.wwMuons_cfi")


### HERE I SET THE SKIMEVENT PRODUCER CONFIGURATION ###
process.wwelmu0 = process.skimEventProducer.clone(hypoType = cms.string("WWELMU"))
process.wwmuel0 = process.skimEventProducer.clone(hypoType = cms.string("WWMUEL"))
process.wwelel0 = process.skimEventProducer.clone(hypoType = cms.string("WWELEL"))
process.wwmumu0 = process.skimEventProducer.clone(hypoType = cms.string("WWMUMU"))
process.wwelmu0.muTag = cms.InputTag("wwMuMatch")
process.wwmuel0.muTag = cms.InputTag("wwMuMatch")
process.wwelel0.muTag = cms.InputTag("wwMuMatch")
process.wwmumu0.muTag = cms.InputTag("wwMuMatch")
process.wwelmu0.elTag = cms.InputTag("wwEleMatch")
process.wwmuel0.elTag = cms.InputTag("wwEleMatch")
process.wwelel0.elTag = cms.InputTag("wwEleMatch")
process.wwmumu0.elTag = cms.InputTag("wwEleMatch")

process.wwelmuIDLHL = process.wwelmu0.clone()
process.wwmuelIDLHL = process.wwmuel0.clone()
process.wwelelIDLHL = process.wwelel0.clone()
process.wwmumuIDLHL = process.wwmumu0.clone()
process.wwelmuIDLHL.muTag = cms.InputTag("wwMuons")
process.wwmuelIDLHL.muTag = cms.InputTag("wwMuons")
process.wwelelIDLHL.muTag = cms.InputTag("wwMuons")
process.wwmumuIDLHL.muTag = cms.InputTag("wwMuons")
process.wwelmuIDLHL.elTag = cms.InputTag("wwEleIDLHL")
process.wwmuelIDLHL.elTag = cms.InputTag("wwEleIDLHL")
process.wwelelIDLHL.elTag = cms.InputTag("wwEleIDLHL")
process.wwmumuIDLHL.elTag = cms.InputTag("wwEleIDLHL")

process.wwelmuISOLHL = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleISOLHL"))
process.wwmuelISOLHL = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleISOLHL"))
process.wwelelISOLHL = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleISOLHL"))
process.wwmumuISOLHL = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleISOLHL"))

process.wwelmuIDLHT = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleIDLHT"))
process.wwmuelIDLHT = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleIDLHT"))
process.wwelelIDLHT = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleIDLHT"))
process.wwmumuIDLHT = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleIDLHT"))

process.wwelmuISOLHT = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleISOLHT"))
process.wwmuelISOLHT = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleISOLHT"))
process.wwelelISOLHT = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleISOLHT"))
process.wwmumuISOLHT = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleISOLHT"))

process.wwelmuIDCBL = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleIDCBL"))
process.wwmuelIDCBL = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleIDCBL"))
process.wwelelIDCBL = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleIDCBL"))
process.wwmumuIDCBL = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleIDCBL"))

process.wwelmuISOCBL = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleISOCBL"))
process.wwmuelISOCBL = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleISOCBL"))
process.wwelelISOCBL = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleISOCBL"))
process.wwmumuISOCBL = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleISOCBL"))

process.wwelmuIDCBT = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleIDCBT"))
process.wwmuelIDCBT = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleIDCBT"))
process.wwelelIDCBT = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleIDCBT"))
process.wwmumuIDCBT = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleIDCBT"))

process.wwelmuISOCBT = process.wwelmuIDLHL.clone(elTag = cms.InputTag("wwEleISOCBT"))
process.wwmuelISOCBT = process.wwmuelIDLHL.clone(elTag = cms.InputTag("wwEleISOCBT"))
process.wwelelISOCBT = process.wwelelIDLHL.clone(elTag = cms.InputTag("wwEleISOCBT"))
process.wwmumuISOCBT = process.wwmumuIDLHL.clone(elTag = cms.InputTag("wwEleISOCBT"))

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

process.skimElMuIDLHL = process.skimElMu0.clone( src = "wwelmuIDLHL" )
process.skimMuElIDLHL = process.skimElMu0.clone( src = "wwmuelIDLHL" )
process.skimElElIDLHL = process.skimElMu0.clone( src = "wwelelIDLHL" )
process.skimMuMuIDLHL = process.skimElMu0.clone( src = "wwmumuIDLHL" )

process.skimElMuISOLHL = process.skimElMu0.clone( src = "wwelmuISOLHL" )
process.skimMuElISOLHL = process.skimElMu0.clone( src = "wwmuelISOLHL" )
process.skimElElISOLHL = process.skimElMu0.clone( src = "wwelelISOLHL" )
process.skimMuMuISOLHL = process.skimElMu0.clone( src = "wwmumuISOLHL" )

process.skimElMuIDLHT = process.skimElMu0.clone( src = "wwelmuIDLHT" )
process.skimMuElIDLHT = process.skimElMu0.clone( src = "wwmuelIDLHT" )
process.skimElElIDLHT = process.skimElMu0.clone( src = "wwelelIDLHT" )
process.skimMuMuIDLHT = process.skimElMu0.clone( src = "wwmumuIDLHT" )

process.skimElMuISOLHT = process.skimElMu0.clone( src = "wwelmuISOLHT" )
process.skimMuElISOLHT = process.skimElMu0.clone( src = "wwmuelISOLHT" )
process.skimElElISOLHT = process.skimElMu0.clone( src = "wwelelISOLHT" )
process.skimMuMuISOLHT = process.skimElMu0.clone( src = "wwmumuISOLHT" )

process.skimElMuIDCBL = process.skimElMu0.clone( src = "wwelmuIDCBL" )
process.skimMuElIDCBL = process.skimElMu0.clone( src = "wwmuelIDCBL" )
process.skimElElIDCBL = process.skimElMu0.clone( src = "wwelelIDCBL" )
process.skimMuMuIDCBL = process.skimElMu0.clone( src = "wwmumuIDCBL" )

process.skimElMuISOCBL = process.skimElMu0.clone( src = "wwelmuISOCBL" )
process.skimMuElISOCBL = process.skimElMu0.clone( src = "wwmuelISOCBL" )
process.skimElElISOCBL = process.skimElMu0.clone( src = "wwelelISOCBL" )
process.skimMuMuISOCBL = process.skimElMu0.clone( src = "wwmumuISOCBL" )

process.skimElMuIDCBT = process.skimElMu0.clone( src = "wwelmuIDCBT" )
process.skimMuElIDCBT = process.skimElMu0.clone( src = "wwmuelIDCBT" )
process.skimElElIDCBT = process.skimElMu0.clone( src = "wwelelIDCBT" )
process.skimMuMuIDCBT = process.skimElMu0.clone( src = "wwmumuIDCBT" )

process.skimElMuISOCBT = process.skimElMu0.clone( src = "wwelmuISOCBT" )
process.skimMuElISOCBT = process.skimElMu0.clone( src = "wwmuelISOCBT" )
process.skimElElISOCBT = process.skimElMu0.clone( src = "wwelelISOCBT" )
process.skimMuMuISOCBT = process.skimElMu0.clone( src = "wwmumuISOCBT" )


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('skimEvent.allHWW160.root'),
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

process.selElMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmu0 * process.skimElMu0)
process.selMuEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuel0 * process.skimMuEl0)
process.selElEl0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelel0 * process.skimElEl0)
process.selMuMu0 = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumu0 * process.skimMuMu0)

process.selElMuIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDLHL * process.skimElMuIDLHL)
process.selMuElIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDLHL * process.skimMuElIDLHL)
process.selElElIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDLHL * process.skimElElIDLHL)
process.selMuMuIDLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDLHL * process.skimMuMuIDLHL)

process.selElMuISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOLHL * process.skimElMuISOLHL)
process.selMuElISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOLHL * process.skimMuElISOLHL)
process.selElElISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOLHL * process.skimElElISOLHL)
process.selMuMuISOLHL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOLHL * process.skimMuMuISOLHL)

process.selElMuIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDLHT * process.skimElMuIDLHT)
process.selMuElIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDLHT * process.skimMuElIDLHT)
process.selElElIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDLHT * process.skimElElIDLHT)
process.selMuMuIDLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDLHT * process.skimMuMuIDLHT)

process.selElMuISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOLHT * process.skimElMuISOLHT)
process.selMuElISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOLHT * process.skimMuElISOLHT)
process.selElElISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOLHT * process.skimElElISOLHT)
process.selMuMuISOLHT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOLHT * process.skimMuMuISOLHT)

process.selElMuIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDCBL * process.skimElMuIDCBL)
process.selMuElIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDCBL * process.skimMuElIDCBL)
process.selElElIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDCBL * process.skimElElIDCBL)
process.selMuMuIDCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDCBL * process.skimMuMuIDCBL)

process.selElMuISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOCBL * process.skimElMuISOCBL)
process.selMuElISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOCBL * process.skimMuElISOCBL)
process.selElElISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOCBL * process.skimElElISOCBL)
process.selMuMuISOCBL = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOCBL * process.skimMuMuISOCBL)

process.selElMuIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuIDCBT * process.skimElMuIDCBT)
process.selMuElIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelIDCBT * process.skimMuElIDCBT)
process.selElElIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelIDCBT * process.skimElElIDCBT)
process.selMuMuIDCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuIDCBT * process.skimMuMuIDCBT)

process.selElMuISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelmuISOCBT * process.skimElMuISOCBT)
process.selMuElISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmuelISOCBT * process.skimMuElISOCBT)
process.selElElISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwelelISOCBT * process.skimElElISOCBT)
process.selMuMuISOCBT = cms.Path(process.wwElectronSequence + process.wwMuonSequence * process.wwmumuISOCBT * process.skimMuMuISOCBT)


process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(
    process.selElMu0     , process.selMuEl0     , process.selElEl0     , process.selMuMu0     ,
    process.selElMuIDLHL , process.selMuElIDLHL , process.selElElIDLHL , process.selMuMuIDLHL , 
    process.selElMuISOLHL, process.selMuElISOLHL, process.selElElISOLHL, process.selMuMuISOLHL, 
    process.selElMuIDLHT , process.selMuElIDLHT , process.selElElIDLHT , process.selMuMuIDLHT , 
    process.selElMuISOLHT, process.selMuElISOLHT, process.selElElISOLHT, process.selMuMuISOLHT, 
    process.selElMuIDCBL , process.selMuElIDCBL , process.selElElIDCBL , process.selMuMuIDCBL , 
    process.selElMuISOCBL, process.selMuElISOCBL, process.selElElISOCBL, process.selMuMuISOCBL, 
    process.selElMuIDCBT , process.selMuElIDCBT , process.selElElIDCBT , process.selMuMuIDCBT , 
    process.selElMuISOCBT, process.selMuElISOCBT, process.selElElISOCBT, process.selMuMuISOCBT,
    process.e
)




